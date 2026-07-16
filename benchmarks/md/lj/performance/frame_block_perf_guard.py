#!/usr/bin/env python3
"""Freeze and compare the two single-rank CPU LJ performance sentinels.

This is the narrow performance guard for the Frame Block cutover.  It drives
the existing ``scaling`` benchmark in both neighbor-reuse and forced-
reneighbor modes, reduces independent process runs to separate medians, and
stores enough provenance to reject comparisons made on a different host or
with a different build configuration.

It deliberately does not cover materialization or MPI communication.  Those
paths have different scaling variables and must be measured separately.

Typical use::

    frame_block_perf_guard.py --binary build/benchmarks/scaling \
        --output /tmp/frame-block-before.json
    frame_block_perf_guard.py --binary build/benchmarks/scaling \
        --baseline /tmp/frame-block-before.json \
        --output /tmp/frame-block-after.json

Each sentinel fails independently when median step latency regresses by more
than 2%; there is no aggregate or geomean escape hatch.  A case whose relative
median absolute deviation exceeds 1% makes the whole run INCONCLUSIVE rather
than a false PASS/FAIL.  The output JSON is written even for an inconclusive or
failed run so the evidence is inspectable.  Schema v1 snapshots are rejected:
this is a hard-cut acceptance format, not a compatibility reader.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


SCHEMA = "atomiverse.frame-block-cpu-lj.v2"
BENCHMARK = "cpu-lj-neighbor-md-two-sentinel-single-rank"
CSV_RE = re.compile(
    r"^(?P<ranks>\d+),(?P<atoms>\d+),(?P<steps>\d+),"
    r"(?P<total>[0-9.eE+-]+),(?P<step>[0-9.eE+-]+),"
    r"(?P<per_rank>\d+),(?P<mode>[^,]+),(?P<rss>\d+)$"
)
PROVENANCE_RE = re.compile(
    r"^# atomiverse-perf xsimd=(?P<xsimd>[0-9.]+) "
    r"xtensor=(?P<xtensor>[0-9.]+) "
    r"storage_batch=(?P<storage_batch>\d+) work_batch=(?P<work_batch>\d+)$"
)
NEIGHBOR_RE = re.compile(
    r"^# atomiverse-neighbor policy=(?P<policy>adaptive|forced) "
    r"rebuilds=(?P<rebuilds>\d+) launches=(?P<launches>\d+)$"
)

CASE_NEIGHBOR_REUSE = "neighbor-reuse"
CASE_FORCED_RENEIGHBOR = "forced-reneighbor"

RELEVANT_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_CXX_COMPILER",
    "CMAKE_CXX_FLAGS",
    "CMAKE_CXX_FLAGS_RELEASE",
    "ATV_BUILD_CUDA",
    "ATV_REAL_STORAGE",
    "ATV_REAL_WORK",
    "ATV_ID_SIZE",
    "ATV_USE_MPI",
)


class GuardError(RuntimeError):
    """An invalid benchmark run or incomparable snapshot."""


def _run_text(command: Sequence[str], cwd: Optional[Path] = None) -> str:
    try:
        proc = subprocess.run(
            list(command),
            cwd=str(cwd) if cwd else None,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        output = getattr(exc, "stdout", "") or ""
        raise GuardError(f"command failed: {' '.join(command)}\n{output}") from exc
    return proc.stdout


def parse_scaling_row(output: str) -> Dict[str, Any]:
    """Return the final machine-readable row from ``scaling`` output."""
    matches = []
    for line in output.splitlines():
        match = CSV_RE.fullmatch(line.strip())
        if match:
            matches.append(match.groupdict())
    if not matches:
        raise GuardError("scaling produced no machine-readable CSV row")

    raw = matches[-1]
    row: Dict[str, Any] = {
        "ranks": int(raw["ranks"]),
        "atoms": int(raw["atoms"]),
        "steps": int(raw["steps"]),
        "time_total_s": float(raw["total"]),
        "time_per_step_ms": float(raw["step"]),
        "atoms_per_rank": int(raw["per_rank"]),
        "mode": raw["mode"],
        "peak_rss_total_kb": int(raw["rss"]),
    }
    if row["ranks"] != 1:
        raise GuardError(f"Frame Block sentinel requires one rank, got {row['ranks']}")
    if row["mode"] != "strong":
        raise GuardError(f"expected strong-scaling fixture, got {row['mode']!r}")
    if row["time_per_step_ms"] <= 0.0:
        raise GuardError("non-positive step latency")
    return row


def parse_provenance(output: str) -> Dict[str, Any]:
    """Return the compile-time dependency/SIMD record emitted by ``scaling``."""
    matches = []
    for line in output.splitlines():
        match = PROVENANCE_RE.fullmatch(line.strip())
        if match:
            matches.append(match.groupdict())
    if len(matches) != 1:
        raise GuardError(
            "scaling must emit exactly one --provenance record; "
            f"observed {len(matches)}"
        )
    raw = matches[0]
    return {
        "xsimd": raw["xsimd"],
        "xtensor": raw["xtensor"],
        "storage_batch": int(raw["storage_batch"]),
        "work_batch": int(raw["work_batch"]),
    }


def parse_neighbor_stats(output: str) -> Dict[str, Any]:
    """Return the measured neighbor-list policy and rebuild count."""
    matches = []
    for line in output.splitlines():
        match = NEIGHBOR_RE.fullmatch(line.strip())
        if match:
            matches.append(match.groupdict())
    if len(matches) != 1:
        raise GuardError(
            "scaling must emit exactly one neighbor stats record; "
            f"observed {len(matches)}"
        )
    raw = matches[0]
    return {
        "policy": raw["policy"],
        "rebuilds": int(raw["rebuilds"]),
        "launches": int(raw["launches"]),
    }


def median_and_rmad(values: Sequence[float]) -> Tuple[float, float]:
    """Median and relative median absolute deviation (MAD / median)."""
    if not values:
        raise GuardError("no timing samples")
    median = float(statistics.median(values))
    mad = float(statistics.median(abs(value - median) for value in values))
    return median, mad / median


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _find_cmake_cache(binary: Path) -> Optional[Path]:
    for parent in (binary.parent, *binary.parents):
        candidate = parent / "CMakeCache.txt"
        if candidate.is_file():
            return candidate
    return None


def _read_cmake_cache(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    result: Dict[str, str] = {}
    wanted = set(RELEVANT_CACHE_KEYS)
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        typed_key, value = line.split("=", 1)
        key = typed_key.split(":", 1)[0]
        if key in wanted:
            result[key] = value
    return result


def _cpu_model() -> str:
    if sys.platform == "darwin":
        for key in ("machdep.cpu.brand_string", "hw.model"):
            try:
                value = _run_text(("sysctl", "-n", key)).strip()
            except GuardError:
                continue
            if value:
                return value
    if sys.platform.startswith("linux"):
        try:
            for line in Path("/proc/cpuinfo").read_text(errors="replace").splitlines():
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
    return platform.processor() or "unknown"


def _git_provenance(repo: Path) -> Dict[str, Any]:
    try:
        sha = _run_text(("git", "rev-parse", "HEAD"), repo).strip()
        status = _run_text(
            ("git", "status", "--porcelain", "--untracked-files=no"), repo
        ).strip()
    except GuardError:
        return {"sha": None, "tracked_dirty": None}
    return {"sha": sha, "tracked_dirty": bool(status)}


def _compiler_version(cache: Mapping[str, str]) -> Optional[str]:
    compiler = cache.get("CMAKE_CXX_COMPILER")
    if not compiler:
        return None
    try:
        return _run_text((compiler, "--version")).splitlines()[0].strip()
    except GuardError:
        return None


def _compile_command(cache_path: Optional[Path]) -> Optional[str]:
    if cache_path is None:
        return None
    database = cache_path.parent / "compile_commands.json"
    if not database.is_file():
        return None
    try:
        entries = json.loads(database.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    for entry in entries:
        source = str(entry.get("file", ""))
        if source.endswith("benchmarks/md/lj/performance/scaling.cpp"):
            if "command" in entry:
                command = str(entry["command"])
            else:
                arguments = entry.get("arguments")
                if not isinstance(arguments, list):
                    return None
                command = " ".join(str(argument) for argument in arguments)

            # A pre-cutover build must live in a separate checkout.  Absolute
            # source/build paths are not compiler provenance and would make an
            # otherwise identical pair incomparable.  Preserve every flag and
            # dependency path shape while replacing only those two roots.
            suffix = Path("benchmarks/md/lj/performance/scaling.cpp")
            source_roots: set[Path] = set()
            source_path = Path(source)
            for candidate in (source_path.absolute(), source_path.resolve()):
                for parent in candidate.parents:
                    try:
                        if candidate.relative_to(parent) == suffix:
                            source_roots.add(parent)
                            break
                    except ValueError:
                        continue
            build_roots = {cache_path.parent.absolute(), cache_path.parent.resolve()}
            directory = entry.get("directory")
            if directory:
                directory_path = Path(str(directory))
                build_roots.update(
                    {directory_path.absolute(), directory_path.resolve()}
                )
            # Build roots are nested under source roots. Replace them first so
            # `/checkout/build` cannot collapse to `<source>/build` merely
            # because `/tmp` and `/private/tmp` name the same macOS directory.
            replacements = [
                *((root, "<build>") for root in build_roots),
                *((root, "<source>") for root in source_roots),
            ]
            for root, marker in replacements:
                command = command.replace(str(root), marker)
            return command
    return None


def environment_fingerprint(
    cache: Mapping[str, str],
    cache_path: Optional[Path],
    benchmark_provenance: Mapping[str, Any],
) -> Dict[str, Any]:
    """Fields that must match before a 2% comparison is meaningful."""
    return {
        "host": platform.node(),
        "system": platform.system(),
        "machine": platform.machine(),
        "cpu_model": _cpu_model(),
        "logical_cpus": os.cpu_count(),
        "compiler": _compiler_version(cache),
        "cmake": dict(cache),
        "compile_command": _compile_command(cache_path),
        "dependencies": dict(benchmark_provenance),
        "thread_environment": {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "close",
            "OMP_PLACES": "cores",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        },
    }


def _repo_from_binary(binary: Path) -> Path:
    cache = _find_cmake_cache(binary)
    start = cache.parent if cache else binary.parent
    try:
        top = _run_text(("git", "rev-parse", "--show-toplevel"), start).strip()
        return Path(top)
    except GuardError:
        return start


def run_fixture(args: argparse.Namespace) -> Dict[str, Any]:
    binary = Path(args.binary).expanduser().resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise GuardError(f"benchmark is not executable: {binary}")

    base_command = [
        str(binary),
        "--atoms",
        str(args.atoms),
        "--steps",
        str(args.steps),
        "--warmup",
        str(args.warmup),
        "--cutoff",
        str(args.cutoff),
        "--skin",
        str(args.skin),
        "--dt",
        str(args.dt),
        "--quiet",
        "--provenance",
    ]
    env = os.environ.copy()
    # Timing instrumentation is opt-in by environment in the benchmark.  The
    # acceptance sentinel measures the production hot path, so inherited local
    # profiling switches must not silently add per-component timers.
    env.pop("ATV_TIMER", None)
    env.pop("ATV_TIMER_VERBOSE", None)
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "OMP_PROC_BIND": "close",
            "OMP_PLACES": "cores",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )

    case_definitions = (
        (CASE_NEIGHBOR_REUSE, args.reuse_cadence, args.reuse_cadence),
        (CASE_FORCED_RENEIGHBOR, 1, 1),
    )
    rows_by_case: Dict[str, List[Dict[str, Any]]] = {
        name: [] for name, _, _ in case_definitions
    }
    benchmark_provenance: Optional[Dict[str, Any]] = None

    for repetition in range(1, args.repetitions + 1):
        ordered_cases = (
            case_definitions
            if repetition % 2 == 1
            else tuple(reversed(case_definitions))
        )
        for name, reneigh, migrate in ordered_cases:
            command = [
                *base_command,
                "--reneigh",
                str(reneigh),
                "--migrate",
                str(migrate),
            ]
            if name == CASE_FORCED_RENEIGHBOR:
                command.append("--force-reneighbor")
            try:
                proc = subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    env=env,
                )
            except (OSError, subprocess.CalledProcessError) as exc:
                output = getattr(exc, "stdout", "") or ""
                raise GuardError(
                    f"{name} repetition {repetition} failed\n{output}"
                ) from exc
            row = parse_scaling_row(proc.stdout)
            observed_provenance = parse_provenance(proc.stdout)
            neighbor_stats = parse_neighbor_stats(proc.stdout)
            expected_policy = "forced" if name == CASE_FORCED_RENEIGHBOR else "adaptive"
            expected_launches = args.warmup + args.steps
            expected_rebuilds = (
                expected_launches + 1 if name == CASE_FORCED_RENEIGHBOR else 1
            )
            if neighbor_stats != {
                "policy": expected_policy,
                "rebuilds": expected_rebuilds,
                "launches": expected_launches,
            }:
                raise GuardError(
                    f"{name} did not exercise its required neighbor path: "
                    f"{neighbor_stats}"
                )
            row["neighbor"] = neighbor_stats
            if benchmark_provenance is None:
                benchmark_provenance = observed_provenance
            elif benchmark_provenance != observed_provenance:
                raise GuardError("benchmark provenance changed between repetitions")
            if row["atoms"] != args.atoms or row["steps"] != args.steps:
                raise GuardError(f"benchmark echoed the wrong fixture: {row}")
            rows_by_case[name].append(row)
            print(
                f"{name:19s} rep {repetition:02d}/{args.repetitions}: "
                f"{row['time_per_step_ms']:.6f} ms/step",
                flush=True,
            )

    if benchmark_provenance is None:
        raise GuardError("benchmark produced no provenance")

    cases: Dict[str, Dict[str, Any]] = {}
    for name, reneigh, migrate in case_definitions:
        rows = rows_by_case[name]
        values = [float(row["time_per_step_ms"]) for row in rows]
        median, rmad = median_and_rmad(values)
        cases[name] = {
            "fixture": {
                "ranks": 1,
                "threads": 1,
                "atoms": args.atoms,
                "steps": args.steps,
                "warmup": args.warmup,
                "cutoff": args.cutoff,
                "skin": args.skin,
                "dt": args.dt,
                "reneigh": reneigh,
                "migrate": migrate,
                "force_reneighbor": name == CASE_FORCED_RENEIGHBOR,
                "repetitions": args.repetitions,
            },
            "samples": values,
            "median_ms": median,
            "relative_mad": rmad,
            "min_ms": min(values),
            "max_ms": max(values),
            "rows": rows,
        }

    cache_path = _find_cmake_cache(binary)
    cache = _read_cmake_cache(cache_path)
    repo = _repo_from_binary(binary)
    return {
        "schema": SCHEMA,
        "benchmark": BENCHMARK,
        "metric": "time_per_step_ms",
        "cases": cases,
        "environment": environment_fingerprint(cache, cache_path, benchmark_provenance),
        "provenance": {
            "binary": str(binary),
            "binary_sha256": _sha256(binary),
            "cmake_cache": str(cache_path) if cache_path else None,
            "git": _git_provenance(repo),
        },
    }


def _comparison_differences(
    baseline: Mapping[str, Any], candidate: Mapping[str, Any]
) -> List[str]:
    differences = []
    for field in ("schema", "benchmark", "metric", "environment"):
        if baseline.get(field) != candidate.get(field):
            differences.append(field)
    baseline_cases = baseline.get("cases")
    candidate_cases = candidate.get("cases")
    if not isinstance(baseline_cases, Mapping) or not isinstance(
        candidate_cases, Mapping
    ):
        differences.append("cases")
        return differences
    if set(baseline_cases) != set(candidate_cases):
        differences.append("case names")
        return differences
    for name in baseline_cases:
        baseline_case = baseline_cases[name]
        candidate_case = candidate_cases[name]
        if not isinstance(baseline_case, Mapping) or not isinstance(
            candidate_case, Mapping
        ):
            differences.append(f"case {name}")
        elif baseline_case.get("fixture") != candidate_case.get("fixture"):
            differences.append(f"case {name} fixture")
    return differences


def compare_case(
    baseline: Mapping[str, Any], candidate: Mapping[str, Any], max_regression: float
) -> Tuple[float, str]:
    """Compare one sentinel; global comparability is checked separately."""
    before = float(baseline["median_ms"])
    after = float(candidate["median_ms"])
    if before <= 0.0 or after <= 0.0:
        raise GuardError("snapshot contains non-positive median latency")
    regression = after / before - 1.0
    exceeds = regression > max_regression and not math.isclose(
        regression, max_regression, rel_tol=1.0e-12, abs_tol=1.0e-12
    )
    return regression, "FAIL" if exceeds else "PASS"


def compare_snapshots(
    baseline: Mapping[str, Any], candidate: Mapping[str, Any], max_regression: float
) -> Dict[str, Dict[str, Any]]:
    differences = _comparison_differences(baseline, candidate)
    if differences:
        raise GuardError(
            "baseline and candidate are not comparable; differing fields: "
            + ", ".join(differences)
        )
    baseline_cases = baseline["cases"]
    candidate_cases = candidate["cases"]
    comparisons: Dict[str, Dict[str, Any]] = {}
    for name in (CASE_NEIGHBOR_REUSE, CASE_FORCED_RENEIGHBOR):
        regression, verdict = compare_case(
            baseline_cases[name], candidate_cases[name], max_regression
        )
        comparisons[name] = {"regression": regression, "verdict": verdict}
    return comparisons


def _write_snapshot(path: Path, snapshot: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")


def _selftest() -> None:
    row = parse_scaling_row(
        "noise\n# atomiverse-perf xsimd=14.2.0 xtensor=0.27.1 "
        "storage_batch=2 work_batch=2\n"
        "1,32768,500,2.0,4.0,32768,strong,12345\n"
    )
    assert row["ranks"] == 1 and row["time_per_step_ms"] == 4.0
    provenance = parse_provenance(
        "# atomiverse-perf xsimd=14.2.0 xtensor=0.27.1 storage_batch=2 work_batch=2\n"
    )
    assert provenance["xsimd"] == "14.2.0" and provenance["work_batch"] == 2
    neighbor = parse_neighbor_stats(
        "# atomiverse-neighbor policy=forced rebuilds=251 launches=250\n"
    )
    assert neighbor == {"policy": "forced", "rebuilds": 251, "launches": 250}
    median, rmad = median_and_rmad([4.0, 4.1, 3.9])
    assert median == 4.0 and abs(rmad - 0.025) < 1.0e-12
    fixture = {"x": 1}
    base = {
        "schema": SCHEMA,
        "benchmark": BENCHMARK,
        "metric": "time_per_step_ms",
        "environment": {"host": "runner"},
        "cases": {
            CASE_NEIGHBOR_REUSE: {"fixture": fixture, "median_ms": 10.0},
            CASE_FORCED_RENEIGHBOR: {"fixture": fixture, "median_ms": 20.0},
        },
    }
    candidate = {
        **base,
        "cases": {
            CASE_NEIGHBOR_REUSE: {"fixture": fixture, "median_ms": 10.2},
            CASE_FORCED_RENEIGHBOR: {"fixture": fixture, "median_ms": 20.4002},
        },
    }
    comparisons = compare_snapshots(base, candidate, 0.02)
    assert comparisons[CASE_NEIGHBOR_REUSE]["verdict"] == "PASS"
    assert comparisons[CASE_FORCED_RENEIGHBOR]["verdict"] == "FAIL"
    print("frame_block_perf_guard selftest: PASS")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", help="Release scaling benchmark executable")
    parser.add_argument("--output", help="JSON snapshot to write")
    parser.add_argument("--baseline", help="prior JSON snapshot to compare")
    parser.add_argument("--atoms", type=int, default=32768)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--warmup", type=int, default=50)
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--skin", type=float, default=0.6)
    parser.add_argument("--dt", type=float, default=0.0125)
    parser.add_argument("--repetitions", type=int, default=11)
    parser.add_argument(
        "--reuse-cadence",
        type=int,
        default=1_000_000,
        help="reneighbor/migrate cadence for the reuse sentinel",
    )
    parser.add_argument(
        "--max-regression",
        type=float,
        default=0.02,
        help="maximum median latency regression (default: 0.02)",
    )
    parser.add_argument(
        "--max-relative-mad",
        type=float,
        default=0.01,
        help="maximum MAD/median before the run is inconclusive (default: 0.01)",
    )
    parser.add_argument("--selftest", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_parser().parse_args(argv)
    if args.selftest:
        _selftest()
        return 0
    if not args.binary or not args.output:
        raise GuardError("--binary and --output are required")
    if args.repetitions < 10:
        raise GuardError("at least ten independent repetitions are required")
    if min(args.atoms, args.steps, args.warmup) < 1:
        raise GuardError("atoms, steps, and warmup must be positive")
    if args.reuse_cadence <= max(args.steps, args.warmup):
        raise GuardError(
            "--reuse-cadence must exceed both --steps and --warmup so the "
            "reuse sentinel has no scheduled rebuild after step zero"
        )
    if not 0.0 <= args.max_regression < 1.0:
        raise GuardError("--max-regression must be in [0, 1)")
    if not 0.0 < args.max_relative_mad < 1.0:
        raise GuardError("--max-relative-mad must be in (0, 1)")

    snapshot = run_fixture(args)
    snapshot["verdict"] = "PENDING"
    _write_snapshot(Path(args.output).expanduser().resolve(), snapshot)
    verdict = "BASELINE"
    exit_code = 0
    noisy_candidate = [
        name
        for name, case in snapshot["cases"].items()
        if float(case["relative_mad"]) > args.max_relative_mad
    ]
    if noisy_candidate:
        verdict = "INCONCLUSIVE"
        exit_code = 2
    elif args.baseline:
        baseline_path = Path(args.baseline).expanduser().resolve()
        baseline = json.loads(baseline_path.read_text())
        if baseline.get("schema") != SCHEMA:
            raise GuardError(
                f"baseline schema must be exactly {SCHEMA!r}; "
                f"got {baseline.get('schema')!r}"
            )
        noisy_baseline = [
            name
            for name, case in baseline.get("cases", {}).items()
            if float(case.get("relative_mad", 1.0)) > args.max_relative_mad
        ]
        if noisy_baseline:
            raise GuardError(
                "baseline noise exceeds --max-relative-mad for: "
                + ", ".join(noisy_baseline)
            )
        comparisons = compare_snapshots(baseline, snapshot, args.max_regression)
        verdict = (
            "FAIL"
            if any(case["verdict"] == "FAIL" for case in comparisons.values())
            else "PASS"
        )
        snapshot["comparison"] = {
            "baseline": str(baseline_path),
            "max_regression": args.max_regression,
            "cases": comparisons,
            "verdict": verdict,
        }
        if verdict == "FAIL":
            exit_code = 1

    snapshot["verdict"] = verdict
    _write_snapshot(Path(args.output).expanduser().resolve(), snapshot)
    print(f"{verdict}: two-sentinel LJ guard")
    for name, case in snapshot["cases"].items():
        line = (
            f"  {name}: median={case['median_ms']:.6f} ms/step, "
            f"relative_MAD={case['relative_mad']:.3%}"
        )
        if "comparison" in snapshot:
            compared = snapshot["comparison"]["cases"][name]
            line += (
                f", regression={compared['regression']:.3%}, "
                f"verdict={compared['verdict']}"
            )
        elif name in noisy_candidate:
            line += ", verdict=INCONCLUSIVE"
        print(line)
    return exit_code


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (GuardError, OSError, json.JSONDecodeError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
