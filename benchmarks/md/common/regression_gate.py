#!/usr/bin/env python3
"""Performance-regression gate for Atomiverse benchmark runs.

Implements ``benchmarks/REGRESSION.md``: joins a *candidate* perf run against a
*baseline* run on the same ``(host, backend, engine, system, mode, ng)`` key and
applies a two-band verdict (WARN ~3%, FAIL ~8%) with a geomean + >=2-key
aggregate guard, so a single noisy sweep point cannot fail the whole job.

It consumes the committed per-run CSVs under

    benchmarks/<domain>/<case>/performance/runs/<sha>_<date>/results_<backend>.csv

whose real schema is ``engine,system,mode,ng,atoms,matom_step_s,vram_gb,rep``,
plus the paired meta JSON for the host/backend guard (``meta_cpu.json`` for the
CPU CSV, ``meta.json`` for the GPU CSV). Repeats (the ``rep`` column) are reduced
to a per-key median; the gate never re-times or re-runs anything.

Rules faithful to REGRESSION.md:
  * Only ``engine == AV`` rows are gated -- LAMMPS/CP2K reference rows are an
    external comparison, never ours to regress (§2).
  * Two runs on different hosts are never compared; the gate skips (exit 0) when
    the candidate and baseline hosts are both known and differ (§2).
  * Metric polarity is normalised by orientation so ``ratio > 1`` is always
    "worse": throughput (``matom_step_s``) is higher-is-better; latency
    (``s_per_iter``, ``wall_s``, ``real_time``) and memory (``vram_gb``) are
    lower-is-better (§1, §4).
  * Aggregate guard (§4): a run FAILs only when the geometric-mean ratio exceeds
    the fail band **or** >= 2 keys individually exceed it; a lone point over the
    fail band is a WARN.

The google-benchmark JSON targets (``bench_md_step``, ``bench_native``) are out
of scope here -- use upstream ``benchmark/tools/compare.py`` per §5.

Stdlib only (argparse, csv, json, statistics); no pandas.

Examples
--------
    regression_gate.py --candidate runs/<sha>_<date>  --baseline runs/<sha0>_<d0>
    regression_gate.py --candidate cand/results_cpu.csv \
                       --baseline base/results_cpu.csv --metric matom_step_s
    regression_gate.py --candidate c.csv --baseline b.csv \
                       --memory-metric vram_gb --warn 0.03 --fail 0.08
    regression_gate.py --selftest        # run the inline helper unit tests
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import statistics
import sys
from typing import Dict, List, NamedTuple, Optional, Tuple

# --- metric polarity ---------------------------------------------------------
# "worse = larger" is the gate's canonical orientation. For higher-is-better
# metrics the worse-ratio is baseline/candidate; for lower-is-better it is
# candidate/baseline (see worse_ratio()).
HIGHER_IS_BETTER = frozenset({"matom_step_s", "items_per_second"})
LOWER_IS_BETTER = frozenset({"s_per_iter", "wall_s", "real_time", "cpu_time", "vram_gb"})

ENGINE_GATED = "AV"  # only our own rows are gated (§2)

# Verdict ordering for "worst of" folding.
_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}


# --- pure, unit-testable helpers --------------------------------------------
def median(vals: List[float]) -> float:
    """Median of a non-empty list of floats (reduces the `rep` repeats)."""
    if not vals:
        raise ValueError("median() of empty sequence")
    return statistics.median(vals)


def geomean(vals: List[float]) -> float:
    """Geometric mean of positive ratios; 1.0 for an empty list.

    Uses log-space accumulation to stay stable across a wide sweep.
    """
    positive = [v for v in vals if v > 0.0]
    if not positive:
        return 1.0
    return math.exp(sum(math.log(v) for v in positive) / len(positive))


def metric_is_higher_better(metric: str, override: Optional[bool] = None) -> bool:
    """Resolve a metric's polarity: explicit override wins, else the table.

    Raises for an unknown metric with no override -- no silent fallback.
    """
    if override is not None:
        return override
    if metric in HIGHER_IS_BETTER:
        return True
    if metric in LOWER_IS_BETTER:
        return False
    raise ValueError(
        f"unknown metric {metric!r}; pass --higher-is-better or --lower-is-better"
    )


def worse_ratio(candidate: float, baseline: float, higher_is_better: bool) -> float:
    """Ratio oriented so >1 always means the candidate got *worse* (§4).

    higher-is-better -> baseline/candidate (throughput dropping is worse);
    lower-is-better  -> candidate/baseline (latency/memory rising is worse).
    """
    if candidate <= 0.0 or baseline <= 0.0:
        raise ValueError("non-positive metric value cannot be gated")
    return baseline / candidate if higher_is_better else candidate / baseline


def classify_ratio(ratio: float, warn_band: float, fail_band: float) -> str:
    """Per-key band: PASS (r<=warn), WARN (warn<r<=fail), FAIL (r>fail)."""
    if ratio <= warn_band:
        return "PASS"
    if ratio <= fail_band:
        return "WARN"
    return "FAIL"


def aggregate_verdict(
    ratios: List[float],
    warn_band: float,
    fail_band: float,
    min_fail_keys: int = 2,
) -> str:
    """Fold per-key ratios into a run verdict with the noise guard (§4).

    FAIL iff geomean(ratios) > fail_band OR >= min_fail_keys keys exceed
    fail_band. Otherwise WARN if any key is above warn_band or geomean exceeds
    it. Otherwise PASS. A lone point over the fail band is only a WARN.
    """
    if not ratios:
        return "PASS"
    gm = geomean(ratios)
    n_fail = sum(1 for r in ratios if r > fail_band)
    if gm > fail_band or n_fail >= min_fail_keys:
        return "FAIL"
    if gm > warn_band or any(r > warn_band for r in ratios):
        return "WARN"
    return "PASS"


def backend_from_csv_name(name: str) -> Optional[str]:
    """'results_cpu.csv'->'cpu', 'results_gpu.csv'->'gpu' (None if unknown)."""
    base = os.path.basename(name)
    if "_cpu" in base:
        return "cpu"
    if "_gpu" in base:
        return "gpu"
    return None


def band_ratio(fraction: float) -> float:
    """Convert a fractional tolerance (0.08) to a worse-ratio band (1.08)."""
    return 1.0 + fraction


# --- data model --------------------------------------------------------------
class Key(NamedTuple):
    """Join key (§2). engine is pinned to AV; host is a run-level guard."""

    backend: str
    engine: str
    system: str
    mode: str
    ng: int


class Cell(NamedTuple):
    """Per-key reduced measurement."""

    metric: float
    memory: Optional[float]
    atoms: Optional[int]


class Run(NamedTuple):
    """A loaded perf run: its cells plus provenance for the host guard."""

    csv_path: str
    backend: str
    host: Optional[str]
    cells: Dict[Key, Cell]


# --- loading -----------------------------------------------------------------
def resolve_csv(path: str, backend: str) -> str:
    """Turn a run dir or CSV path into a concrete results CSV path."""
    if os.path.isfile(path):
        return path
    if os.path.isdir(path):
        if backend != "auto":
            cand = os.path.join(path, f"results_{backend}.csv")
            if os.path.isfile(cand):
                return cand
            raise FileNotFoundError(f"no results_{backend}.csv under {path}")
        found = sorted(
            p for p in glob.glob(os.path.join(path, "results_*.csv"))
            if "_min" not in os.path.basename(p)
        )
        if len(found) == 1:
            return found[0]
        if not found:
            raise FileNotFoundError(f"no results_*.csv under {path}")
        raise ValueError(
            f"multiple result CSVs under {path}; pass --backend cpu|gpu: "
            + ", ".join(os.path.basename(p) for p in found)
        )
    raise FileNotFoundError(path)


def resolve_meta(csv_path: str, backend: str) -> Dict[str, object]:
    """Load the sibling meta JSON (meta_cpu.json for cpu, meta.json else)."""
    d = os.path.dirname(os.path.abspath(csv_path))
    name = "meta_cpu.json" if backend == "cpu" else "meta.json"
    meta_path = os.path.join(d, name)
    if os.path.isfile(meta_path):
        try:
            with open(meta_path) as fh:
                return json.load(fh)
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def load_run(
    path: str,
    metric: str,
    memory_metric: Optional[str],
    backend_arg: str,
    host_override: Optional[str],
) -> Run:
    """Load a run (dir or CSV) into reduced per-key cells (§5 steps 1-3, 5)."""
    csv_path = resolve_csv(path, backend_arg)
    backend = (
        backend_arg
        if backend_arg != "auto"
        else (backend_from_csv_name(csv_path) or "unknown")
    )
    meta = resolve_meta(csv_path, backend)
    host = host_override or (str(meta["host"]) if "host" in meta else None)
    meta_backend = meta.get("backend")
    if isinstance(meta_backend, str) and backend_arg == "auto":
        # meta_cpu.json carries an explicit "cpu"; trust it over the filename.
        backend = meta_backend

    # Group raw rows by key, collect the per-rep values, then reduce to median.
    raw: Dict[Key, Dict[str, List[float]]] = {}
    atoms_of: Dict[Key, Optional[int]] = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("engine") != ENGINE_GATED:
                continue  # LAMMPS/CP2K reference rows are not gated (§2, step 3)
            try:
                key = Key(
                    backend=backend,
                    engine=ENGINE_GATED,
                    system=row["system"],
                    mode=row["mode"],
                    ng=int(row["ng"]),
                )
                metric_val = float(row[metric])
            except (KeyError, ValueError):
                continue
            bucket = raw.setdefault(key, {"metric": [], "memory": []})
            bucket["metric"].append(metric_val)
            if memory_metric and memory_metric in row:
                try:
                    bucket["memory"].append(float(row[memory_metric]))
                except ValueError:
                    pass
            try:
                atoms_of[key] = int(row["atoms"])
            except (KeyError, ValueError):
                atoms_of.setdefault(key, None)

    cells: Dict[Key, Cell] = {}
    for key, bucket in raw.items():
        if not bucket["metric"]:
            continue
        mem = median(bucket["memory"]) if bucket["memory"] else None
        cells[key] = Cell(
            metric=median(bucket["metric"]),
            memory=mem,
            atoms=atoms_of.get(key),
        )
    return Run(csv_path=csv_path, backend=backend, host=host, cells=cells)


# --- reporting ---------------------------------------------------------------
class Comparison(NamedTuple):
    key: Key
    baseline: float
    candidate: float
    ratio: float
    band: str


def compare_runs(
    candidate: Run,
    baseline: Run,
    higher_is_better: bool,
    warn_band: float,
    fail_band: float,
) -> Tuple[List[Comparison], List[Key], List[Key]]:
    """Inner-join and classify (§5 steps 4-6). Returns (rows, only_cand, only_base)."""
    shared = sorted(set(candidate.cells) & set(baseline.cells))
    only_cand = sorted(set(candidate.cells) - set(baseline.cells))
    only_base = sorted(set(baseline.cells) - set(candidate.cells))
    rows: List[Comparison] = []
    for key in shared:
        cand_v = candidate.cells[key].metric
        base_v = baseline.cells[key].metric
        try:
            ratio = worse_ratio(cand_v, base_v, higher_is_better)
        except ValueError:
            continue
        rows.append(
            Comparison(key, base_v, cand_v, ratio, classify_ratio(ratio, warn_band, fail_band))
        )
    return rows, only_cand, only_base


def _fmt_key(key: Key) -> str:
    return f"{key.backend}/{key.system}/{key.mode}/ng={key.ng}"


def print_table(
    title: str,
    rows: List[Comparison],
    higher_is_better: bool,
) -> None:
    """Human-readable per-key table (§5 step 8)."""
    print(f"\n{title}")
    if not rows:
        print("  (no shared keys)")
        return
    kw = max(len(_fmt_key(r.key)) for r in rows)
    kw = max(kw, len("key"))
    hdr = f"  {'key':<{kw}}  {'baseline':>12}  {'candidate':>12}  {'ratio':>7}  {'delta%':>7}  band"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        # delta% in the metric's natural direction (positive = improvement).
        if higher_is_better:
            delta = (r.candidate - r.baseline) / r.baseline * 100.0
        else:
            delta = (r.baseline - r.candidate) / r.baseline * 100.0
        print(
            f"  {_fmt_key(r.key):<{kw}}  {r.baseline:>12.4g}  {r.candidate:>12.4g}"
            f"  {r.ratio:>7.4f}  {delta:>+7.2f}  {r.band}"
        )


def run_gate(args: argparse.Namespace) -> int:
    """Execute the gate; return the process exit code (0 unless FAIL)."""
    higher = metric_is_higher_better(
        args.metric,
        True if args.higher_is_better else (False if args.lower_is_better else None),
    )
    warn_band = band_ratio(args.warn)
    fail_band = band_ratio(args.fail)

    candidate = load_run(args.candidate, args.metric, args.memory_metric, args.backend, args.host)
    baseline = load_run(args.baseline, args.metric, args.memory_metric, args.backend, args.host)

    print("Atomiverse performance-regression gate")
    print(f"  metric        : {args.metric} ({'higher' if higher else 'lower'} is better)")
    print(f"  candidate     : {candidate.csv_path}")
    print(f"                  host={candidate.host or '?'} backend={candidate.backend}")
    print(f"  baseline      : {baseline.csv_path}")
    print(f"                  host={baseline.host or '?'} backend={baseline.backend}")
    print(f"  bands         : warn>{args.warn:.0%}  fail>{args.fail:.0%}"
          f"  (aggregate: geomean or >={args.min_fail_keys} keys over fail)")

    # Host guard (§2): never compare across hardware.
    if candidate.host and baseline.host and candidate.host != baseline.host:
        print(
            f"\nSKIP: candidate host {candidate.host!r} != baseline host "
            f"{baseline.host!r} -- perf numbers are not comparable across "
            "hardware; no same-host baseline. (exit 0)"
        )
        return 0

    rows, only_cand, only_base = compare_runs(candidate, baseline, higher, warn_band, fail_band)
    print_table(f"Primary metric: {args.metric}", rows, higher)

    if only_base:
        print("\n  keys in baseline only (removed / not re-run): "
              + ", ".join(_fmt_key(k) for k in only_base))
    if only_cand:
        print("  keys in candidate only (new): "
              + ", ".join(_fmt_key(k) for k in only_cand))

    ratios = [r.ratio for r in rows]
    primary_verdict = aggregate_verdict(ratios, warn_band, fail_band, args.min_fail_keys)
    gm = geomean(ratios)
    n_fail = sum(1 for r in ratios if r > fail_band)
    n_warn = sum(1 for r in rows if r.band == "WARN")
    print(f"\n  primary: {len(rows)} keys  geomean={gm:.4f}"
          f"  warn={n_warn}  fail={n_fail}  -> {primary_verdict}")

    overall = primary_verdict

    # Secondary memory gate (§4) with its own looser band and guard.
    if args.memory_metric:
        mrows, _, _ = _memory_rows(candidate, baseline, band_ratio(args.warn),
                                   band_ratio(args.memory_fail))
        print_table(f"Memory metric: {args.memory_metric} (looser fail>{args.memory_fail:.0%})",
                    mrows, higher_is_better=False)
        m_ratios = [r.ratio for r in mrows]
        mem_verdict = aggregate_verdict(
            m_ratios, band_ratio(args.warn), band_ratio(args.memory_fail), args.min_fail_keys
        )
        print(f"\n  memory: {len(mrows)} keys  geomean={geomean(m_ratios):.4f} -> {mem_verdict}")
        if _ORDER[mem_verdict] > _ORDER[overall]:
            overall = mem_verdict

    print(f"\nVERDICT: {overall}")
    if overall == "FAIL" and args.warn_only:
        print("  (--warn-only: not failing the job; band-calibration mode, §7)")
        return 0
    return 1 if overall == "FAIL" else 0


def _memory_rows(
    candidate: Run,
    baseline: Run,
    warn_band: float,
    fail_band: float,
) -> Tuple[List[Comparison], List[Key], List[Key]]:
    """Memory keys are lower-is-better; reuse the join with the memory cell."""
    shared = sorted(set(candidate.cells) & set(baseline.cells))
    rows: List[Comparison] = []
    for key in shared:
        cm = candidate.cells[key].memory
        bm = baseline.cells[key].memory
        if cm is None or bm is None:
            continue
        try:
            ratio = worse_ratio(cm, bm, higher_is_better=False)
        except ValueError:
            continue
        rows.append(Comparison(key, bm, cm, ratio, classify_ratio(ratio, warn_band, fail_band)))
    return rows, [], []


# --- inline unit tests -------------------------------------------------------
def _selftest() -> int:
    """Assert the pure helpers; run with --selftest (no external deps)."""
    assert median([3.0, 1.0, 2.0]) == 2.0
    assert median([2.0, 4.0]) == 3.0
    assert abs(geomean([1.0, 4.0]) - 2.0) < 1e-12
    assert geomean([]) == 1.0
    # throughput dropped 10% -> worse-ratio ~1.111 > 1
    r = worse_ratio(candidate=90.0, baseline=100.0, higher_is_better=True)
    assert abs(r - 100.0 / 90.0) < 1e-12 and r > 1.0
    # latency rose 10% -> worse-ratio 1.10
    assert abs(worse_ratio(110.0, 100.0, higher_is_better=False) - 1.10) < 1e-12
    # identical run -> ratio 1.0 -> PASS
    assert worse_ratio(100.0, 100.0, True) == 1.0
    assert classify_ratio(1.00, 1.03, 1.08) == "PASS"
    assert classify_ratio(1.05, 1.03, 1.08) == "WARN"
    assert classify_ratio(1.20, 1.03, 1.08) == "FAIL"
    # one lone point over fail -> WARN (noise guard); two -> FAIL
    assert aggregate_verdict([1.00, 1.00, 1.20], 1.03, 1.08) == "WARN"
    assert aggregate_verdict([1.20, 1.20, 1.00], 1.03, 1.08) == "FAIL"
    # geomean over fail even with <2 individual fails -> FAIL
    assert aggregate_verdict([1.09, 1.09], 1.03, 1.08) == "FAIL"
    assert aggregate_verdict([1.00, 1.01], 1.03, 1.08) == "PASS"
    assert metric_is_higher_better("matom_step_s") is True
    assert metric_is_higher_better("s_per_iter") is False
    assert metric_is_higher_better("anything", override=True) is True
    assert backend_from_csv_name("x/results_gpu.csv") == "gpu"
    assert backend_from_csv_name("results_cpu_min.csv") == "cpu"
    assert abs(band_ratio(0.08) - 1.08) < 1e-12
    print("selftest OK")
    return 0


# --- CLI ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="regression_gate.py",
        description="Performance-regression gate for Atomiverse perf runs "
                    "(benchmarks/REGRESSION.md). Compares a candidate run "
                    "against a baseline on the same host/backend/case key.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Exit status: 0 on PASS/WARN/skip, 1 on FAIL.",
    )
    p.add_argument("--candidate", "--current", dest="candidate",
                   help="candidate run: a run dir or a results_*.csv path")
    p.add_argument("--baseline",
                   help="baseline run: a run dir or a results_*.csv path")
    p.add_argument("--backend", choices=["cpu", "gpu", "auto"], default="auto",
                   help="which results_<backend>.csv to read (default: auto-detect)")
    p.add_argument("--metric", default="matom_step_s",
                   help="metric column to gate (default: matom_step_s)")
    pol = p.add_mutually_exclusive_group()
    pol.add_argument("--higher-is-better", action="store_true",
                     help="override polarity: larger metric is better (throughput)")
    pol.add_argument("--lower-is-better", action="store_true",
                     help="override polarity: smaller metric is better (latency)")
    p.add_argument("--warn", type=float, default=0.03,
                   help="warn band as a fraction (default: 0.03 = 3%%)")
    p.add_argument("--fail", type=float, default=0.08,
                   help="fail band as a fraction (default: 0.08 = 8%%)")
    p.add_argument("--memory-metric", default=None,
                   help="optional secondary lower-is-better metric, e.g. vram_gb")
    p.add_argument("--memory-fail", type=float, default=0.15,
                   help="looser fail band for --memory-metric (default: 0.15)")
    p.add_argument("--min-fail-keys", type=int, default=2,
                   help="min individual FAIL keys for a run FAIL (default: 2, §4)")
    p.add_argument("--host", default=None,
                   help="override host (when a run has no meta JSON)")
    p.add_argument("--warn-only", action="store_true",
                   help="never exit non-zero; band-calibration mode (§7)")
    p.add_argument("--selftest", action="store_true",
                   help="run the inline helper unit tests and exit")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.selftest:
        return _selftest()
    if not args.candidate or not args.baseline:
        build_parser().error("--candidate and --baseline are required "
                             "(or pass --selftest)")
    return run_gate(args)


if __name__ == "__main__":
    sys.exit(main())
