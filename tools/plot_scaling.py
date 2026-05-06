#!/usr/bin/env python3
"""tools/plot_scaling.py — strong + weak scaling plots from scaling.cpp CSV.

Consumes one or more CSV files emitted by `benchmarks/md/scaling.cpp`
(columns: ranks,atoms,steps,time_total_s,time_per_step_ms,
atoms_per_rank,mode) and produces two PNG figures:

  * scaling_strong.png — log-log time/step vs ranks, with ideal-scaling
    dashed line; strong-scaling efficiency annotated at 8 / 32 ranks.
  * scaling_weak.png   — efficiency vs ranks at fixed atoms/rank; ideal
    line at y=1.

Strong-scaling rows are CSV rows with `mode == "strong"` (or detected
heuristically as "ranks varies, atoms is constant"). Weak-scaling rows
have `mode == "weak"` (or "ranks varies, atoms_per_rank constant").

Pure matplotlib (no seaborn). Refuses to plot if no rows of the
respective mode are found, but still emits the other plot and reports.

Usage:
    python tools/plot_scaling.py runs/scaling_*.csv \\
        --out-dir benchmarks/md/reference/polyethylene/

CSV may have a header row (the script tolerates either).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")  # no display required
import matplotlib.pyplot as plt  # noqa: E402


# CSV columns emitted by scaling.cpp (see file header docstring).
COLUMNS = (
    "ranks",
    "atoms",
    "steps",
    "time_total_s",
    "time_per_step_ms",
    "atoms_per_rank",
    "mode",
)


def _read_rows(paths: Iterable[Path]) -> list[dict]:
    rows: list[dict] = []
    for p in paths:
        with p.open() as fh:
            sample = fh.read(4096)
            fh.seek(0)
            has_header = "ranks" in sample.split("\n", 1)[0]
            reader = (
                csv.DictReader(fh)
                if has_header
                else csv.DictReader(fh, fieldnames=list(COLUMNS))
            )
            for raw in reader:
                if raw is None:
                    continue
                # Coerce numerics; tolerate stray whitespace.
                try:
                    rows.append({
                        "ranks":            int(float(raw["ranks"])),
                        "atoms":            int(float(raw["atoms"])),
                        "steps":            int(float(raw["steps"])),
                        "time_total_s":     float(raw["time_total_s"]),
                        "time_per_step_ms": float(raw["time_per_step_ms"]),
                        "atoms_per_rank":   int(float(raw["atoms_per_rank"])),
                        "mode":             (raw.get("mode") or "").strip().lower(),
                    })
                except (KeyError, TypeError, ValueError) as exc:
                    print(f"  [warn] skipping bad row in {p}: {exc}",
                          file=sys.stderr)
    return rows


def _classify(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split rows into (strong, weak) by mode column or heuristic fallback.

    Heuristic fallback: among rows of unset mode, if `atoms` is constant
    across distinct rank counts, classify as strong; if `atoms_per_rank`
    is constant, classify as weak; otherwise drop.
    """
    strong = [r for r in rows if r["mode"] == "strong"]
    weak   = [r for r in rows if r["mode"] == "weak"]

    untagged = [r for r in rows if r["mode"] not in ("strong", "weak")]
    if untagged:
        atoms_set = {r["atoms"] for r in untagged}
        apr_set   = {r["atoms_per_rank"] for r in untagged}
        if len(atoms_set) == 1 and len(apr_set) > 1:
            strong.extend(untagged)
        elif len(apr_set) == 1 and len(atoms_set) > 1:
            weak.extend(untagged)
        else:
            print(f"  [warn] {len(untagged)} rows have no `mode` and ambiguous "
                  "atoms/atoms_per_rank — skipped", file=sys.stderr)

    strong.sort(key=lambda r: r["ranks"])
    weak.sort(key=lambda r: r["ranks"])
    return strong, weak


def plot_strong(rows: list[dict], out_path: Path) -> None:
    if not rows:
        print(f"  [skip] no strong-scaling rows; not writing {out_path.name}")
        return
    ranks = [r["ranks"] for r in rows]
    times = [r["time_per_step_ms"] for r in rows]
    if ranks[0] != 1:
        print(f"  [warn] strong-scaling does not start at rank 1; ideal line "
              f"will be relative to rank={ranks[0]}", file=sys.stderr)
    t1 = times[0]
    ideal = [t1 / (n / ranks[0]) for n in ranks]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(ranks, times, "o-", label="measured", linewidth=1.6, markersize=7)
    ax.loglog(ranks, ideal, "k--", label="ideal", linewidth=1.0)
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("time / MD step (ms)")
    ax.set_title(f"Strong scaling — atoms={rows[0]['atoms']}")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)

    # Annotate efficiency at 8 / 32 ranks if present.
    for marker in (8, 32):
        for r in rows:
            if r["ranks"] == marker:
                eff = (t1 / r["time_per_step_ms"]) / (marker / ranks[0])
                ax.annotate(f"{eff*100:.0f}% eff @ {marker}r",
                            xy=(marker, r["time_per_step_ms"]),
                            xytext=(8, 8), textcoords="offset points",
                            fontsize=9)
                break
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  [ok]   wrote {out_path}")


def plot_weak(rows: list[dict], out_path: Path) -> None:
    if not rows:
        print(f"  [skip] no weak-scaling rows; not writing {out_path.name}")
        return
    ranks = [r["ranks"] for r in rows]
    times = [r["time_per_step_ms"] for r in rows]
    t1 = times[0]
    eff = [t1 / t for t in times]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(ranks, eff, "o-", label="measured", linewidth=1.6, markersize=7)
    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.0, label="ideal")
    ax.set_xlabel("MPI ranks")
    ax.set_ylabel("weak-scaling efficiency  (t1 / t_n)")
    ax.set_ylim(0, 1.2)
    ax.set_title(f"Weak scaling — atoms/rank={rows[0]['atoms_per_rank']}")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  [ok]   wrote {out_path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", type=Path,
                    help="One or more scaling.cpp CSV files.")
    ap.add_argument("--out-dir", type=Path, default=Path("."),
                    help="Output directory for PNGs (default: cwd).")
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = _read_rows(args.csv)
    if not rows:
        print("error: no rows parsed from input CSV(s).", file=sys.stderr)
        return 1
    strong, weak = _classify(rows)
    print(f"  parsed {len(rows)} rows; {len(strong)} strong + {len(weak)} weak")

    plot_strong(strong, args.out_dir / "scaling_strong.png")
    plot_weak(weak,     args.out_dir / "scaling_weak.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
