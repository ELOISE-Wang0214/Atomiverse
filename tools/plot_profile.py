#!/usr/bin/env python3
"""tools/plot_profile.py — stacked-bar plot of per-component MD profile.

Consumes a CSV emitted by `benchmarks/md/pe_profile_single.cpp`
(columns: component,total_ms,pct_of_step,n_calls,mean_us) and produces
a single PNG sibling at `<input.csv>.png` (or `--out`) showing each
component's share of total MD-step wall time as a stacked horizontal
bar. Components below `--threshold` percent of step time collapse into
"other" so the chart stays legible.

Usage:
    python tools/plot_profile.py benchmarks/md/reference/polyethylene/profile_8.csv
    python tools/plot_profile.py profile_1.csv --out my_plot.png --threshold 1.0
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _read_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            try:
                rows.append({
                    "component":   raw["component"].strip(),
                    "total_ms":    float(raw["total_ms"]),
                    "pct_of_step": float(raw["pct_of_step"]),
                    "n_calls":     int(float(raw["n_calls"])),
                    "mean_us":     float(raw["mean_us"]),
                })
            except (KeyError, TypeError, ValueError) as exc:
                print(f"  [warn] skipping bad row in {path}: {exc}",
                      file=sys.stderr)
    return rows


def _collapse(rows: list[dict], threshold_pct: float) -> list[dict]:
    """Roll up rows below `threshold_pct` into one 'other' bucket."""
    big = [r for r in rows if r["pct_of_step"] >= threshold_pct]
    small = [r for r in rows if r["pct_of_step"] < threshold_pct]
    if small:
        big.append({
            "component":   f"other ({len(small)} comps < {threshold_pct:g}%)",
            "total_ms":    sum(r["total_ms"] for r in small),
            "pct_of_step": sum(r["pct_of_step"] for r in small),
            "n_calls":     0,
            "mean_us":     0.0,
        })
    big.sort(key=lambda r: r["pct_of_step"], reverse=True)
    return big


def plot(rows: list[dict], out_path: Path, source_name: str) -> None:
    if not rows:
        print(f"  [skip] no rows; not writing {out_path.name}")
        return

    labels = [r["component"] for r in rows]
    pcts   = [r["pct_of_step"] for r in rows]

    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.barh(range(len(labels)), pcts,
                    color=plt.cm.tab20(range(len(labels))))
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()  # largest at top
    ax.set_xlabel("% of MD-step wall time")
    ax.set_title(f"Per-component breakdown — {source_name}")
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)

    for bar, pct, row in zip(bars, pcts, rows):
        if pct < 0.5:
            continue
        ax.annotate(f"{pct:.1f}% ({row['mean_us']:.1f} µs/call)",
                    xy=(bar.get_width(), bar.get_y() + bar.get_height() / 2),
                    xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  [ok]   wrote {out_path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", type=Path, help="profile_*.csv from pe_profile_single.")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output PNG path (default: <input>.png).")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="Collapse components below this %% (default: 0.5).")
    args = ap.parse_args(argv)

    if not args.csv.exists():
        print(f"error: {args.csv} not found", file=sys.stderr)
        return 1

    out_path = args.out
    if out_path is None:
        out_path = args.csv.with_suffix(".png")

    rows = _read_rows(args.csv)
    if not rows:
        print("error: no rows parsed.", file=sys.stderr)
        return 1
    rows = _collapse(rows, args.threshold)

    plot(rows, out_path, args.csv.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
