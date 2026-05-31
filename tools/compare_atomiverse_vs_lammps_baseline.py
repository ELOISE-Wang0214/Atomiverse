#!/usr/bin/env python3
"""Atomiverse-CPU vs LAMMPS-CPU per-stage equivalence report for the MSU PE
baseline scan.

Reads paired runs at matching MPI rank counts from:
  - <ATV_RUNS>/n<R>r/thermo.csv        (raw Atomiverse 17-col extended schema)
  - <LMP_RUNS>/n<R>r/thermo.csv        (raw LAMMPS per-stage extended schema)
and, when present, the LAMMPS stage_summary.csv (per-stage wall time / ns/day).

For each stage in {1,2,3,4,5} the driver computes the Flyvbjerg-Petersen
block-average mean ± 1σ for:
  T (K), P (bar), Pxx/Pyy/Pzz (bar), E_pot (eV), box_x/y/z (Å)
on both engines (after unit-converting LAMMPS `real` units to the Atomiverse
canonical eV/bar/Å), then renders one markdown table per stage.

Stress-strain (stage 5) is overlaid via tools/plot_msu_compare.py downstream;
this driver only emits the numeric comparison.

Outputs:
  - <OUT_DIR>/comparison.md          (the full markdown report)
  - <OUT_DIR>/per_stage_deltas.csv   (machine-readable delta table)

Pure stdlib (csv + argparse + math + statistics). Idempotent — re-running
overwrites the two output files and nothing else.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Dict, List, Tuple

# ── unit conversions (LAMMPS `real` → Atomiverse canonical) ───────────────
ATM_TO_BAR = 1.01325            # 1 atm = 1.01325 bar
KCAL_MOL_TO_EV = 0.04336411     # 1 kcal/mol = 0.04336411 eV

STAGES = [1, 2, 3, 4, 5]

# Columns of interest (canonical/Atomiverse-side names + units)
METRICS = [
    ("T", "T_K", "K"),
    ("P", "P_bar", "bar"),
    ("Pxx", "Pxx", "bar"),
    ("Pyy", "Pyy", "bar"),
    ("Pzz", "Pzz", "bar"),
    ("E_pot", "E_pot_eV", "eV"),
    ("box_x", "box_x", "Å"),
    ("box_y", "box_y", "Å"),
    ("box_z", "box_z", "Å"),
]


# ── readers ───────────────────────────────────────────────────────────────


def read_atomiverse_thermo(path: str) -> Dict[int, Dict[str, List[float]]]:
    """Group Atomiverse thermo.csv by stage column.

    Atomiverse schema (msu-pe-cuda-components-01 §G.2/G.3):
      step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,
      E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,
      E_kin_eV,box_x,box_y,box_z
    """
    by_stage: Dict[int, Dict[str, List[float]]] = {s: {} for s in STAGES}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"empty atomiverse thermo: {path}")

    # Identify the per-axis stress column key — Atomiverse writes "Pxx".
    for s in STAGES:
        for _, atv_col, _ in METRICS:
            by_stage[s][atv_col] = []

    for row in rows:
        try:
            stage = int(row.get("stage", "0"))
        except ValueError:
            continue
        if stage not in STAGES:
            continue
        for _, atv_col, _ in METRICS:
            v = row.get(atv_col)
            if v in (None, ""):
                continue
            try:
                by_stage[stage][atv_col].append(float(v))
            except ValueError:
                pass
    return by_stage


def read_lammps_thermo(path: str) -> Dict[int, Dict[str, List[float]]]:
    """Group LAMMPS molpy-reader thermo.csv by stage_index.

    Schema (from `tools/lammps_baseline_log_to_csv.py` union output):
      stage_index,Step,Temp,Press,E_pair,E_bond,E_angle,E_dihed,
      Pxx,Pyy,Pzz,Lx,Ly,Lz

    Convert in-place to Atomiverse canonical units (bar / eV / Å).
    """
    by_stage: Dict[int, Dict[str, List[float]]] = {
        s: {atv_col: [] for _, atv_col, _ in METRICS} for s in STAGES
    }
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"empty lammps thermo: {path}")

    for row in rows:
        try:
            stage = int(float(row.get("stage_index", "0")))
        except ValueError:
            continue
        if stage not in STAGES:
            continue
        # T (K) — same units
        if (v := row.get("Temp")) not in (None, ""):
            by_stage[stage]["T_K"].append(float(v))
        # P (atm → bar)
        if (v := row.get("Press")) not in (None, ""):
            by_stage[stage]["P_bar"].append(float(v) * ATM_TO_BAR)
        # Pxx/Pyy/Pzz (atm → bar)
        for col in ("Pxx", "Pyy", "Pzz"):
            if (v := row.get(col)) not in (None, ""):
                by_stage[stage][col].append(float(v) * ATM_TO_BAR)
        # E_pot = epair + ebond + eangle + edihed (kcal/mol → eV)
        terms = []
        for c in ("E_pair", "E_bond", "E_angle", "E_dihed"):
            v = row.get(c)
            if v in (None, ""):
                terms = []
                break
            try:
                terms.append(float(v))
            except ValueError:
                terms = []
                break
        if terms:
            by_stage[stage]["E_pot_eV"].append(sum(terms) * KCAL_MOL_TO_EV)
        # box (Lx/Ly/Lz → box_x/y/z, Å in both)
        for src, dst in (("Lx", "box_x"), ("Ly", "box_y"), ("Lz", "box_z")):
            if (v := row.get(src)) not in (None, ""):
                by_stage[stage][dst].append(float(v))
    return by_stage


def read_lammps_stage_summary(path: str) -> Dict[int, Dict[str, float]]:
    out: Dict[int, Dict[str, float]] = {}
    if not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                s = int(float(row.get("stage_index", "0")))
            except ValueError:
                continue
            out[s] = {k: row.get(k, "") for k in row}
    return out


# ── Flyvbjerg-Petersen block-average σ ────────────────────────────────────


def block_average(series: List[float]) -> Tuple[float, float, int]:
    """Mean and 1σ via Flyvbjerg & Petersen 1989 block-renormalization.

    Doubles block size until the per-block σ plateaus, then returns the
    converged σ. Falls back to plain σ/√N when the series is too short
    (< 4 blocks at any size).
    """
    n = len(series)
    if n == 0:
        return float("nan"), float("nan"), 0
    if n < 4:
        m = sum(series) / n
        if n < 2:
            return m, 0.0, n
        var = sum((x - m) ** 2 for x in series) / (n - 1)
        return m, math.sqrt(var / n), n

    mean = sum(series) / n
    # Renormalize block size until under 4 blocks remain.
    cur = list(series)
    last_sigma_block = 0.0
    while len(cur) >= 4:
        nb = len(cur)
        m_b = sum(cur) / nb
        var_b = sum((x - m_b) ** 2 for x in cur) / (nb - 1)
        sigma_block = math.sqrt(var_b / nb)
        last_sigma_block = sigma_block
        cur = [(cur[2 * i] + cur[2 * i + 1]) * 0.5 for i in range(nb // 2)]
    return mean, last_sigma_block, n


# ── reporting ─────────────────────────────────────────────────────────────


def _fmt(v: float) -> str:
    if not math.isfinite(v):
        return "—"
    a = abs(v)
    if a >= 100 or (a < 0.01 and a > 0):
        return f"{v:.3e}"
    return f"{v:.4g}"


def render_stage_table(rank: int, stage: int,
                       atv: Dict[str, List[float]],
                       lmp: Dict[str, List[float]]) -> Tuple[str, List[dict]]:
    lines = []
    lines.append(f"### Stage {stage}")
    lines.append("")
    lines.append("| metric | unit | Atomiverse μ±σ | LAMMPS μ±σ | Δμ | Δμ / 3√(σ_a²+σ_l²) | N_a | N_l |")
    lines.append("|---|---|---|---|---|---|---|---|")
    rows = []
    for name, col, unit in METRICS:
        a_series = atv.get(col, [])
        l_series = lmp.get(col, [])
        ma, sa, na = block_average(a_series)
        ml, sl, nl = block_average(l_series)
        d = ma - ml
        sig = math.sqrt(sa * sa + sl * sl)
        z = d / (3.0 * sig) if sig > 0 and math.isfinite(sig) else float("nan")
        lines.append(
            f"| {name} | {unit} "
            f"| {_fmt(ma)} ± {_fmt(sa)} "
            f"| {_fmt(ml)} ± {_fmt(sl)} "
            f"| {_fmt(d)} "
            f"| {_fmt(z)} "
            f"| {na} | {nl} |"
        )
        rows.append({
            "rank": rank, "stage": stage, "metric": name, "unit": unit,
            "atv_mean": ma, "atv_sigma": sa, "atv_n": na,
            "lmp_mean": ml, "lmp_sigma": sl, "lmp_n": nl,
            "delta_mean": d, "z_3sigma": z,
        })
    lines.append("")
    return "\n".join(lines), rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    repo_default = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p.add_argument("--lammps-runs",
                   default=os.path.join(repo_default,
                       "benchmarks/md/reference/polyethylene/lammps_baseline/runs"))
    p.add_argument("--atv-runs",
                   default=os.path.join(repo_default,
                       "benchmarks/md/reference/polyethylene/atomiverse_baseline/runs"))
    p.add_argument("--ranks", default="1,32,64,128",
                   help="comma-separated rank counts to compare (default: 1,32,64,128)")
    p.add_argument("--out-dir", default=None,
                   help="output directory (default: <atv-runs>/../comparison_<TS>)")
    args = p.parse_args(argv)

    try:
        ranks = [int(x) for x in args.ranks.split(",") if x.strip()]
    except ValueError:
        print(f"error: invalid --ranks: {args.ranks}", file=sys.stderr)
        return 2

    out_dir = args.out_dir
    if out_dir is None:
        import datetime as _dt
        ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = os.path.join(
            os.path.dirname(args.atv_runs.rstrip("/")),
            f"comparison_{ts}")
    os.makedirs(out_dir, exist_ok=True)
    md_path = os.path.join(out_dir, "comparison.md")
    csv_path = os.path.join(out_dir, "per_stage_deltas.csv")

    delta_rows: List[dict] = []
    sections: List[str] = []
    sections.append("# Atomiverse-CPU vs LAMMPS-CPU MSU PE baseline comparison")
    sections.append("")
    sections.append(f"- LAMMPS runs:    `{args.lammps_runs}`")
    sections.append(f"- Atomiverse runs:`{args.atv_runs}`")
    sections.append(f"- Ranks:          {ranks}")
    sections.append(f"- Generated:      {__import__('datetime').datetime.now().isoformat(timespec='seconds')}")
    sections.append("")
    sections.append("All LAMMPS values converted to Atomiverse canonical units "
                    "(atm→bar ×1.01325, kcal/mol→eV ×0.04336411). Block-average "
                    "σ via Flyvbjerg & Petersen 1989. `Δμ / 3√(σ_a²+σ_l²)` ≤ 1 "
                    "means agreement within 3σ.")
    sections.append("")

    timing_lines = ["## Wall-time per stage", "",
                    "| ranks | stage | LAMMPS wall_s | LAMMPS ns/day | LAMMPS ts/s |",
                    "|---|---|---|---|---|"]

    for R in ranks:
        lmp_dir = os.path.join(args.lammps_runs, f"n{R}r")
        atv_dir = os.path.join(args.atv_runs, f"n{R}r")
        lmp_thermo = os.path.join(lmp_dir, "thermo.csv")
        atv_thermo = os.path.join(atv_dir, "thermo.csv")
        if not os.path.isfile(lmp_thermo) or not os.path.isfile(atv_thermo):
            sections.append(f"## ranks = {R}")
            sections.append("")
            sections.append(f"**MISSING** — lmp_thermo={os.path.isfile(lmp_thermo)} "
                            f"atv_thermo={os.path.isfile(atv_thermo)}")
            sections.append("")
            continue

        atv = read_atomiverse_thermo(atv_thermo)
        lmp = read_lammps_thermo(lmp_thermo)
        lmp_sum = read_lammps_stage_summary(os.path.join(lmp_dir, "stage_summary.csv"))

        sections.append(f"## ranks = {R}")
        sections.append("")
        for s in STAGES:
            tbl, rows = render_stage_table(R, s, atv[s], lmp[s])
            sections.append(tbl)
            delta_rows.extend(rows)
            if s in lmp_sum:
                timing_lines.append(
                    f"| {R} | {s} "
                    f"| {lmp_sum[s].get('wall_seconds', '')} "
                    f"| {lmp_sum[s].get('ns_per_day', '')} "
                    f"| {lmp_sum[s].get('timesteps_per_second', '')} |")
        sections.append("")

    sections.append("")
    sections.extend(timing_lines)
    sections.append("")
    sections.append("Atomiverse per-stage timing isn't currently emitted by "
                    "polyethylene_tension's thermo.csv; the only Atomiverse "
                    "wall-time signal is total `wall_seconds` in each "
                    "run_env.txt. Per-stage Atomiverse timing requires a "
                    "follow-up extractor or instrumenting the benchmark to "
                    "log stage_start/end timestamps.")

    with open(md_path, "w") as fh:
        fh.write("\n".join(sections) + "\n")
    with open(csv_path, "w", newline="") as fh:
        if delta_rows:
            writer = csv.DictWriter(fh, fieldnames=list(delta_rows[0].keys()))
            writer.writeheader()
            for r in delta_rows:
                writer.writerow(r)

    print(f"[compare] wrote {md_path}")
    print(f"[compare] wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
