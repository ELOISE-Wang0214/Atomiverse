#!/usr/bin/env python3
"""Convert a `lammps_baseline/` run's `log.lammps` into two CSVs.

Inputs
------
  --log <path>    LAMMPS log file (default: log.lammps in --run-dir)
  --run-dir <D>   Run directory under benchmarks/md/reference/polyethylene/lammps_baseline/runs/.
                  CSVs are written into the same directory:
                    thermo.csv        — full per-step thermo, one row per
                                        thermo print, with stage_index column.
                    stage_summary.csv — one row per LAMMPS `run`: stage index,
                                        steps, wall-clock seconds, ns/day,
                                        timesteps/s, MPI ranks, OMP threads,
                                        per-section timing breakdown (Pair, Bond,
                                        Neigh, Comm, Output, Modify, Other) as
                                        percentages.

The reader is `molpy.io.log.lammps.read_LAMMPS_log` (molcrafts molpy).

Usage
-----
  tools/lammps_baseline_log_to_csv.py --run-dir <run_dir>
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# molpy lives outside the Atomiverse repo (molcrafts/molpy). The salamander
# venv at ~/salamander_2025/250501/.venv ships molpy + molrs; the user
# activates it before invoking this script. We just import.
try:
    from molpy.io.log.lammps import read_LAMMPS_log
except ModuleNotFoundError as exc:
    sys.stderr.write(
        "error: cannot import molpy.io.log.lammps -- activate the salamander\n"
        "       venv first: source ~/salamander_2025/250501/.venv/bin/activate\n"
        f"       ({exc})\n"
    )
    sys.exit(2)


def _section_pct(timing, section_name):
    """Return percent_total for the named section in a timing breakdown."""
    if timing is None:
        return None
    for row in timing.rows:
        if row.section.strip().lower() == section_name.lower():
            return row.percent_total
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="Run directory containing log.lammps; CSVs are written here.",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Override log file path (default: <run-dir>/log.lammps).",
    )
    args = p.parse_args(argv)

    log_path = args.log if args.log is not None else (args.run_dir / "log.lammps")
    if not log_path.is_file():
        sys.stderr.write(f"error: log file not found: {log_path}\n")
        return 1

    log = read_LAMMPS_log(log_path)

    if not log.runs:
        sys.stderr.write(f"error: no `run` blocks parsed from {log_path}\n")
        return 1

    # ── thermo.csv ────────────────────────────────────────────────────
    # Stages share the same thermo columns under the baseline deck.
    # Union of columns across runs (preserving first-seen order) makes
    # the schema robust to per-stage thermo_style customization.
    columns_union: list[str] = []
    for run in log.runs:
        if run.thermo is None:
            continue
        for col in run.thermo.columns:
            if col not in columns_union:
                columns_union.append(col)

    thermo_path = args.run_dir / "thermo.csv"
    with thermo_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["stage_index", *columns_union])
        for stage_index, run in enumerate(log.runs, start=1):
            if run.thermo is None:
                continue
            cols = run.thermo.columns
            for record in run.thermo.data.tolist():
                row_by_col = dict(zip(cols, record))
                writer.writerow(
                    [stage_index]
                    + [row_by_col.get(col, "") for col in columns_union]
                )

    # ── stage_summary.csv ─────────────────────────────────────────────
    summary_path = args.run_dir / "stage_summary.csv"
    summary_cols = [
        "stage_index",
        "n_steps",
        "n_atoms",
        "wall_seconds",
        "ns_per_day",
        "timesteps_per_second",
        "mpi_tasks",
        "omp_threads",
        "cpu_use_percent",
        "pct_pair",
        "pct_bond",
        "pct_neigh",
        "pct_comm",
        "pct_output",
        "pct_modify",
        "pct_other",
    ]
    with summary_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(summary_cols)
        for stage_index, run in enumerate(log.runs, start=1):
            loop = run.loop_time
            perf = run.performance
            cpu = run.CPU_use
            timing = run.MPI_task_timing
            writer.writerow(
                [
                    stage_index,
                    loop.steps if loop else "",
                    loop.atoms if loop else "",
                    loop.seconds if loop else "",
                    perf.ns_per_day if perf else "",
                    perf.timesteps_per_second if perf else "",
                    cpu.MPI_tasks if cpu else "",
                    cpu.OMP_threads if cpu else "",
                    cpu.percent if cpu else "",
                    _section_pct(timing, "Pair"),
                    _section_pct(timing, "Bond"),
                    _section_pct(timing, "Neigh"),
                    _section_pct(timing, "Comm"),
                    _section_pct(timing, "Output"),
                    _section_pct(timing, "Modify"),
                    _section_pct(timing, "Other"),
                ]
            )

    print(f"[lammps_baseline_log_to_csv] wrote:")
    print(f"  {thermo_path}")
    print(f"  {summary_path}")
    print(f"  {len(log.runs)} stage(s), schema cols = {columns_union}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
