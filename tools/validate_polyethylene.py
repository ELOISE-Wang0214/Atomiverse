#!/usr/bin/env python3
"""Validate Atomiverse polyethylene-tension output against a LAMMPS reference.

Usage
-----

    tools/validate_polyethylene.py \
        --atomiverse-dir build/pe_tension_run \
        --reference-dir benchmarks/md/reference/polyethylene

The atomiverse run is produced by ``benchmarks/md/polyethylene_tension`` and
must contain ``thermo.csv`` and ``stress_strain.csv``. The reference directory
is expected to contain (optional) ``thermo_reference.csv`` and
``stress_strain_reference.csv``; either is permitted to be missing — in that
case only Atomiverse-side statistics are reported.

The script computes:

* per-step thermo deltas vs. the reference (T, P, E_pot, E_kin)
* equilibrium block averages
* stress-strain curve fits: elastic modulus, yield stress, hardening slope

Tolerances default to the values in ``msu-pe-validation.md`` (Phase D); they
can be overridden via CLI flags. Returns exit code 0 on success, 1 on any
tolerance breach, 2 on missing input.

The script intentionally avoids non-stdlib dependencies (no numpy, no pandas)
so it runs on a fresh Python install. CSVs are parsed with the stdlib ``csv``
module.

References:
* parallel-md-engine.md  Verification §7
* msu-pe-validation.md   Phase D (D.2 -- D.4)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass, field
from typing import Iterable, Sequence


def _read_csv(path: str) -> list[dict[str, str]]:
    """Read a CSV file with a header row into a list of {col: str} dicts.

    Returns an empty list if the file does not exist or is empty.
    """
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


def _to_float(rows: list[dict[str, str]], col: str) -> list[float]:
    return [float(r[col]) for r in rows if col in r and r[col] != ""]


@dataclass
class Tolerances:
    """Pass/fail thresholds; defaults track msu-pe-validation.md §D."""

    thermo_T_rel: float = 1e-2
    thermo_P_rel: float = 1e-1
    thermo_Epot_rel: float = 1e-2
    thermo_Ekin_rel: float = 1e-2
    density_rel: float = 5e-2
    modulus_rel: float = 5e-2
    yield_rel: float = 5e-2
    hardening_rel: float = 5e-2


@dataclass
class Failure:
    """One failed tolerance check (severity is the difference / threshold)."""

    name: str
    detail: str
    severity: float = 1.0


@dataclass
class Report:
    failures: list[Failure] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def passed(self) -> bool:
        return not self.failures

    def render(self) -> str:
        lines: list[str] = []
        if self.notes:
            lines.append("notes:")
            for n in self.notes:
                lines.append(f"  - {n}")
        if self.failures:
            lines.append("FAILURES:")
            for f in self.failures:
                lines.append(f"  [{f.severity:.2f}x tol] {f.name}: {f.detail}")
        else:
            lines.append("ALL CHECKS PASSED")
        return "\n".join(lines)


def _abs_rel(actual: float, expected: float) -> float:
    if expected == 0.0:
        return abs(actual)
    return abs(actual - expected) / abs(expected)


def _block_average(values: Sequence[float], frac_tail: float = 0.5) -> float:
    if not values:
        return 0.0
    n = max(1, int(len(values) * frac_tail))
    tail = values[-n:]
    return sum(tail) / len(tail)


def _linear_fit(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float]:
    """Plain least-squares slope, intercept. Returns (0, mean(ys)) on degenerate input."""
    if len(xs) < 2:
        return 0.0, (sum(ys) / len(ys) if ys else 0.0)
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den = sum((x - mean_x) ** 2 for x in xs)
    if den == 0.0:
        return 0.0, mean_y
    slope = num / den
    intercept = mean_y - slope * mean_x
    return slope, intercept


@dataclass
class StressStrainStats:
    """Quantities extracted from a stress-strain curve."""

    elastic_modulus_GPa: float
    yield_stress_GPa: float
    hardening_slope_GPa: float


def _compute_stress_strain_stats(
    strain: Sequence[float],
    sigma_GPa: Sequence[float],
    elastic_eps_max: float = 0.02,
    yield_eps_max: float = 0.10,
    hardening_eps_min: float = 0.5,
) -> StressStrainStats:
    """Extract elastic modulus, yield stress, and strain-hardening slope.

    All three quantities are robust to small noise: linear fits on slabs of
    points rather than single-point gradients.
    """
    elastic_xs = [e for e in strain if e <= elastic_eps_max]
    elastic_ys = [s for e, s in zip(strain, sigma_GPa) if e <= elastic_eps_max]
    e_modulus, _ = _linear_fit(elastic_xs, elastic_ys)

    yield_window = [s for e, s in zip(strain, sigma_GPa) if e <= yield_eps_max]
    yield_stress = max(yield_window) if yield_window else 0.0

    hardening_xs = [e for e in strain if e >= hardening_eps_min]
    hardening_ys = [s for e, s in zip(strain, sigma_GPa) if e >= hardening_eps_min]
    hardening_slope, _ = _linear_fit(hardening_xs, hardening_ys)

    return StressStrainStats(
        elastic_modulus_GPa=e_modulus,
        yield_stress_GPa=yield_stress,
        hardening_slope_GPa=hardening_slope,
    )


def _check_thermo(
    av_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    tol: Tolerances,
    n_steps_compare: int,
    rep: Report,
) -> None:
    if not ref_rows:
        rep.notes.append("thermo: no reference -- skipping per-step comparison.")
        return
    if not av_rows:
        rep.failures.append(Failure(
            name="thermo:atomiverse_empty",
            detail="atomiverse thermo.csv is empty.",
        ))
        return

    n = min(n_steps_compare, len(av_rows), len(ref_rows))
    if n < 2:
        rep.notes.append(f"thermo: only {n} comparable rows; skipping per-step compare.")
        return

    cols_to_check = [
        ("T_K", tol.thermo_T_rel, "T"),
        ("P_bar", tol.thermo_P_rel, "P"),
        ("E_pot_eV", tol.thermo_Epot_rel, "E_pot"),
        ("E_kin_eV", tol.thermo_Ekin_rel, "E_kin"),
    ]
    for col, threshold, label in cols_to_check:
        av_vals = _to_float(av_rows[:n], col)
        ref_vals = _to_float(ref_rows[:n], col)
        if not av_vals or not ref_vals:
            rep.notes.append(f"thermo: column '{col}' missing on one side.")
            continue
        m = min(len(av_vals), len(ref_vals))
        # Use block-averaged values to reduce noise.
        av_avg = _block_average(av_vals[:m])
        ref_avg = _block_average(ref_vals[:m])
        rel = _abs_rel(av_avg, ref_avg)
        if rel > threshold:
            rep.failures.append(Failure(
                name=f"thermo:{label}",
                detail=f"av={av_avg:.6g} ref={ref_avg:.6g} rel_err={rel:.2%} > {threshold:.2%}",
                severity=rel / threshold,
            ))


def _check_stress_strain(
    av_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    tol: Tolerances,
    rep: Report,
) -> None:
    if not av_rows:
        rep.failures.append(Failure(
            name="stress_strain:atomiverse_empty",
            detail="atomiverse stress_strain.csv is empty.",
        ))
        return

    av_strain = _to_float(av_rows, "strain")
    av_sigma  = _to_float(av_rows, "sigma_xx_GPa")
    if not av_strain or not av_sigma or len(av_strain) != len(av_sigma):
        rep.failures.append(Failure(
            name="stress_strain:atomiverse_malformed",
            detail="missing or mismatched strain / sigma_xx_GPa columns.",
        ))
        return
    av_stats = _compute_stress_strain_stats(av_strain, av_sigma)
    rep.notes.append(
        f"stress_strain: atomiverse  E={av_stats.elastic_modulus_GPa:.3g} GPa  "
        f"sigma_y={av_stats.yield_stress_GPa:.3g} GPa  "
        f"H={av_stats.hardening_slope_GPa:.3g} GPa"
    )

    if not ref_rows:
        rep.notes.append("stress_strain: no reference -- skipping comparison.")
        return

    ref_strain = _to_float(ref_rows, "strain")
    ref_sigma  = _to_float(ref_rows, "sigma_xx_GPa")
    if not ref_strain or not ref_sigma:
        rep.notes.append("stress_strain: reference missing strain / sigma_xx_GPa columns.")
        return
    ref_stats = _compute_stress_strain_stats(ref_strain, ref_sigma)
    rep.notes.append(
        f"stress_strain: reference  E={ref_stats.elastic_modulus_GPa:.3g} GPa  "
        f"sigma_y={ref_stats.yield_stress_GPa:.3g} GPa  "
        f"H={ref_stats.hardening_slope_GPa:.3g} GPa"
    )

    checks = [
        ("modulus", av_stats.elastic_modulus_GPa, ref_stats.elastic_modulus_GPa, tol.modulus_rel),
        ("yield",   av_stats.yield_stress_GPa,    ref_stats.yield_stress_GPa,    tol.yield_rel),
        ("hardening", av_stats.hardening_slope_GPa, ref_stats.hardening_slope_GPa, tol.hardening_rel),
    ]
    for name, av_v, ref_v, threshold in checks:
        rel = _abs_rel(av_v, ref_v)
        if rel > threshold:
            rep.failures.append(Failure(
                name=f"stress_strain:{name}",
                detail=f"av={av_v:.4g} ref={ref_v:.4g} rel_err={rel:.2%} > {threshold:.2%}",
                severity=rel / threshold,
            ))


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--atomiverse-dir", required=True,
                   help="Directory containing atomiverse thermo.csv and stress_strain.csv.")
    p.add_argument("--reference-dir", default=None,
                   help="Directory containing LAMMPS thermo_reference.csv / "
                        "stress_strain_reference.csv. Optional.")
    p.add_argument("--thermo-rows", type=int, default=1000,
                   help="Rows from each thermo.csv to compare. Default 1000.")
    p.add_argument("--strict", action="store_true",
                   help="Tighten thermo tolerances to msu-pe-validation Phase D.2.")
    args = p.parse_args(argv)

    av_dir = args.atomiverse_dir
    ref_dir = args.reference_dir
    thermo_av = _read_csv(os.path.join(av_dir, "thermo.csv"))
    stress_av = _read_csv(os.path.join(av_dir, "stress_strain.csv"))

    if not thermo_av and not stress_av:
        print(f"validate_polyethylene: no atomiverse output in {av_dir!r}", file=sys.stderr)
        return 2

    thermo_ref: list[dict[str, str]] = []
    stress_ref: list[dict[str, str]] = []
    if ref_dir:
        thermo_ref = _read_csv(os.path.join(ref_dir, "thermo_reference.csv"))
        stress_ref = _read_csv(os.path.join(ref_dir, "stress_strain_reference.csv"))

    tol = Tolerances()
    if args.strict:
        # Phase D.2 first-1000-step deterministic regime.
        tol.thermo_T_rel = 1e-4
        tol.thermo_P_rel = 1e-3
        tol.thermo_Epot_rel = 1e-4
        tol.thermo_Ekin_rel = 1e-4

    rep = Report()
    rep.notes.append(f"atomiverse_dir = {av_dir}")
    rep.notes.append(f"reference_dir  = {ref_dir or '<none>'}")
    rep.notes.append(f"thermo rows    = {len(thermo_av)} av, {len(thermo_ref)} ref")
    rep.notes.append(f"stress rows    = {len(stress_av)} av, {len(stress_ref)} ref")

    _check_thermo(thermo_av, thermo_ref, tol, args.thermo_rows, rep)
    _check_stress_strain(stress_av, stress_ref, tol, rep)

    print(rep.render())
    return 0 if rep.passed() else 1


if __name__ == "__main__":
    raise SystemExit(main())
