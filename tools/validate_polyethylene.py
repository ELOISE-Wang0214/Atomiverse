#!/usr/bin/env python3
"""Validate Atomiverse polyethylene-tension output against a LAMMPS reference.

This is the Phase D comparison tool for ``.claude/specs/msu-pe-validation.md``.
It is pure stdlib Python (no numpy / pandas) so it runs on a fresh install.

Checks (selected by ``--check``)
--------------------------------

* ``thermo`` — Phase D.2 per-step thermo parity. Over the first ``--thermo-rows``
  steps of the deterministic regime, compare the Atomiverse and LAMMPS canonical
  thermo CSVs step-by-step: relative agreement of 1e-4 on ``T``/``E_pot``/``E_kin``
  and 1e-3 on ``Pxx``/``Pyy``/``Pzz``. Runs for both the 1-rank and 8-rank pairs.
* ``all`` — the ``thermo`` check plus the legacy stress-strain curve fits
  (elastic modulus, yield stress, hardening slope). Phase D.3 (equilibrium) and
  D.4 (the unified ``--check`` dispatch) extend this file further.

Canonical thermo CSV schema (produced by Phase D.1's reference runs)::

    step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz

Reference CSVs live under ``--reference-dir`` (default
``benchmarks/md/reference/polyethylene``) with the engine- and rank-suffixed
names ``thermo_<engine>_<N>r.csv`` written by ``tools/run_msu_<engine>.sh``.

Exit codes
----------

* ``0`` — every selected check passed.
* ``1`` — a tolerance was breached.
* ``2`` — a required input file was missing.

References:
* msu-pe-validation.md          Phase D
* msu-pe-phase-d-01-reference-runs  (produces the reference CSVs)
* msu-pe-phase-d-02-thermo-parity   (this per-step thermo check)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass, field
from typing import Sequence


# ── CSV helpers ──────────────────────────────────────────────────────────────


def _read_csv(path: str) -> list[dict[str, str]]:
    """Read a header CSV into a list of ``{column: value}`` dicts.

    Returns an empty list when the file is absent or empty.
    """
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return [row for row in reader]


def _to_float(rows: list[dict[str, str]], col: str) -> list[float]:
    """Extract one numeric column from CSV rows, skipping blank cells."""
    return [float(r[col]) for r in rows if col in r and r[col] != ""]


# ── Tolerances ───────────────────────────────────────────────────────────────


@dataclass
class Tolerances:
    """Pass/fail thresholds for the non-thermo checks (stress-strain, density).

    Thermo tolerances live in the separate frozen :class:`ThermoTolerances`;
    defaults track ``msu-pe-validation.md`` Phase D.
    """

    density_rel: float = 1e-2
    modulus_rel: float = 5e-2
    yield_rel: float = 5e-2
    hardening_rel: float = 5e-2


@dataclass(frozen=True)
class ThermoTolerances:
    """Per-column relative thresholds for the Phase D.2 per-step thermo check.

    Defaults are the Phase D.2 contract: 1e-4 on the smoothly-varying state
    quantities ``T``/``E_pot``/``E_kin`` and the looser 1e-3 on the
    virial-derived pressure components ``Pxx``/``Pyy``/``Pzz``. All values are
    dimensionless relative-error bounds.
    """

    T: float = 1e-4
    E_pot: float = 1e-4
    E_kin: float = 1e-4
    Pxx: float = 1e-3
    Pyy: float = 1e-3
    Pzz: float = 1e-3

    def threshold(self, column: str) -> float:
        """Relative-error threshold for one canonical thermo column."""
        return getattr(self, column)


#: Canonical thermo columns compared per step, in report order.
THERMO_COLUMNS: tuple[str, ...] = ("T", "E_pot", "E_kin", "Pxx", "Pyy", "Pzz")


@dataclass(frozen=True)
class StepDeviation:
    """The single worst per-step breach for one thermo column.

    ``column`` is a canonical thermo column, or ``"step_misalignment"`` when a
    step is present in only one of the two CSVs (``rel_err`` is then infinite).
    """

    column: str
    step: int
    rel_err: float


# ── Report plumbing ──────────────────────────────────────────────────────────


@dataclass
class Failure:
    """One failed tolerance check; ``severity`` is ``rel_err / threshold``."""

    name: str
    detail: str
    severity: float = 1.0


@dataclass
class Report:
    """Accumulated failures + informational notes for one validation run."""

    failures: list[Failure] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def passed(self) -> bool:
        """True when no tolerance was breached."""
        return not self.failures

    def render(self) -> str:
        """Format the report for stdout."""
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
    """Relative error |actual - expected| / |expected|; absolute when expected==0."""
    if expected == 0.0:
        return abs(actual)
    return abs(actual - expected) / abs(expected)


def _block_average(values: Sequence[float], frac_tail: float = 0.5) -> float:
    """Mean of the trailing ``frac_tail`` fraction of a series (0.0 when empty)."""
    if not values:
        return 0.0
    n = max(1, int(len(values) * frac_tail))
    tail = values[-n:]
    return sum(tail) / len(tail)


def _linear_fit(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float]:
    """Least-squares ``(slope, intercept)``; ``(0, mean(ys))`` on degenerate input."""
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


# ── Phase D.2 — per-step thermo parity ───────────────────────────────────────


def _rows_by_step(rows: list[dict[str, str]], n_steps: int) -> dict[int, dict[str, str]]:
    """Index thermo rows by integer ``step``, keeping only ``step <= n_steps``."""
    indexed: dict[int, dict[str, str]] = {}
    for r in rows:
        if "step" not in r or r["step"] == "":
            continue
        step = int(float(r["step"]))
        if step <= n_steps:
            indexed[step] = r
    return indexed


def compare_thermo_per_step(
    av_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    tol: ThermoTolerances,
    n_steps: int,
) -> list[StepDeviation]:
    """Compare two canonical thermo CSVs step-by-step over the first ``n_steps``.

    The two CSVs are aligned by the explicit ``step`` column, not by row index,
    so a mismatched logging cadence (a step present in one CSV only) is itself
    reported as a breach rather than silently skewing the comparison.

    @param av_rows   Atomiverse thermo rows (``{column: str}`` dicts).
    @param ref_rows  LAMMPS reference thermo rows.
    @param tol       Per-column relative-error thresholds.
    @param n_steps   Inclusive upper bound on the ``step`` window to compare.
    @return  One :class:`StepDeviation` per breached column (its worst step),
             plus one per misaligned step. Empty list ⇒ per-step parity holds.
    """
    av_by_step = _rows_by_step(av_rows, n_steps)
    ref_by_step = _rows_by_step(ref_rows, n_steps)

    deviations: list[StepDeviation] = []

    # Step misalignment — a step in one CSV but not the other is a failure.
    for step in sorted(set(av_by_step) ^ set(ref_by_step)):
        deviations.append(StepDeviation("step_misalignment", step, float("inf")))

    # Per-column worst breach over the common steps.
    common = sorted(set(av_by_step) & set(ref_by_step))
    for col in THERMO_COLUMNS:
        threshold = tol.threshold(col)
        worst: StepDeviation | None = None
        for step in common:
            av_cell = av_by_step[step].get(col, "")
            ref_cell = ref_by_step[step].get(col, "")
            if av_cell == "" or ref_cell == "":
                continue
            rel = _abs_rel(float(av_cell), float(ref_cell))
            if rel > threshold and (worst is None or rel > worst.rel_err):
                worst = StepDeviation(col, step, rel)
        if worst is not None:
            deviations.append(worst)

    return deviations


def run_thermo_check(
    pairs: Sequence[tuple[str, str, str]],
    tol: ThermoTolerances,
    n_steps: int,
) -> Report:
    """Run the per-step thermo comparison over a set of labelled CSV pairs.

    @param pairs    Sequence of ``(label, av_csv_path, ref_csv_path)`` — the
                    label is a rank tag such as ``"1rank"`` / ``"8rank"``.
    @param tol      Per-column relative-error thresholds.
    @param n_steps  Inclusive ``step`` window passed to
                    :func:`compare_thermo_per_step`.
    @return  A :class:`Report`; ``passed()`` is true when every pair agrees.
    """
    rep = Report()
    for label, av_csv, ref_csv in pairs:
        av_rows = _read_csv(av_csv)
        ref_rows = _read_csv(ref_csv)
        deviations = compare_thermo_per_step(av_rows, ref_rows, tol, n_steps)
        if not deviations:
            rep.notes.append(
                f"thermo[{label}]: per-step parity OK over first {n_steps} steps"
            )
            continue
        for dev in deviations:
            if dev.column == "step_misalignment":
                rep.failures.append(Failure(
                    name=f"thermo:{label}:step_misalignment",
                    detail=f"step {dev.step} present in only one CSV",
                    severity=float("inf"),
                ))
            else:
                threshold = tol.threshold(dev.column)
                rep.failures.append(Failure(
                    name=f"thermo:{label}:{dev.column}",
                    detail=(f"worst step {dev.step}: rel_err={dev.rel_err:.3e} "
                            f"> {threshold:.1e}"),
                    severity=dev.rel_err / threshold,
                ))
    return rep


# ── Stress-strain curve fits (legacy; Phase D.4 owns the comparison) ─────────


@dataclass
class StressStrainStats:
    """Quantities extracted from a stress-strain curve (all in GPa)."""

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

    All three are robust to small noise: linear fits over slabs of points
    rather than single-point gradients.
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


def _check_stress_strain(
    av_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    tol: Tolerances,
    rep: Report,
) -> None:
    """Compare Atomiverse vs reference stress-strain descriptors into ``rep``."""
    if not av_rows:
        rep.failures.append(Failure(
            name="stress_strain:atomiverse_empty",
            detail="atomiverse stress_strain.csv is empty.",
        ))
        return

    av_strain = _to_float(av_rows, "strain")
    av_sigma = _to_float(av_rows, "sigma_xx_GPa")
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
    ref_sigma = _to_float(ref_rows, "sigma_xx_GPa")
    if not ref_strain or not ref_sigma:
        rep.notes.append("stress_strain: reference missing strain / sigma_xx_GPa.")
        return
    ref_stats = _compute_stress_strain_stats(ref_strain, ref_sigma)
    rep.notes.append(
        f"stress_strain: reference  E={ref_stats.elastic_modulus_GPa:.3g} GPa  "
        f"sigma_y={ref_stats.yield_stress_GPa:.3g} GPa  "
        f"H={ref_stats.hardening_slope_GPa:.3g} GPa"
    )

    checks = [
        ("modulus", av_stats.elastic_modulus_GPa,
         ref_stats.elastic_modulus_GPa, tol.modulus_rel),
        ("yield", av_stats.yield_stress_GPa,
         ref_stats.yield_stress_GPa, tol.yield_rel),
        ("hardening", av_stats.hardening_slope_GPa,
         ref_stats.hardening_slope_GPa, tol.hardening_rel),
    ]
    for name, av_v, ref_v, threshold in checks:
        rel = _abs_rel(av_v, ref_v)
        if rel > threshold:
            rep.failures.append(Failure(
                name=f"stress_strain:{name}",
                detail=(f"av={av_v:.4g} ref={ref_v:.4g} "
                        f"rel_err={rel:.2%} > {threshold:.2%}"),
                severity=rel / threshold,
            ))


# ── CLI ──────────────────────────────────────────────────────────────────────

#: Rank tags and the (label, suffix) used to resolve canonical thermo CSVs.
_RANK_TAGS: tuple[tuple[str, str], ...] = (("1rank", "1r"), ("8rank", "8r"))


def _thermo_pairs(ref_dir: str) -> list[tuple[str, str, str]]:
    """Resolve the (label, atomiverse_csv, lammps_csv) thermo pairs in ``ref_dir``."""
    pairs: list[tuple[str, str, str]] = []
    for label, suffix in _RANK_TAGS:
        pairs.append((
            label,
            os.path.join(ref_dir, f"thermo_atomiverse_{suffix}.csv"),
            os.path.join(ref_dir, f"thermo_lammps_{suffix}.csv"),
        ))
    return pairs


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. See the module docstring for the exit-code contract."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--check", choices=["thermo", "all"], default="all",
                   help="Which check(s) to run (default: all).")
    p.add_argument("--reference-dir",
                   default="benchmarks/md/reference/polyethylene",
                   help="Directory holding the canonical thermo_<engine>_<N>r.csv "
                        "reference CSVs.")
    p.add_argument("--atomiverse-dir", default=None,
                   help="Directory with atomiverse thermo.csv / stress_strain.csv. "
                        "Required for --check all.")
    p.add_argument("--thermo-rows", type=int, default=1000,
                   help="Inclusive step window for the per-step thermo compare "
                        "(default 1000 -- the deterministic regime).")
    args = p.parse_args(argv)

    tol_thermo = ThermoTolerances()
    tol = Tolerances()

    # ── Phase D.2 per-step thermo parity ──
    pairs = _thermo_pairs(args.reference_dir)
    for label, av_csv, ref_csv in pairs:
        if not os.path.isfile(av_csv) or not os.path.isfile(ref_csv):
            print(f"validate_polyethylene: missing thermo reference CSV for "
                  f"{label} under {args.reference_dir!r}", file=sys.stderr)
            return 2
    rep = run_thermo_check(pairs, tol_thermo, args.thermo_rows)

    # ── legacy stress-strain (only under --check all) ──
    if args.check == "all":
        if not args.atomiverse_dir:
            print("validate_polyethylene: --check all requires --atomiverse-dir",
                  file=sys.stderr)
            return 2
        stress_av = _read_csv(os.path.join(args.atomiverse_dir,
                                           "stress_strain.csv"))
        stress_ref = _read_csv(os.path.join(args.reference_dir,
                                            "stress_strain_lammps_1r.csv"))
        _check_stress_strain(stress_av, stress_ref, tol, rep)

    print(rep.render())
    return 0 if rep.passed() else 1


if __name__ == "__main__":
    raise SystemExit(main())
