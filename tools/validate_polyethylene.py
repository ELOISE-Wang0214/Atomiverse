#!/usr/bin/env python3
"""Validate Atomiverse polyethylene-tension output against a LAMMPS reference.

This is the Phase D comparison tool for ``.claude/specs/msu-pe-validation.md``.
It is pure stdlib Python (no numpy / pandas) so it runs on a fresh install.

Checks (selected by ``--check``)
--------------------------------

* ``thermo`` — statistical block-average thermo parity (supersedes the old
  per-step D.2 check). The MSU protocol's stage 1 is a stochastic Langevin
  thermostat, so the two engines run independent RNG streams and never share a
  trajectory — only an *ensemble*. Each thermo column is therefore reduced to
  contiguous block averages (Flyvbjerg & Petersen 1989) and the two engines are
  required to agree on the ensemble means of ``T`` (K) and ``E_pot`` (eV) to
  within ~2% relative. Runs for both the 1-rank and 8-rank pairs.
* ``equilibrium`` — Phase D.3 equilibrium-state observables of the stage-4
  (100 K NPT) glass: density ρ (1%), ⟨E_pot⟩ (0.5%), RDF first-peak position
  (1 bin / 0.05 Å), RDF first-peak height (3%), and end-to-end ⟨R²⟩ (10%).
* ``stress-strain`` — legacy stress-strain curve fits (elastic modulus, yield
  stress, hardening slope). Phase D.4 owns the unified ``--check`` dispatch.
* ``all`` — every check above.

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
* msu-pe-validation.md                 Phase D
* msu-pe-phase-d-01-reference-runs       (produces the reference CSVs)
* msu-pe-align-03-statistical-thermo     (statistical block-average thermo check)
* msu-pe-phase-d-03-equilibrium          (equilibrium-observables check)
* msu-pe-phase-d-04-stress-strain        (stress-strain check + run_checks dispatch)
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
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
    """Pass/fail thresholds for the stress-strain and equilibrium checks.

    Statistical thermo tolerances live in the separate frozen
    :class:`StatTolerances`; defaults track ``msu-pe-validation.md`` Phase D.
    """

    # stress-strain descriptors (Phase D.4) — relative
    modulus_rel: float = 5e-2
    yield_rel: float = 5e-2
    hardening_rel: float = 5e-2
    # equilibrium observables (Phase D.3)
    eq_density_rel: float = 1e-2            # relative
    eq_epot_rel: float = 5e-3              # relative
    eq_rdf_peak_pos_abs: float = 5e-2      # absolute, Å (1 histogram bin)
    eq_rdf_peak_height_rel: float = 3e-2  # relative
    eq_r2_rel: float = 1e-1               # relative


@dataclass(frozen=True)
class StatTolerances:
    """Per-column relative thresholds for the statistical block-average check.

    The MSU stage-1 protocol is a stochastic Langevin thermostat, so the two
    engines run uncorrelated RNG streams — only their *ensemble means* can be
    compared. Defaults are the contract bar: a ~2% relative-error bound on the
    block-averaged ensemble means of ``T`` (K) and ``E_pot`` (eV), loose enough
    to absorb finite-sample block-mean scatter, tight enough that a genuine
    force-field or thermostat-setpoint discrepancy fails. Both values are
    dimensionless relative-error bounds.
    """

    T: float = 2e-2
    E_pot: float = 2e-2

    def threshold(self, column: str) -> float:
        """Relative-error threshold for one compared thermo column."""
        return getattr(self, column)


#: Thermo columns whose ensemble mean is physically fixed and thus compared.
THERMO_COLUMNS: tuple[str, ...] = ("T", "E_pot")

#: (amu / Å³) → (g / cm³); numerator and denominator both carry a 1e-24 factor.
AMU_PER_A3_TO_G_PER_CM3 = 1.66053906660


@dataclass(frozen=True)
class BlockStat:
    """Block-averaged statistics for one thermo column.

    @field column    Canonical thermo column name (e.g. ``"T"``, ``"E_pot"``).
    @field mean      Mean of the per-block means — the ensemble-mean estimate
                     (K for ``T``, eV for ``E_pot``).
    @field sem       Standard error of the block means: ``stdev / sqrt(n)``,
                     the statistical uncertainty on ``mean`` (same unit).
    @field n_blocks  Number of full blocks; a trailing partial block is dropped.
    """

    column: str
    mean: float
    sem: float
    n_blocks: int


@dataclass(frozen=True)
class MeanGap:
    """Cross-engine block-averaged mean gap for one thermo column.

    @field column    Canonical thermo column name.
    @field av_mean   Atomiverse block-averaged ensemble mean.
    @field ref_mean  LAMMPS block-averaged ensemble mean.
    @field rel_gap   Relative gap ``|av - ref| / |ref|`` between the two means.
    """

    column: str
    av_mean: float
    ref_mean: float
    rel_gap: float


# ── Report plumbing ──────────────────────────────────────────────────────────


@dataclass
class Failure:
    """One failed tolerance check; ``severity`` is ``err / threshold``."""

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


# ── Statistical block-average thermo parity ──────────────────────────────────


def block_average(values: Sequence[float], block_size: int) -> BlockStat:
    """Reduce a thermo series to block-averaged statistics.

    Block averaging (Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461)
    partitions a time series into contiguous, equal-length blocks, takes each
    block's mean, and treats the block means as the sampling units. The mean of
    the block means estimates the ensemble mean; their standard error gives a
    variance estimate that accounts for serial correlation in the raw series.

    A trailing partial block shorter than ``block_size`` is dropped so every
    block carries equal weight — hence ``n_blocks == len(values) // block_size``.

    @param values      Numeric thermo series for one column (units per column:
                       K for ``T``, eV for ``E_pot``).
    @param block_size  Number of consecutive samples per block (>= 1).
    @return  A :class:`BlockStat` whose ``column`` is left empty (the caller
             tags it). ``sem`` is ``0.0`` when fewer than two full blocks exist
             (a single block has no spread to estimate).
    @raises ValueError  When ``block_size < 1``.
    """
    if block_size < 1:
        raise ValueError(f"block_size must be >= 1, got {block_size}")
    n_blocks = len(values) // block_size
    if n_blocks == 0:
        return BlockStat(column="", mean=0.0, sem=0.0, n_blocks=0)
    block_means = [
        statistics.fmean(values[b * block_size:(b + 1) * block_size])
        for b in range(n_blocks)
    ]
    mean = statistics.fmean(block_means)
    if n_blocks < 2:
        sem = 0.0
    else:
        sem = statistics.stdev(block_means) / (n_blocks ** 0.5)
    return BlockStat(column="", mean=mean, sem=sem, n_blocks=n_blocks)


def compare_thermo_statistical(
    av_rows: list[dict[str, str]],
    ref_rows: list[dict[str, str]],
    tol: StatTolerances,
    block_size: int,
) -> list[MeanGap]:
    """Compare two thermo CSVs by block-averaged ensemble means.

    The two engines run independent RNG streams under the same stochastic
    thermostat, so they sample the *same ensemble* but never the same
    trajectory. The physically meaningful parity statement is therefore
    agreement of the ensemble means, not of instantaneous per-step values.
    Each compared column (``T``, ``E_pot``) is block-averaged for both engines
    and the relative gap between the two means is gated at its
    :class:`StatTolerances` threshold. The comparison does not align engines
    step-by-step, so a differing logging cadence is tolerated provided each
    series has enough samples to form at least one block.

    @param av_rows     Atomiverse thermo rows (``{column: str}`` dicts).
    @param ref_rows    LAMMPS reference thermo rows.
    @param tol         Per-column relative-error thresholds on the means.
    @param block_size  Samples per block, passed to :func:`block_average`.
    @return  One :class:`MeanGap` per column whose relative mean gap exceeds
             its threshold. Empty list ⇒ statistical parity holds.
    """
    breaches: list[MeanGap] = []
    for col in THERMO_COLUMNS:
        av_stat = block_average(_to_float(av_rows, col), block_size)
        ref_stat = block_average(_to_float(ref_rows, col), block_size)
        rel_gap = _abs_rel(av_stat.mean, ref_stat.mean)
        if rel_gap > tol.threshold(col):
            breaches.append(MeanGap(
                column=col,
                av_mean=av_stat.mean,
                ref_mean=ref_stat.mean,
                rel_gap=rel_gap,
            ))
    return breaches


def run_thermo_check(
    pairs: Sequence[tuple[str, str, str]],
    tol: StatTolerances,
    block_size: int,
) -> Report:
    """Run the statistical block-average thermo comparison over labelled pairs.

    @param pairs       Sequence of ``(label, av_csv_path, ref_csv_path)`` — the
                       label is a rank tag such as ``"1rank"`` / ``"8rank"``.
    @param tol         Per-column relative-error thresholds on the means.
    @param block_size  Samples per block, passed to
                       :func:`compare_thermo_statistical`.
    @return  A :class:`Report`; ``passed()`` is true when every pair agrees on
             the block-averaged ensemble means of ``T`` and ``E_pot``.
    """
    rep = Report()
    for label, av_csv, ref_csv in pairs:
        av_rows = _read_csv(av_csv)
        ref_rows = _read_csv(ref_csv)
        breaches = compare_thermo_statistical(av_rows, ref_rows, tol, block_size)
        if not breaches:
            rep.notes.append(
                f"thermo[{label}]: block-averaged <T>/<E_pot> agree within "
                f"tolerance (block_size={block_size})"
            )
            continue
        for gap in breaches:
            threshold = tol.threshold(gap.column)
            rep.failures.append(Failure(
                name=f"thermo:{label}:{gap.column}",
                detail=(f"<{gap.column}> av={gap.av_mean:.4g} "
                        f"ref={gap.ref_mean:.4g} rel_gap={gap.rel_gap:.2%} "
                        f"> {threshold:.2%}"),
                severity=gap.rel_gap / threshold,
            ))
    return rep


# ── Phase D.3 — equilibrium-state observables ────────────────────────────────


@dataclass(frozen=True)
class EquilibriumStats:
    """Equilibrium-state observables of the stage-4 (100 K NPT) glass.

    @field density_g_cm3    Mass density, g/cm³.
    @field epot_mean_eV     Block-averaged potential energy, eV.
    @field rdf_peak_pos_A   First RDF peak position r*, Å.
    @field rdf_peak_height  First RDF peak height g(r*), dimensionless.
    @field r2_mean_A2       Mean squared end-to-end chain distance ⟨R²⟩, Å².
    """

    density_g_cm3: float
    epot_mean_eV: float
    rdf_peak_pos_A: float
    rdf_peak_height: float
    r2_mean_A2: float


def _read_xyz(path: str) -> list[list[tuple[float, float, float]]]:
    """Parse a multi-frame plain-XYZ trajectory into a list of frames.

    Each frame is a list of ``(x, y, z)`` tuples in Å. The expected per-frame
    layout (as written by ``cpu::trajectory_append``) is::

        <n_atoms>
        <comment line>
        <type> <x> <y> <z>
        ...   (n_atoms coordinate lines)

    @raises ValueError  On a truncated frame or a non-numeric coordinate.
    """
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    frames: list[list[tuple[float, float, float]]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() == "":
            i += 1
            continue
        try:
            n_atoms = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(
                f"{path}: expected an atom count on line {i + 1}") from exc
        if i + 2 + n_atoms > len(lines):
            raise ValueError(
                f"{path}: truncated frame at line {i + 1} "
                f"(declared {n_atoms} atoms)")
        atoms: list[tuple[float, float, float]] = []
        for j in range(n_atoms):
            parts = lines[i + 2 + j].split()
            if len(parts) < 4:
                raise ValueError(
                    f"{path}: malformed atom line {i + 3 + j}")
            atoms.append((float(parts[1]), float(parts[2]), float(parts[3])))
        frames.append(atoms)
        i += 2 + n_atoms
    return frames


def _rdf_first_peak(
    r_values: Sequence[float],
    g_values: Sequence[float],
) -> tuple[float, float]:
    """First interior local maximum of g(r), skipping the r≈0 noise floor.

    @return  ``(r_star, g_star)`` of the first peak, or ``(0.0, 0.0)`` when no
             interior peak exists — a degenerate g(r) the caller turns into a
             tolerance breach rather than a silent zero.
    """
    n = min(len(r_values), len(g_values))
    for i in range(1, n - 1):
        if (g_values[i] > 0.0
                and g_values[i] > g_values[i - 1]
                and g_values[i] >= g_values[i + 1]):
            return (r_values[i], g_values[i])
    return (0.0, 0.0)


def _end_to_end_r2(
    frames: Sequence[Sequence[tuple[float, float, float]]],
    chain_ranges: Sequence[tuple[int, int]],
) -> float:
    """Mean squared end-to-end distance ⟨R²⟩ in Å², over chains and frames.

    @param chain_ranges  ``(first_atom, last_atom)`` inclusive index pairs.
    @raises ValueError   When a chain index exceeds a frame's atom count.
    """
    total = 0.0
    count = 0
    for frame in frames:
        for first, last in chain_ranges:
            if first >= len(frame) or last >= len(frame):
                raise ValueError(
                    f"chain range ({first},{last}) outside a frame of "
                    f"{len(frame)} atoms")
            x0, y0, z0 = frame[first]
            x1, y1, z1 = frame[last]
            total += (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2
            count += 1
    return total / count if count else 0.0


def compute_equilibrium_stats(
    thermo_rows: list[dict[str, str]],
    rdf_table: tuple[Sequence[float], Sequence[float]],
    frames: Sequence[Sequence[tuple[float, float, float]]],
    chain_ranges: Sequence[tuple[int, int]],
    total_mass_amu: float,
) -> EquilibriumStats:
    """Assemble the five equilibrium observables for one engine.

    @param thermo_rows     Stage-4 thermo rows; must carry ``volume`` (Å³) and
                           ``E_pot`` (eV) columns.
    @param rdf_table       ``(r[Å], g(r))`` finalized RDF curve.
    @param frames          Stage-4 trajectory frames (Å).
    @param chain_ranges    Per-chain ``(first, last)`` terminal-bead indices.
    @param total_mass_amu  Total system mass, amu.
    """
    v_mean = _block_average(_to_float(thermo_rows, "volume"))
    density = (total_mass_amu / v_mean * AMU_PER_A3_TO_G_PER_CM3
               if v_mean > 0.0 else 0.0)
    epot_mean = _block_average(_to_float(thermo_rows, "E_pot"))
    r_star, g_star = _rdf_first_peak(rdf_table[0], rdf_table[1])
    r2 = _end_to_end_r2(frames, chain_ranges)
    return EquilibriumStats(
        density_g_cm3=density,
        epot_mean_eV=epot_mean,
        rdf_peak_pos_A=r_star,
        rdf_peak_height=g_star,
        r2_mean_A2=r2,
    )


def _check_equilibrium(
    av: EquilibriumStats,
    ref: EquilibriumStats,
    tol: Tolerances,
    rep: Report,
) -> None:
    """Gate the five Phase D.3 equilibrium observables into ``rep``.

    Density / ⟨E_pot⟩ / RDF peak height / ⟨R²⟩ use relative tolerances; the
    RDF first-peak *position* uses an absolute tolerance (one histogram bin).
    """
    rep.notes.append(
        f"equilibrium: atomiverse  rho={av.density_g_cm3:.4g} g/cm3  "
        f"E_pot={av.epot_mean_eV:.4g} eV  r*={av.rdf_peak_pos_A:.3g} A  "
        f"g(r*)={av.rdf_peak_height:.3g}  <R2>={av.r2_mean_A2:.4g} A2")
    rep.notes.append(
        f"equilibrium: reference   rho={ref.density_g_cm3:.4g} g/cm3  "
        f"E_pot={ref.epot_mean_eV:.4g} eV  r*={ref.rdf_peak_pos_A:.3g} A  "
        f"g(r*)={ref.rdf_peak_height:.3g}  <R2>={ref.r2_mean_A2:.4g} A2")
    checks = [
        ("density", av.density_g_cm3, ref.density_g_cm3,
         tol.eq_density_rel, "rel"),
        ("epot", av.epot_mean_eV, ref.epot_mean_eV,
         tol.eq_epot_rel, "rel"),
        ("rdf_peak_pos", av.rdf_peak_pos_A, ref.rdf_peak_pos_A,
         tol.eq_rdf_peak_pos_abs, "abs"),
        ("rdf_peak_height", av.rdf_peak_height, ref.rdf_peak_height,
         tol.eq_rdf_peak_height_rel, "rel"),
        ("r2", av.r2_mean_A2, ref.r2_mean_A2,
         tol.eq_r2_rel, "rel"),
    ]
    for name, av_v, ref_v, threshold, kind in checks:
        err = abs(av_v - ref_v) if kind == "abs" else _abs_rel(av_v, ref_v)
        if err > threshold:
            rep.failures.append(Failure(
                name=f"equilibrium:{name}",
                detail=(f"av={av_v:.4g} ref={ref_v:.4g} "
                        f"{kind}_err={err:.3g} > {threshold:.3g}"),
                severity=err / threshold,
            ))


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


def _stress_sigma(rows: list[dict[str, str]]) -> list[float]:
    """Stress column from a stress-strain CSV (GPa).

    Accepts the canonical reference schema (``sigma_xx``) and the benchmark's
    native schema (``sigma_xx_GPa``).
    """
    for col in ("sigma_xx", "sigma_xx_GPa"):
        values = _to_float(rows, col)
        if values:
            return values
    return []


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
    av_sigma = _stress_sigma(av_rows)
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
    ref_sigma = _stress_sigma(ref_rows)
    if not ref_strain or not ref_sigma:
        rep.notes.append("stress_strain: reference missing strain / sigma_xx.")
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

#: Filenames the equilibrium check reads from each engine's directory.
_EQ_THERMO_FILE = "equilibrium_thermo.csv"
_EQ_RDF_FILE = "rdf.csv"
_EQ_TRAJ_FILE = "trajectory.xyz"


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


def _load_equilibrium(
    directory: str,
    label: str,
    n_chains: int,
    chain_length: int,
    bead_mass_amu: float,
    rep: Report,
) -> EquilibriumStats | None:
    """Load one engine's equilibrium inputs; ``None`` when an input is absent.

    A *missing* input file is a graceful skip (recorded as a note); a file that
    is present but *malformed* is a reported failure.
    """
    thermo_path = os.path.join(directory, _EQ_THERMO_FILE)
    rdf_path = os.path.join(directory, _EQ_RDF_FILE)
    traj_path = os.path.join(directory, _EQ_TRAJ_FILE)
    for path in (thermo_path, rdf_path, traj_path):
        if not os.path.isfile(path):
            rep.notes.append(
                f"equilibrium[{label}]: missing {os.path.basename(path)} "
                f"-- skipping equilibrium check")
            return None

    rdf_rows = _read_csv(rdf_path)
    rdf_table = (_to_float(rdf_rows, "r"), _to_float(rdf_rows, "g"))
    chain_ranges = [
        (c * chain_length, c * chain_length + chain_length - 1)
        for c in range(n_chains)
    ]
    total_mass = float(n_chains * chain_length) * bead_mass_amu
    try:
        frames = _read_xyz(traj_path)
        return compute_equilibrium_stats(
            _read_csv(thermo_path), rdf_table, frames, chain_ranges, total_mass)
    except ValueError as exc:
        rep.failures.append(Failure(
            name=f"equilibrium:{label}:malformed_input",
            detail=str(exc)))
        return None


def run_stress_strain_check(
    ref_dir: str,
    ranks: int,
    tol: Tolerances,
    rep: Report,
) -> int | None:
    """Phase D.4 stress-strain comparison for one rank count.

    Reads the canonical ``stress_strain_<engine>_<N>r.csv`` curves from
    ``ref_dir`` and gates the three mechanical descriptors (elastic modulus,
    yield stress, hardening slope) at ±5%.

    @return  Exit code ``2`` when the LAMMPS reference CSV is absent (an
             operator error for this gate); otherwise ``None`` with any
             tolerance failures accumulated into ``rep``.
    """
    av_csv = os.path.join(ref_dir, f"stress_strain_atomiverse_{ranks}r.csv")
    ref_csv = os.path.join(ref_dir, f"stress_strain_lammps_{ranks}r.csv")
    if not os.path.isfile(ref_csv):
        print(f"validate_polyethylene: missing stress-strain reference "
              f"{os.path.basename(ref_csv)} under {ref_dir!r}", file=sys.stderr)
        return 2
    _check_stress_strain(_read_csv(av_csv), _read_csv(ref_csv), tol, rep)
    return None


def run_checks(
    check: str,
    args: argparse.Namespace,
    tol: Tolerances,
    tol_thermo: StatTolerances,
    rep: Report,
) -> int | None:
    """Dispatch the selected ``--check`` value(s) to their check functions.

    Routes ``thermo`` → statistical block-average parity, ``equilibrium`` →
    D.3, ``stress-strain`` → D.4, and ``all`` → every check, accumulating
    failures into ``rep``.

    @return  Exit code ``2`` when a required input file is missing; otherwise
             ``None`` (``rep`` then carries any tolerance failures).
    """
    if check in ("thermo", "all"):
        pairs = _thermo_pairs(args.reference_dir)
        for label, av_csv, ref_csv in pairs:
            if not os.path.isfile(av_csv) or not os.path.isfile(ref_csv):
                print(f"validate_polyethylene: missing thermo reference CSV "
                      f"for {label} under {args.reference_dir!r}",
                      file=sys.stderr)
                return 2
        thermo_rep = run_thermo_check(pairs, tol_thermo, args.thermo_block_size)
        rep.failures.extend(thermo_rep.failures)
        rep.notes.extend(thermo_rep.notes)

    if check in ("stress-strain", "all"):
        rc = run_stress_strain_check(args.reference_dir, args.ranks or 1,
                                     tol, rep)
        if rc is not None:
            return rc

    if check in ("equilibrium", "all"):
        if not args.atomiverse_dir:
            print("validate_polyethylene: --check {equilibrium,all} requires "
                  "--atomiverse-dir", file=sys.stderr)
            return 2
        av_eq = _load_equilibrium(args.atomiverse_dir, "atomiverse",
                                  args.n_chains, args.chain_length,
                                  args.bead_mass, rep)
        ref_eq = _load_equilibrium(args.reference_dir, "reference",
                                   args.n_chains, args.chain_length,
                                   args.bead_mass, rep)
        if av_eq is not None and ref_eq is not None:
            _check_equilibrium(av_eq, ref_eq, tol, rep)

    return None


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. See the module docstring for the exit-code contract."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--check",
                   choices=["all", "thermo", "stress-strain", "equilibrium"],
                   default="all",
                   help="Which check(s) to run (default: all).")
    p.add_argument("--reference-dir",
                   default="benchmarks/md/reference/polyethylene",
                   help="Directory holding the canonical reference CSVs / "
                        "equilibrium inputs.")
    p.add_argument("--atomiverse-dir", default=None,
                   help="Directory with atomiverse equilibrium inputs. "
                        "Required for --check {equilibrium,all}.")
    p.add_argument("--thermo-block-size", type=int, default=10,
                   help="Samples per block for the statistical thermo compare "
                        "(default 10). Each thermo column is reduced to "
                        "contiguous block averages and the engines' ensemble "
                        "means are compared. Raise it for longer reference "
                        "series to better average out serial correlation.")
    p.add_argument("--ranks", type=int, choices=[1, 8], default=None,
                   help="Rank count this run covers; recorded in the report "
                        "and selects the stress-strain CSV suffix (default 1).")
    p.add_argument("--n-chains", type=int, default=10,
                   help="Chains in the system, for end-to-end R^2 (default 10).")
    p.add_argument("--chain-length", type=int, default=1000,
                   help="United atoms per chain (default 1000).")
    p.add_argument("--bead-mass", type=float, default=14.027,
                   help="Per-bead mass in amu (default 14.027 -- CH2).")
    args = p.parse_args(argv)

    rep = Report()
    if args.ranks is not None:
        rep.notes.append(f"rank count: {args.ranks}")

    rc = run_checks(args.check, args, Tolerances(), StatTolerances(), rep)
    if rc is not None:
        return rc

    print(rep.render())
    return 0 if rep.passed() else 1


if __name__ == "__main__":
    raise SystemExit(main())
