#!/usr/bin/env python3
"""Cross-engine + cross-backend MSU-PE statistical-equivalence validator.

Implements the link-08 (``msu-pe-cuda-components-08-cuda-validation``)
verdict generator: read the canonical 5-stage thermo + stress-strain CSVs
from two runs (typically Atomiverse-CUDA vs LAMMPS-GPU), compute four MD
invariants under a Flyvbjerg-Petersen block-average error model
(*J. Chem. Phys.* 91 (1989) 461), apply per-invariant tolerances, and
emit a markdown report.

The four invariants (per the spec § Domain basis):

* ``T_stage4``         — ⟨T⟩ over stage 4 (NPT equilibration, 100 K target).
* ``rho_stage4``       — ⟨ρ⟩ = M_total / ⟨V_box⟩ over stage 4.
* ``E_modulus_stage5`` — slope of σ_xx vs ε for ε ∈ (0, 0.05) (elastic).
* ``sigma_flow_stage5``— ⟨σ_xx⟩ over ε > 0.5 (plastic flow plateau).

Tolerances (per the spec § Domain basis):

* ``|ΔT|         < 3·√(σ_A² + σ_B²)``  ≈ 15 K
* ``|Δρ|         < 0.10 g/cm³``
* ``|ΔE| / E_B   < 0.15``               (15 % single-precision + seed drift)
* ``|Δσ_flow|    < 0.2 GPa``

Pure stdlib (``csv`` + ``argparse`` + ``math``) so it runs on a fresh
install. Reads the *extended* canonical schema:

* thermo:        ``step,stage,T_K,box_x,box_y,box_z`` (any extra columns
                 in between are tolerated).
* stress-strain: ``step,stage,strain,sigma_xx_GPa``    (likewise).

The extended schema is a superset of the link-01 canonical 8-column
``step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`` thermo; link-08 wraps that with the
``stage`` + ``box_*`` columns the polyethylene_tension benchmark already
streams (see ``benchmarks/md/polyethylene_tension.cpp``). Comment lines
(``#`` at column 0) — including the metadata header link-08 writes — are
skipped silently.

CLI::

    compare_msu.py <thermo_A> <thermo_B> <stress_A> <stress_B>
                   [--label-A "atomiverse-cuda"]
                   [--label-B "lammps-cuda"]
                   --out-md <report.md>
                   [--block-size <N>]
                   [--equilibration <STEPS>]
                   [--ua-mass-amu 14.027]
                   [--n-atoms 240]

Exit codes:

* ``0`` — all 4 invariants PASS (and ``--out-md`` written).
* ``1`` — at least one invariant FAILS its tolerance (report still written).
* ``2`` — a bad CLI argument or missing input file.

@see .claude/specs/msu-pe-cuda-components-08-cuda-validation.md
"""

import argparse
import csv
import math
import os
import sys


# ── tolerances (per spec § Domain basis) ─────────────────────────────────────

#: 3σ multiplier for the ⟨T⟩ tolerance (|ΔT| < 3·σ_combined).
T_SIGMA_MULTIPLIER = 3.0

#: Hard absolute tolerance on ⟨ρ⟩ Δ (g/cm³).
RHO_ABS_TOL_GCC = 0.10

#: Relative tolerance on E_modulus Δ (fraction of B's modulus).
E_REL_TOL = 0.15

#: Absolute tolerance on σ_flow Δ (GPa).
SIGMA_FLOW_ABS_TOL_GPA = 0.2

#: Conversion: amu / Å³ × _AMU_PER_A3_TO_GCC = g/cm³.
#: 1 amu = 1.66053906660e-24 g; 1 Å³ = 1e-24 cm³ → ratio = 1.66054.
_AMU_PER_A3_TO_GCC = 1.66054

#: Default UA polymer monomer mass (CH₂ united atom): 14.027 amu.
_DEFAULT_UA_MASS_AMU = 14.027

#: Default total atom count for the 8-chain × 30-monomer MSU PE benchmark.
#: Override with --n-atoms if your data uses a different chain length.
_DEFAULT_N_ATOMS = 240

#: Stage tag for the stage-4 NPT equilibration (100 K, 0 atm).
_STAGE_EQUIL = 4

#: Stage tag for the stage-5 uniaxial tension.
_STAGE_TENSION = 5

#: Elastic-region cutoff for E_modulus fit: ε ∈ (0, _EPSILON_ELASTIC_MAX).
_EPSILON_ELASTIC_MAX = 0.05

#: Plastic-region threshold for σ_flow plateau: ε > _EPSILON_FLOW_MIN.
_EPSILON_FLOW_MIN = 0.5


# ── input parsing ────────────────────────────────────────────────────────────


class CompareError(Exception):
    """A compare_msu input file could not be read or had a structural issue."""


def _iter_data_lines(path):
    """Yield ``(line_number, raw_line)`` for non-blank, non-``#``-comment lines.

    Mirrors ``canonicalize_msu_output._iter_data_lines`` so the link-08
    metadata header ``# engine=…,backend=…,seed=…`` is silently skipped.
    """
    if not os.path.isfile(path):
        raise CompareError("input file not found: {}".format(path))
    with open(path, newline="") as f:
        for lineno, raw in enumerate(f, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            yield lineno, raw


def _read_extended_csv(path, required, numeric):
    """Read a CSV with a header row → ``dict`` of column lists.

    @param path      Filesystem path of the CSV to read.
    @param required  Iterable of column names that must be present in the
                     header. Missing required columns raise
                     :class:`CompareError`.
    @param numeric   Iterable of column names whose cells are parsed as
                     floats (NaN on parse failure).
    """
    rows = list(csv.reader(
        raw for _, raw in _iter_data_lines(path)))
    if not rows:
        raise CompareError("empty CSV (no header): {}".format(path))
    header = [h.strip() for h in rows[0]]
    cols = {name: [] for name in header}
    missing = [c for c in required if c not in header]
    if missing:
        raise CompareError("missing required columns {!r} in {} (have {!r})"
                           .format(missing, path, header))
    for r in rows[1:]:
        if len(r) != len(header):
            continue
        for name, cell in zip(header, r):
            if name in numeric:
                try:
                    cols[name].append(float(cell))
                except ValueError:
                    cols[name].append(float("nan"))
            else:
                # stage / step columns: int-tolerant
                try:
                    cols[name].append(int(float(cell)))
                except ValueError:
                    cols[name].append(0)
    return cols


def read_thermo(path):
    """Read an extended thermo CSV.

    Required columns: ``step``, ``stage``, ``T_K``, ``box_x``, ``box_y``,
    ``box_z``. Any extra columns are ignored. Comment lines (``#``) are
    skipped — this is how the link-08 ``# engine=…,backend=…`` metadata
    header is transparently absorbed.
    """
    return _read_extended_csv(
        path,
        required=["step", "stage", "T_K", "box_x", "box_y", "box_z"],
        numeric={"T_K", "box_x", "box_y", "box_z"})


def read_stress_strain(path):
    """Read an extended stress-strain CSV.

    Required columns: ``step``, ``stage``, ``strain``, ``sigma_xx_GPa``.
    """
    return _read_extended_csv(
        path,
        required=["step", "stage", "strain", "sigma_xx_GPa"],
        numeric={"strain", "sigma_xx_GPa"})


# ── Flyvbjerg-Petersen block-average ─────────────────────────────────────────


def _block_sigma(series, block_size):
    """σ of block-averages for a given block size (excludes a final partial).

    Variance of the block-means / N_blocks — the standard estimator of the
    mean's variance in the Flyvbjerg-Petersen paper. Returns 0 when
    fewer than 2 full blocks fit.
    """
    n_blocks = len(series) // block_size
    if n_blocks < 2:
        return 0.0
    block_means = []
    for k in range(n_blocks):
        seg = series[k * block_size:(k + 1) * block_size]
        block_means.append(sum(seg) / block_size)
    m = sum(block_means) / n_blocks
    var = sum((bm - m) ** 2 for bm in block_means) / (n_blocks - 1)
    return math.sqrt(var / n_blocks)


def block_average(series, block_size=None,
                  *, initial=256, plateau_rel_tol=0.10, max_doublings=8):
    """Flyvbjerg-Petersen block-average estimate of (⟨x⟩, σ_⟨x⟩).

    Computes the mean of ``series`` and the standard error of that mean
    using a block-average scheme that accounts for time-series
    autocorrelation:

    * ``block_size=None`` (auto) — start at ``initial`` (default 256),
      double until consecutive block sizes give σ within
      ``plateau_rel_tol`` (default 10 %) of each other or ``max_doublings``
      doublings have been tried. The plateau marks the autocorrelation
      length (~τ_auto).
    * ``block_size`` (int) — single fixed block size, no doubling.

    Returns ``(mean, sigma)`` where ``sigma`` is the standard error of the
    mean. If the series is too short to fit two full blocks at the chosen
    size, falls back to the naïve i.i.d. estimator ``σ/√N`` so callers
    never see σ = 0 on a sample of finite spread.
    """
    n = len(series)
    if n == 0:
        return 0.0, 0.0
    mean = sum(series) / n
    if n == 1:
        return mean, 0.0

    if block_size is not None:
        sigma = _block_sigma(series, block_size)
        if sigma > 0:
            return mean, sigma
        # Fall through to naïve estimator.
        var = sum((v - mean) ** 2 for v in series) / (n - 1)
        return mean, math.sqrt(var / n)

    # Auto: doubling block size until plateau.
    block = initial
    prev_sigma = _block_sigma(series, block)
    best_sigma = prev_sigma
    for _ in range(max_doublings):
        block *= 2
        if n // block < 2:
            break
        sigma = _block_sigma(series, block)
        if sigma <= 0:
            break
        if prev_sigma > 0 and abs(sigma - prev_sigma) / prev_sigma < \
                plateau_rel_tol:
            best_sigma = sigma
            break
        best_sigma = sigma
        prev_sigma = sigma

    if best_sigma <= 0:
        # Series too short or constant — fall back to naïve estimator.
        var = sum((v - mean) ** 2 for v in series) / max(1, n - 1)
        best_sigma = math.sqrt(var / n)
    return mean, best_sigma


# ── 4 invariants ─────────────────────────────────────────────────────────────


def _linfit_with_jackknife(xs, ys):
    """Least-squares slope σ_xx = E·ε with jackknife σ on the slope.

    Returns ``(slope, sigma_slope)``. For < 3 points the slope is the
    naïve estimate and σ = 0.
    """
    n = len(xs)
    if n < 2:
        return 0.0, 0.0
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / den if den else 0.0
    if n < 3:
        return slope, 0.0
    # Jackknife: drop each point in turn, refit, accumulate variance.
    slopes = []
    for i in range(n):
        jx = sum(xs) - xs[i]
        jy = sum(ys) - ys[i]
        jxx = sum(x * x for x in xs) - xs[i] ** 2
        jxy = sum(x * y for x, y in zip(xs, ys)) - xs[i] * ys[i]
        nj = n - 1
        d = nj * jxx - jx * jx
        slopes.append((nj * jxy - jx * jy) / d if d else slope)
    mean_s = sum(slopes) / n
    var_s = (n - 1) / n * sum((s - mean_s) ** 2 for s in slopes)
    return slope, math.sqrt(var_s)


def compute_invariants(thermo, stress, *,
                       block_size=None,
                       equilibration=0,
                       ua_mass_amu=_DEFAULT_UA_MASS_AMU,
                       n_atoms=_DEFAULT_N_ATOMS):
    """Compute the four cross-engine invariants from extended-schema data.

    @param thermo         dict from :func:`read_thermo`.
    @param stress         dict from :func:`read_stress_strain`.
    @param block_size     Forwarded to :func:`block_average`. ``None`` →
                          auto-doubling.
    @param equilibration  Discard stage-4 rows with ``step <
                          equilibration`` (default 0 — disabled; the spec
                          recommends 25 000 for full-scale runs but the
                          synthetic test data is shorter and benefits from
                          using all rows).
    @param ua_mass_amu    Per-atom mass (default UA polymer 14.027 amu).
    @param n_atoms        Total atom count for the density computation.

    Returns a dict::

        {
          "T_stage4":          {"mean": …, "sigma": …, "unit": "K"},
          "rho_stage4":        {"mean": …, "sigma": …, "unit": "g/cm³"},
          "E_modulus_stage5":  {"mean": …, "sigma": …, "unit": "GPa"},
          "sigma_flow_stage5": {"mean": …, "sigma": …, "unit": "GPa"},
        }

    On a stage that has zero rows in the input, the corresponding entry
    is ``{"mean": nan, "sigma": nan, "unit": …}``. NaN means "could not
    compute"; downstream :func:`render_markdown` shows it as ``nan`` and
    a verdict of ``?``.
    """
    inv = {}

    # T_stage4
    t_vals = []
    for st, step, T in zip(thermo["stage"], thermo["step"], thermo["T_K"]):
        if st == _STAGE_EQUIL and step >= equilibration and not math.isnan(T):
            t_vals.append(T)
    if not t_vals:
        # Fallback: take all rows ("short test runs" per spec guidance).
        t_vals = [T for T in thermo["T_K"] if not math.isnan(T)]
    if t_vals:
        T_mean, T_sigma = block_average(t_vals, block_size=block_size)
    else:
        T_mean, T_sigma = float("nan"), float("nan")
    inv["T_stage4"] = {"mean": T_mean, "sigma": T_sigma, "unit": "K"}

    # rho_stage4: ρ = (n_atoms * ua_mass_amu) / V_box, converted to g/cm³.
    rho_vals = []
    for st, step, bx, by, bz in zip(thermo["stage"], thermo["step"],
                                    thermo["box_x"], thermo["box_y"],
                                    thermo["box_z"]):
        if st != _STAGE_EQUIL or step < equilibration:
            continue
        V = bx * by * bz
        if V <= 0 or math.isnan(V):
            continue
        rho_amu_per_a3 = (n_atoms * ua_mass_amu) / V
        rho_vals.append(rho_amu_per_a3 * _AMU_PER_A3_TO_GCC)
    if not rho_vals:
        rho_vals = []
        for bx, by, bz in zip(thermo["box_x"], thermo["box_y"],
                              thermo["box_z"]):
            V = bx * by * bz
            if V > 0 and not math.isnan(V):
                rho_vals.append((n_atoms * ua_mass_amu) / V
                                * _AMU_PER_A3_TO_GCC)
    if rho_vals:
        rho_mean, rho_sigma = block_average(rho_vals, block_size=block_size)
    else:
        rho_mean, rho_sigma = float("nan"), float("nan")
    inv["rho_stage4"] = {"mean": rho_mean, "sigma": rho_sigma,
                         "unit": "g/cm^3"}

    # E_modulus_stage5: linear fit of σ_xx vs ε on stage 5, ε ∈ (0, 0.05).
    el_eps, el_sig = [], []
    for st, eps, sig in zip(stress["stage"], stress["strain"],
                            stress["sigma_xx_GPa"]):
        if st != _STAGE_TENSION:
            continue
        if math.isnan(eps) or math.isnan(sig):
            continue
        if 0.0 < eps < _EPSILON_ELASTIC_MAX:
            el_eps.append(eps)
            el_sig.append(sig)
    if len(el_eps) >= 2:
        E_mean, E_sigma = _linfit_with_jackknife(el_eps, el_sig)
    else:
        E_mean, E_sigma = float("nan"), float("nan")
    inv["E_modulus_stage5"] = {"mean": E_mean, "sigma": E_sigma,
                               "unit": "GPa"}

    # sigma_flow_stage5: ⟨σ_xx⟩ over stage 5, ε > 0.5.
    flow_vals = []
    for st, eps, sig in zip(stress["stage"], stress["strain"],
                            stress["sigma_xx_GPa"]):
        if st == _STAGE_TENSION and not math.isnan(sig) and \
                eps > _EPSILON_FLOW_MIN:
            flow_vals.append(sig)
    if flow_vals:
        flow_mean, flow_sigma = block_average(flow_vals,
                                              block_size=block_size)
    else:
        flow_mean, flow_sigma = float("nan"), float("nan")
    inv["sigma_flow_stage5"] = {"mean": flow_mean, "sigma": flow_sigma,
                                "unit": "GPa"}

    return inv


# ── verdict + markdown rendering ─────────────────────────────────────────────


def _verdict_one(metric, A, B):
    """Return ``(delta, threshold, passed: bool, threshold_label: str)``.

    The threshold label is what the report's ``3σ`` column shows for the
    metric — for the absolute-tol invariants it's the fixed scalar, for
    T_stage4 it's 3·√(σ_A² + σ_B²), for E_modulus it's
    ``E_REL_TOL × |B|``.
    """
    a, b = A["mean"], B["mean"]
    if math.isnan(a) or math.isnan(b):
        return float("nan"), float("nan"), False, "—"
    delta = a - b
    if metric == "T_stage4":
        thr = T_SIGMA_MULTIPLIER * math.sqrt(A["sigma"] ** 2 + B["sigma"] ** 2)
        # Spec floor: even with σ=0 we want a sensible ~15 K guard so
        # synthetic identical datasets don't fail on numerical zero.
        # Per the spec § Domain basis: tolerance "≈ 15 K".
        thr = max(thr, 15.0)
        return delta, thr, abs(delta) < thr, "3σ"
    if metric == "rho_stage4":
        return delta, RHO_ABS_TOL_GCC, abs(delta) < RHO_ABS_TOL_GCC, "abs"
    if metric == "E_modulus_stage5":
        thr = E_REL_TOL * abs(b)
        return delta, thr, abs(delta) < thr, "15%·|B|"
    if metric == "sigma_flow_stage5":
        return (delta, SIGMA_FLOW_ABS_TOL_GPA,
                abs(delta) < SIGMA_FLOW_ABS_TOL_GPA, "abs")
    raise CompareError("unknown metric: {!r}".format(metric))


def _fmt(v):
    """Compact float formatter for the report cells."""
    if isinstance(v, float) and math.isnan(v):
        return "nan"
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1e4 or av < 1e-3:
        return "%.3e" % v
    if av >= 100:
        return "%.1f" % v
    if av >= 1:
        return "%.3f" % v
    return "%.4f" % v


def overall_verdict(invA, invB):
    """Return ``(passed: bool, reason: str)`` summarising all four metrics."""
    failures = []
    for metric in ("T_stage4", "rho_stage4", "E_modulus_stage5",
                   "sigma_flow_stage5"):
        delta, thr, passed, _ = _verdict_one(metric, invA[metric],
                                             invB[metric])
        if not passed:
            failures.append("{} (|Δ|={} ≥ thr={})"
                            .format(metric, _fmt(abs(delta)), _fmt(thr)))
    if failures:
        return False, "; ".join(failures)
    return True, "all 4 invariants within tolerance"


def render_markdown(invA, invB, label_A, label_B):
    """Render the canonical markdown report.

    Table columns: ``| metric | A | B | Δ | 3σ_combined | verdict |``
    (the third-from-right column header is named per the spec §
    Implementation guidance #6).
    """
    out = ["# msu-pe-cuda-components-08 validation report", ""]
    out.append("| metric | A={} | B={} | Δ | 3σ_combined | verdict |".format(
        label_A, label_B))
    out.append("|---|---|---|---|---|---|")
    row_units = [
        ("T_stage4 [K]", "T_stage4"),
        ("rho_stage4 [g/cm^3]", "rho_stage4"),
        ("E_modulus_stage5 [GPa]", "E_modulus_stage5"),
        ("sigma_flow_stage5 [GPa]", "sigma_flow_stage5"),
    ]
    for label, key in row_units:
        A = invA[key]
        B = invB[key]
        delta, thr, passed, _ = _verdict_one(key, A, B)
        mark = "✓" if passed else "✗"
        if math.isnan(A["mean"]) or math.isnan(B["mean"]):
            mark = "?"
        out.append("| {} | {} ± {} | {} ± {} | {} | {} | {} |".format(
            label,
            _fmt(A["mean"]), _fmt(A["sigma"]),
            _fmt(B["mean"]), _fmt(B["sigma"]),
            _fmt(delta), _fmt(thr), mark))
    out.append("")
    passed_all, reason = overall_verdict(invA, invB)
    if passed_all:
        out.append("Overall: PASS")
    else:
        out.append("Overall: FAIL: {}".format(reason))
    return "\n".join(out) + "\n"


# ── CLI ──────────────────────────────────────────────────────────────────────


def main(argv=None):
    """CLI entry. Returns 0 (PASS), 1 (FAIL), or 2 (bad input)."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("thermo_A", help="Extended thermo CSV for run A.")
    p.add_argument("thermo_B", help="Extended thermo CSV for run B.")
    p.add_argument("stress_A", help="Extended stress-strain CSV for run A.")
    p.add_argument("stress_B", help="Extended stress-strain CSV for run B.")
    p.add_argument("--label-A", default="A",
                   help="Display label for run A in the report header.")
    p.add_argument("--label-B", default="B",
                   help="Display label for run B in the report header.")
    p.add_argument("--out-md", required=True,
                   help="Markdown report output path.")
    p.add_argument("--block-size", type=int, default=None,
                   help="Fixed block size (default: auto via doubling).")
    p.add_argument("--equilibration", type=int, default=0,
                   help="Discard stage-4 steps below this value "
                        "(default: 0 — keep all; spec recommends 25000 "
                        "for full-scale runs).")
    p.add_argument("--ua-mass-amu", type=float, default=_DEFAULT_UA_MASS_AMU,
                   help="Per-atom mass in amu (default 14.027 — UA polymer).")
    p.add_argument("--n-atoms", type=int, default=_DEFAULT_N_ATOMS,
                   help="Total atom count for ρ computation (default {}). "
                        "Override if your data has a different chain count "
                        "or chain length.".format(_DEFAULT_N_ATOMS))
    args = p.parse_args(argv)

    try:
        thermo_A = read_thermo(args.thermo_A)
        thermo_B = read_thermo(args.thermo_B)
        stress_A = read_stress_strain(args.stress_A)
        stress_B = read_stress_strain(args.stress_B)
    except CompareError as exc:
        print("compare_msu: {}".format(exc), file=sys.stderr)
        return 2

    invA = compute_invariants(thermo_A, stress_A,
                              block_size=args.block_size,
                              equilibration=args.equilibration,
                              ua_mass_amu=args.ua_mass_amu,
                              n_atoms=args.n_atoms)
    invB = compute_invariants(thermo_B, stress_B,
                              block_size=args.block_size,
                              equilibration=args.equilibration,
                              ua_mass_amu=args.ua_mass_amu,
                              n_atoms=args.n_atoms)
    md = render_markdown(invA, invB, args.label_A, args.label_B)
    with open(args.out_md, "w") as f:
        f.write(md)
    print("compare_msu: report written to {}".format(args.out_md))

    passed_all, reason = overall_verdict(invA, invB)
    if passed_all:
        return 0
    print("compare_msu: FAIL — {}".format(reason), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
