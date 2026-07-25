# Acceptance: `gpw-variable-cell-ice`

Each criterion is observable + tolerance vs the offline CP2K pin (numbers
hard-coded in `regressions/gpw/gpw-variable-cell-ice/REFERENCE.md` and the
regression binary; CP2K is never launched from any registered test).

## A. Capability (unit tier, `tests/`)

- [x] A1 `AnisotropicCellRelax` with `relax_y = relax_z = false` leaves the
      y/z box edges bit-identical over a Minimize run and changes only x.
- [x] A2 `AnisotropicCellRelax` with `stress_is_virial = false` interprets
      `kStress` as intensive [eV/Å³]: with a synthetic constant stress σ_xx
      the first scale step is `1 + clamp(alpha0·σ_xx/bulk_mod, ±max_scale)`
      (no volume division), asserted to 1e-12.
- [x] A3 `StressConverged{0, stol, x-only}` verdict: Converged when
      |σ_xx| < stol regardless of masked σ_yy/σ_zz; Continue otherwise;
      name() == "stress_converged".
- [x] A4 `gpw_stress` axes mask: masked axes return exactly 0 and are not
      FD-solved; unmasked axis on the existing tiny-system stress unit tier
      matches the full-axes call to 1e-12 (warm start changes cost, not
      values beyond SCF tol).
- [x] A5 Variable cell: `GpwDft` cell-member change between runs is
      well-defined (unit `GpwMethod.VariableCellRunRebuildsFromCellMember`);
      `GPWEvaluator` box-change launch no longer throws — exercised by the
      regression relax loop (frame box rescaled every iteration).

## B. Physics vs CP2K pin (regression, `regressions/gpw/gpw-variable-cell-ice/`)

Twin: periodic water chain, **PADE LDA / O SZV-GTH-PADE + H DZV-GTH-PADE /
GTH-PADE** (measured basis pivot — periodic MOLOPT image sums cost > 10 min
per relax launch, see spec.md and REFERENCE.md), 300 Ry, fixed fractional
coordinates, relax a only from a₀ = 3.34 Å (b = c = 8 Å fixed).

- [x] B1 Starting from a₀ = a\*(CP2K) × (1 + δ) with δ = 4%, the Minimize +
      AnisotropicCellRelax + GPWEvaluator loop **converges** (criterion
      `stress_converged`, not budget exhaustion). Verified 2026-07-25.
- [x] B2 Relaxed chain lattice constant: |a_ATV − a\*_CP2K(CELL_OPT)| ≤
      **0.02 Å** (CP2K-internal consensus spread 0.008 Å + strain-FD h +
      grid differences). **Measured 1.0×10⁻⁴ Å** (3.211249 vs 3.211148).
- [x] B3 Residual stress at the relaxed cell: |σ_xx| ≤ 1.2×stol,
      stol = 5×10⁻⁴ eV/Å³ (= 2.7×10⁻⁶ Ha/Bohr³ ≤ the 2×10⁻⁵ Ha/Bohr³
      ceiling). **Measured 4.4×10⁻⁷ eV/Å³.**
- [x] B4 Energy at the relaxed cell within **1×10⁻³ Ha** of the CP2K CELL_OPT
      final energy. **Measured 3.5×10⁻⁵ Ha.**
- [x] B5 Masked axes: final b and c exactly 8 Å (bit-identical). Verified.
- [x] B6 Regression wall ≤ 60 s (single rank, 4 OMP): **5.1 s**; public
      headers only; hard-coded pins; registered via
      `atv_add_regression(gpw gpw-variable-cell-ice)`.

## C. Process

- [x] C1 REFERENCE.md cites CP2K 2026.1 build string, input echo pointer
      (`tools/cp2k_oracle/inputs/gpw-variable-cell-ice/`), scan table,
      stress pins, CELL_OPT cross-check, and walls.
- [x] C2 Scoped gate green (2026-07-25: 3/3 ctest in 2.4 s; run_reg PASS 5.1 s): `--target atomiverse`, touched test binaries,
      `reg_gpw-variable-cell-ice`; ctest -R exactly the touched test names;
      regression via `--target run_reg_gpw-variable-cell-ice`.
- [x] C3 (ratio 1.29 — CELL_OPT twin 3.9 s vs ATV relax 5.1 s) Perf: fair ratio T_atv/T_cp2k recorded in REFERENCE.md;
      within_10x true or perf_exempt with reason.
- [x] C4 Discovered debt (kStress convention divergence; `copy_gpw` dropped
      fields) named in the final summary; `copy_gpw` fixed in-slug.
