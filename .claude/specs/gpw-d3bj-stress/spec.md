---
slug: gpw-d3bj-stress
status: code-complete (all 6 acceptance criteria verified 2026-07-31; close on merge)
created: 2026-07-31
revised: 2026-07-31
domain: gpw / cell stress / dispersion
queue: cp2k-method-loop
deferred: false
depends_on: []
unblocks: [gpw-d3bj-stress]
grilled: false
sibling: gpw-d3-stress (shared D3 stress fold — patches must stay disjoint)
---

# Spec: `gpw-d3bj-stress` — GPW analytic fixed-P cell stress with DFT-D3(BJ)

**Slug:** `gpw-d3bj-stress`
**Domain:** GPW (CPU), analytic fixed-P cell stress, DFT-D3(BJ) dispersion fold
**Queue row:** `gpw-d3bj-stress` (GapCensus wave 13, re-verified wave 15)
**Oracle:** CP2K 2026.1 `TYPE DFTD3(BJ)` + `STRESS_TENSOR NUMERICAL`, offline
via `tools/cp2k_oracle` — never at ctest time.

## Gap

The Γ / multi-k analytic fixed-P stress (`gpw_stress` →
`internal::gpw_stress_analytic_impl`) already **folds** the geometry-only
D3(BJ) dispersion energy — `dftd3_pbe_bj` is rebuilt at every strained cell in
`energy_fixed_P` / `energy_fixed_P_k` (`src/cpu/gpw/internal/gpw_stress.cpp`).
But the fold has **no dedicated offline oracle pin and no e2e regression**:
the D4 twin ships `reg_gpw-d4-stress`, the D3(BJ) *forces* ship
`reg_gpw-d3bj-forces`, while the D3(BJ) *stress* is only reachable through the
generic strain-FD twin. Foothold ≠ capability closed.

## Deliverable

1. **Offline CP2K pin** (new inputs under
   `tools/cp2k_oracle/inputs/gpw-d3bj-stress/`): two CP2K 2026.1 runs on the
   shipped `gpw-d3bj-forces` twin geometry (periodic water dimer, 12 Å cubic
   cell, PBE, DZVP-MOLOPT-SR-GTH, GTH-PBE, CUTOFF 300) with
   `STRESS_TENSOR NUMERICAL`:
   - `water_dimer_d3bj_stress.inp` — PBE-D3(BJ), CP2K PBE BJ parameters,
     `CALCULATE_C9_TERM .FALSE.`, `R_CUTOFF 20.0`;
   - `water_dimer_plain_stress.inp` — identical, no `&VDW_POTENTIAL`.

   D3 is geometry-only (never enters the KS potential), so for every strained
   cell the two SCF densities are identical and

   \[ \sigma^{\mathrm{CP2K}}(\mathrm{D3BJ}) - \sigma^{\mathrm{CP2K}}(\mathrm{plain})
      = \sigma_{\mathrm{disp}}^{\mathrm{D3BJ}} \]

   up to shared FD truncation — the large electronic stress cancels exactly.
   The pin is this **difference of two numbers CP2K printed** plus the CP2K
   "Dispersion energy" line (re-runnable: binary + version + input + line).

2. **Regression** `regressions/gpw/gpw-d3bj-stress/` (REFERENCE.md +
   `reg_gpw-d3bj-stress.cpp`, registered `atv_add_regression(gpw
   gpw-d3bj-stress THREADS 8)`), out-of-library consumer, gates:
   - **E_disp pin:** standalone `dftd3_pbe_bj` on the twin nuclei vs the CP2K
     printed dispersion energy (cutoff-independent).
   - **Dispersion-stress pin:** in-regression strain-FD of the geometry-only
     `dftd3_pbe_bj` energy (scaled positions + cell, same convention as the
     fold) vs the CP2K numerical-stress **difference**, diagonal components.
   - **Dual gate (decisive):** product analytic `GPWEvaluator{d3bj=true,
     compute_stress=true, stress_use_fd=false}` vs multi-SCF strain-FD
     `gpw_stress_fd` on the same `d3bj=true` method.
   - **Wiring:** product analytic σ(d3bj) minus plain analytic σ equals the
     standalone dispersion stress (same code both sides) — only if wall budget
     allows a second SCF; else the dual gate plus the `d3bj`-moves-the-diagonal
     guard covers it.
   - Frame `kStress` [eV/Å³] unit conversion; off-diagonal ≈ 0; diag floor.

3. **Unit test:** `GpwStress.AnalyticMatchesFdWithD3Bj` on cheap MT H₂
   (mirrors `AnalyticMatchesFdWithD4`) in
   `tests/cpu/gpw/test_cpu_gpw_stress.cpp` — fast tier, no oracle.

## Non-goals

- No change to the D3(BJ) energy/CN implementation (pinned by
  `gpw-d3bj-forces` / `gpw-pbe-d3bj-ice`).
- No ATM/C9 term (that is `gpw-d3-atm-forces` / D3 ATM rows).
- No k-path stress work (fold already present multi-k; Γ pin here).
- Zero-damping D3 stress is the sibling `gpw-d3-stress` — do not touch its
  surface (`dftd3_pbe`), keep patches disjoint.

## Tolerances (spec-set, finalized in REFERENCE.md from measured floors)

| Gate | Target |
| --- | --- |
| \|E_disp − pin\| | ≤ 1e−8 Ha |
| diag \|σ_disp^ATV − Δσ^CP2K\| | ≤ 5e−9 Ha/Bohr³ (signal ~1e−7; CP2K FD-difference noise floor measured before pinning) |
| dual gate max diag \|σ_an − σ_fd\| | ≤ 2e−5 Ha/Bohr³ (D4 twin measured 7.4e−7) |
| off-diag product | ≤ 1e−8 Ha/Bohr³ |
| wall | ≤ 60 s target; perf ratio vs CP2K stress run ≤ 10 or perf_exempt |
