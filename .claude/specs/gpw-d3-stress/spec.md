---
slug: gpw-d3-stress
status: approved
created: 2026-07-31
domain: gpw / dispersion / cell stress
queue: cp2k-method-loop
deferred: false
depends_on: []
unblocks: [gpw-d3-stress]
grilled: false
sibling: gpw-d3bj-stress (shared D3 dispersion stress fold — serialize)
---

# Spec: `gpw-d3-stress` — GPW analytic fixed-P cell stress with DFT-D3 zero damping

**Slug:** `gpw-d3-stress`
**Domain:** GPW (Quickstep) cell stress, geometry-only dispersion (CPU)
**Queue row:** `docs/developer/cp2k-method-queue.md` → `gpw-d3-stress`
**Oracle:** CP2K 2026.1, `TYPE DFTD3` (zero damping) + `STRESS_TENSOR NUMERICAL`,
offline pin via `tools/cp2k_oracle/` (never launched from ctest)

## Summary

`gpw_stress` (analytic fixed-\(P\) Hellmann–Feynman path) already folds the
D3-zero dispersion energy into the strained non-self-consistent Lagrangian on
both the Γ and multi-k paths (`src/cpu/gpw/internal/gpw_stress.cpp`, `cfg.d3`
branch, `kD3Cp2kPairCutoffBohr` pair cutoff). What is **missing** is any gate
on that fold:

- no dedicated offline CP2K pin (wave-15 census: commit `3ee28243` added D3
  stress on the **DFTB2** path only — the GPW pin is still absent);
- no unit test that `d3 = true` reaches the stress Lagrangian (D4 and MBD have
  `GpwStress.AnalyticMatchesFdWithD4` / `WithMbd`; D3-zero has none);
- no `regressions/gpw/gpw-d3-stress/` case.

This spec closes the capability with a **pin + dual gate**, following the
shipped `gpw-d4-stress` pattern. Implementation changes are expected to be
zero-to-small (gates may surface a divergence; fix if local).

## System (consistent with the D3-zero force twin)

Periodic water dimer — the same geometry family as
`tools/cp2k_oracle/inputs/gpw-d3-forces/water_dimer_d3_periodic_forces.inp`
and `regressions/gpw/gpw-d4-stress`:

| Item | Value |
|------|-------|
| Atoms | O H H O H H (6), O at origin |
| Cell | 12 × 12 × 12 Å, `PERIODIC XYZ` |
| Basis / PP | `DZVP-MOLOPT-SR-GTH`; O `GTH-PBE-q6`, H `GTH-PBE-q1` (`n_occ = 8`) |
| XC | PBE |
| Dispersion | DFT-D3 **zero damping**, PBE (\(s_6=1, s_8=0.722, s_{r,6}=1.217\)), `CALCULATE_C9_TERM .FALSE.` |
| Pair cutoff | `kD3Cp2kPairCutoffBohr` = 2 × `R_CUTOFF 20.0` Å (CP2K half-cutoff convention) |
| Cutoff | 200 Ry (ATV single grid; CP2K `CUTOFF 200, REL_CUTOFF 40, NGRIDS 4`) — matches `gpw-d4-stress` cost budget; the pinned quantities (E_disp, Δσ_disp) are geometry-only and cutoff-independent |

## Oracle pin (offline, difference of two runs)

D3 is geometry-only: it does not enter the Kohn–Sham potential, so the
converged density is identical with and without it. Therefore

\[
\sigma^{\mathrm{CP2K}}(\mathrm{D3}) - \sigma^{\mathrm{CP2K}}(\mathrm{plain})
  = \sigma_{\mathrm{disp}}^{\mathrm{D3}}
\]

with the (much larger, grid-noise-limited) electronic stress cancelling —
same isolation trick as the `gpw-d3-forces` force pin. Two CP2K
`STRESS_TENSOR NUMERICAL` runs on the identical geometry
(`tools/cp2k_oracle/inputs/gpw-d3-stress/water_dimer_d3_stress.inp` and
`water_dimer_plain_stress.inp`) give:

1. **E_disp** (hard pin) — CP2K `Dispersion energy` line, expected to
   reproduce the existing periodic D3-zero value ≈ −8.20894676620e−04 Ha;
2. **Δσ diagonal** (hard pin) — the D3-zero dispersion stress;
3. **σ(D3) total** (soft pin) — bounded by single-grid-vs-multigrid
   discretisation; tolerance justified in REFERENCE.md from measurement.

## Product path (public API only)

`cpu::gpw::GPWEvaluator{kinds}` with `d3 = true`, `compute_stress = true`,
`stress_use_fd = false` → Frame `kStress` [eV/Å³], `last_stress_ha_bohr3()`
[Ha/Bohr³], `last_result().e_disp` [Ha]. ATV's dispersion stress is isolated
the same way the pin is: analytic product with `d3 = true` minus analytic
product with `d3 = false` (warm-started; density identical).

## Dual gate

1. **Correctness:** E_disp pin (1e−8 Ha); Δσ_disp diagonal vs CP2K NUMERICAL
   difference (tol set from measured noise floor in REFERENCE.md); analytic
   vs multi-SCF `gpw_stress_fd` on the same `d3 = true` method (≤ 2e−5
   Ha/Bohr³, `gpw-d4-stress` precedent); off-diagonal ≈ 0; diag floor;
   `kStress` unit conversion.
2. **Speed:** T_atv (analytic product) / T_cp2k (NUMERICAL stress twin) ≤ 10,
   measured on the same host; else `perf_exempt` with reason.

## Deliverables

- `tools/cp2k_oracle/inputs/gpw-d3-stress/{water_dimer_d3_stress.inp, water_dimer_plain_stress.inp}`
- `tests/cpu/gpw/test_cpu_gpw_stress.cpp`: `AnalyticMatchesFdWithD3` +
  `D3ChangesGammaDiagonal` (cheap MT H₂ fixture, mirrors D4/MBD)
- `regressions/gpw/gpw-d3-stress/{REFERENCE.md, reg_gpw-d3-stress.cpp}` +
  CMake registration (`atv_add_regression(gpw gpw-d3-stress THREADS 8)`)
- `docs/tutorial/methods/gpw-d3-stress.md`
- queue row → done; roadmap dispersion-stress residual prose updated

## Out of scope

- D3(BJ) stress pin (`gpw-d3bj-stress` — sibling queue row, serialized)
- ATM/C9 stress (held behind ATM force pins per census)
- Shear / off-diagonal stress (orthorhombic design)
- Multi-k D3 stress pin (Γ pin here; the multi-k fold shares the same
  `e_disp` branch and stays gated by the Γ pin + unit FD identity)
