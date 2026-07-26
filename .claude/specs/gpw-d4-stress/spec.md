---
slug: gpw-d4-stress
status: acceptance_met
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [gpw-d4-water, gpw-analytic-stress]
---

# Spec: `gpw-d4-stress` — GPW analytic cell stress with DFT-D4

**Slug:** `gpw-d4-stress`  
**Domain:** GPW / periodic DFT (CPU)  
**Queue:** `cp2k-method-loop`  
**Depends on:** `gpw-d4-water` ✓ (E_disp product), `gpw-analytic-stress` ✓ (fixed-P stress)

## Summary

Productize **Hellmann–Feynman fixed-\(P\) analytic cell stress** when
`GpwDft::d4` / `GPWEvaluator::d4` is set (Grimme/Caldeweyher DFT-D4 BJ,
two-body, PBE parameters). Today `gpw_stress` / `energy_fixed_P` **throws**
on D4 (Γ and multi-k). Energy already folds D4 into \(E_{\mathrm{tot}}\) /
\(e_{\mathrm{disp}}\) (`gpw-d4-water`).

## Domain basis

D4 is a **geometry-only** additive London term (same class as free-atom D3
in the fixed-\(P\) path). Under orthorhombic diagonal strain at fixed
fractional coordinates, rebuild positions + cell and re-evaluate:

\[
E_{\mathrm{disp}}^{\mathrm{D4}}(\{\mathbf{R}\},\mathbf{h})
  = \mathrm{dftd4\_energy}(Z, \mathbf{R}[\varepsilon], \mathbf{h}[\varepsilon])
\]

and include it in the fixed-\(P\) Lagrangian \(E_L[P;\varepsilon]\)
alongside \(T\), \(V_{\mathrm{pp}}\), \(E_{\mathrm{es}}\), \(E_{\mathrm{xc}}\),
and \(-\mathrm{Tr}(WS)\). Then

\[
\sigma_{aa} = -\Omega^{-1}\partial E_L/\partial\varepsilon_{aa}
\quad\text{(CP2K STRESS\_TENSOR sign; Ha/Bohr³)}.
\]

Oracle twin: multi-SCF strain-FD (`gpw_stress_fd`) on the same `d4=true`
method (re-converged SCF includes D4 automatically via `fold_dispersion`).

ATM three-body remains **off** (`s9=0`), matching `gpw-d4-water`.

## Acceptance

1. `gpw_stress` / analytic evaluator stress **does not throw** when `d4=true`
   (Γ RKS product path at minimum; multi-k same fold).
2. Dual-gate: \(\max_a|\sigma_a^{\mathrm{an}}-\sigma_a^{\mathrm{fd}}|
   \le\) stated tol (Ha/Bohr³) on the regression twin.
3. Soft energy: \(e_{\mathrm{disp}}\) within 1e-8 Ha of the offline dftd4
   3.7.0 pin from `gpw-d4-water` on the same geometry.
4. Off-diagonal product stress remains ~0 (orthorhombic).
5. Unit test: cheap fixture analytic vs FD with `d4=true` does not throw and
   residual is unit-tier loose.
6. Regression under `regressions/gpw/gpw-d4-stress/` with `REFERENCE.md`.
7. Dual gate green; `perf_exempt` allowed (Arrhenius CP2K lacks DFTD4 —
   same policy as `gpw-d4-water`).

## Out of scope

- D4 analytic nuclear **forces** (separate slug if queued).
- ATM three-body / MBD-scale terms.
- Shear / off-diagonal stress.
- Live CP2K or live dftd4 in ctest.

## Files (expected)

- `src/cpu/gpw/internal/gpw_stress.cpp` — remove D4 throw; fold `dftd4_energy`
- `include/atomiverse/cpu/gpw/internal/gpw_stress.h` — doc note
- `tests/cpu/gpw/test_cpu_gpw_stress.cpp` — D4 dual-gate unit
- `regressions/gpw/gpw-d4-stress/` — e2e dual-gate + E_disp pin
- `tools/cp2k_oracle/inputs/gpw-d4-stress/` — optional twin note (reuse d4-water)
