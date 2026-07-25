# Spec — `gpw-uks-pbe`

## Goal

Ship **Γ-point UKS with PBE-GGA** on the GPW path. Today Γ-UKS is PADE-LSD only
(`reg_gpw-uks-o2-gamma`); RKS already has closed-shell PBE/PBEsol
(`xc_grid_gga`). This slug closes the `gpw_solver.h` follow-on:

> UKS PBE-GGA … are follow-on.

## Scope

| In | Out |
|----|-----|
| Spin-resolved GGA XC on the GPW grid (`xc_grid_gga_spin` / PBE) | k-point UKS (`gpw-kpoints-uks`) |
| Wire UKS SCF loop for `XcKind::PBE` (and PBEsol if free) | Analytic UKS forces / stress |
| Offline CP2K pin O₂ triplet UKS-PBE, MT / MOLOPT | OT + UKS |
| Regression `regressions/gpw/gpw-uks-pbe` | Multigrid product collocation for UKS |

## Science

Unrestricted KS with PBE (Perdew–Burke–Ernzerhof) GGA:

\[
v_\sigma = v_{\rho_\sigma} - \nabla\cdot\bigl(2\,v_{\sigma_{\sigma\sigma}}\nabla\rho_\sigma
  + v_{\sigma_{\alpha\beta}}\nabla\rho_{\bar\sigma}\bigr)
\]

Densities \(\rho_\alpha,\rho_\beta\) from \(P_\alpha,P_\beta\); spectral
gradient/divergence (same as RKS GGA). Hartree from total \(\rho=\rho_\alpha+\rho_\beta\).

## API (no public Config / factories)

```cpp
cpu::gpw::GPWEvaluator eval({O});
eval.xc = cpu::gpw::XcKind::PBE;
eval.n_alpha = 7;
eval.n_beta = 5;
eval.nonperiodic = true;
cpu::single_point(frame, eval);
// energy on frame.meta [eV]; last_energy_hartree / last_s2 for pins
```

## Acceptance

1. Unit: spin-PBE collapses to closed-shell PBE when \(\rho_\alpha=\rho_\beta=\rho/2\);
   UKS+PBE accepted at `build()`; tiny high-spin H₂ converges with PBE.
2. Regression vs offline CP2K 2026.1 O₂ UKS-PBE pin: \(|\Delta E|\le 1\,\mathrm{mHa}\),
   \(|\Delta\langle S^2\rangle|\le 0.15\); public energy still eV on Frame.
3. Perf: \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) or documented `perf_exempt`.

## Twin system

Same geometry as `gpw-uks-o2-gamma`: O₂ R=1.2075 Å, 12 Å cubic, MT / `PERIODIC NONE`,
DZVP-MOLOPT-GTH + GTH-PBE-q6, cutoff 400 Ry, multiplicity 3 → Nα=7 Nβ=5,
`XC_FUNCTIONAL PBE`.
