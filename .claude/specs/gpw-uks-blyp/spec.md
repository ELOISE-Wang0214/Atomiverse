# Spec — `gpw-uks-blyp`

## Goal

Ship **Γ-point UKS with BLYP spin-GGA** on the GPW path. RKS already has
closed-shell BLYP (`xc_grid_blyp` / `reg_gpw-blyp-water`); Γ-UKS spin-GGA
today supports only PBE/PBEsol (`xc_grid_gga_spin` wrappers + UKS gates).
This slug closes the residual:

> UKS spin-GGA wrappers only PBE/PBEsol; extend spin path + solver for BLYP.

## Scope

| In | Out |
|----|-----|
| `xc_grid_blyp_spin` thin wrapper over `xc_grid_gga_spin(BLYP)` | New XC functional math (reuse `evaluate_xc`) |
| Wire Γ-UKS SCF / forces / stress (+ k-UKS dispatch) for `XcKind::BLYP` | OT + UKS, multigrid product UKS |
| Offline CP2K pin O₂ triplet UKS-BLYP, MT / MOLOPT | Molecular post-HF, PW spin-GGA |
| Regression `regressions/gpw/gpw-uks-blyp` | Analytic force re-tune beyond energy pin |

## Science

Unrestricted KS with BLYP (Becke 1988 exchange + Lee–Yang–Parr correlation;
libxc `GGA_X_B88` + `GGA_C_LYP` / CP2K `XC_FUNCTIONAL BLYP`):

\[
v_\sigma = v_{\rho_\sigma} - \nabla\cdot\bigl(2\,v_{\sigma_{\sigma\sigma}}\nabla\rho_\sigma
  + v_{\sigma_{\alpha\beta}}\nabla\rho_{\bar\sigma}\bigr)
\]

Densities \(\rho_\alpha,\rho_\beta\) from \(P_\alpha,P_\beta\); spectral
gradient/divergence (same engine as RKS GGA / UKS-PBE). Hartree from total
\(\rho=\rho_\alpha+\rho_\beta\).

## API (no public Config / factories)

```cpp
cpu::gpw::GPWEvaluator eval({O});
eval.xc = cpu::gpw::XcKind::BLYP;
eval.n_alpha = 7;
eval.n_beta = 5;
eval.nonperiodic = true;
cpu::single_point(frame, eval);
// energy on frame.meta [eV]; last_energy_hartree / last_s2 for pins
```

## Acceptance

1. Unit: spin-BLYP collapses to closed-shell BLYP when \(\rho_\alpha=\rho_\beta=\rho/2\);
   UKS+BLYP accepted at `build()`; tiny high-spin H₂ converges with BLYP.
2. Regression vs offline CP2K 2026.1 O₂ UKS-BLYP pin: \(|\Delta E|\le 1\,\mathrm{mHa}\),
   \(|\Delta\langle S^2\rangle|\le 0.15\); public energy still eV on Frame.
3. Dual-gate: \(|E_{\mathrm{BLYP}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha on same system
   (activation; PBE twin from `gpw-uks-pbe` setup).
4. Perf: \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) or documented `perf_exempt`.

## Twin system

Same geometry as `gpw-uks-pbe` / `gpw-uks-o2-gamma`: O₂ R=1.2075 Å, 12 Å cubic,
MT / `PERIODIC NONE`, DZVP-MOLOPT-GTH + GTH-PBE-q6, cutoff 400 Ry, multiplicity 3
→ Nα=7 Nβ=5, `XC_FUNCTIONAL BLYP`.
