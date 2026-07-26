---
slug: gpw-uks-rpbe
status: approved
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [gpw-uks-revpbe, gpw-rpbe-water]
---

# Spec: `gpw-uks-rpbe` — Γ-UKS RPBE spin-GGA

**Slug:** `gpw-uks-rpbe`  
**Domain:** GPW / unrestricted KS  
**Queue:** CP2K method loop  
**Oracle:** offline CP2K Quickstep (never live in ctest)

## Summary

Productize **Γ-point collinear UKS** with **RPBE** (Hammer, Hansen, Nørskov
1999 — exponential exchange enhancement \(F_x = 1 + \kappa(1-e^{-\mu s^2/\kappa})\),
PBE \(\kappa/\mu\), **PBE correlation**; libxc `GGA_X_RPBE` + `GGA_C_PBE` /
CP2K `XC_FUNCTIONAL RPBE`). Distinct from **revPBE** (Zhang–Yang κ-revision of
the rational PBE form).

Today: `XcKind::RPBE` and closed-shell `xc_grid_rpbe` are product (RKS /
`gpw-rpbe-water`). Γ-UKS spin-GGA covers PADE / PBE / PBEsol / RevPBE / BLYP
via `xc_grid_*_spin`, but **RPBE is rejected at build/SCF**
(`"UKS supports … RevPBE, or BLYP"`) and there is no `xc_grid_rpbe_spin`
wrapper. This slug wires the thin RPBE spin path and pins O₂ triplet energy
against offline CP2K.

## Acceptance

1. `XcKind::RPBE` accepted for Γ-UKS (`GpwDft` / `GPWEvaluator` with
   `n_alpha`/`n_beta`).
2. `xc_grid_rpbe_spin` exists and collapses to closed-shell `xc_grid_rpbe`
   when \(\rho_\alpha=\rho_\beta=\rho/2\).
3. Dual KS SCF converges on centered O₂ (MT, DZVP-MOLOPT-GTH, GTH-PBE-q6,
   multiplicity 3, 400 Ry).
4. Offline CP2K pin for \(E_{\mathrm{tot}}\) and \(\langle S^2\rangle\) in
   `regressions/gpw/gpw-uks-rpbe/REFERENCE.md`.
5. Dual-gate: \(|E_{\mathrm{RPBE}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha on the same
   geometry.
6. Public Frame energy remains eV; internal Ha cited in REFERENCE.
7. Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or `perf_exempt` with reason.
8. No public `*Config` / factories; public API =
   `GPWEvaluator` + `single_point` / members only.

## Non-goals

- Multi-k UKS RPBE stress/forces depth (may wire for consistency with other
  GGAs, but regression is Γ energy only).
- Hybrid / meta-GGA, OT-specific RPBE pin, analytic forces regression for this
  slug.

## Implementation plan

1. Add `xc_grid_rpbe_spin` → `xc_grid_gga_spin(..., XCType::RPBE)`.
2. Extend UKS XcKind allow-lists and `xc_spin` dispatch in:
   - `gpw_solver.cpp` (Γ-UKS, k-UKS, `GpwDft::build`)
   - `gpw_forces.cpp` (Γ + k)
   - `gpw_stress.cpp` (Γ analytic UKS)
3. Unit tests in `test_cpu_gpw_uks.cpp` (collapse + build accept + toy converge).
4. Regression `reg_gpw-uks-rpbe` + CP2K input twin + REFERENCE.md.

## System (oracle twin)

| Item | Value |
|------|--------|
| Molecule | O₂, centered, R = 1.2075 Å |
| Cell | 12³ Å, `PERIODIC NONE`, MT |
| Basis / PP | DZVP-MOLOPT-GTH / GTH-PBE-q6 |
| Spin | LSD, multiplicity 3 → Nα=7, Nβ=5 |
| XC | RPBE |
| Cutoff | 400 Ry |

## Citations

- B. Hammer, L. B. Hansen, J. K. Nørskov, *Phys. Rev. B* **59**, 7413 (1999).
- G. Lippert, J. Hutter, M. Parrinello, *Mol. Phys.* **92**, 477 (1997) — GPW.
- G. J. Martyna, M. E. Tuckerman, *J. Chem. Phys.* **110**, 2810 (1999) — MT.
