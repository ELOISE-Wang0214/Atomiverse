# Spec — `gpw-optb88-vdw`

## Goal

Ship **optB88-vdW** on Γ-point GPW: Klimeš–Bowler–Michaelides optimized
Becke-88 exchange + PW92 LDA correlation + nonlocal **DRSLL** (vdW-DF)
correlation, aligned with the CP2K twin
`XC_GGA_X_OPTB88_VDW` + `PW92` + `NON_LOCAL TYPE DRSLL`.

## Science

| Piece | Form | Source |
| ----- | ---- | ------ |
| Exchange | B88 enhancement with \(\beta=0.003368671161\), \(\gamma=6.981317007977318\) | libxc `GGA_X_OPTB88_VDW` (id 139); Klimeš et al., *J. Phys.: Condens. Matter* **22**, 022201 (2010) |
| Local correlation | PW92 LDA only (no GGA-C) | CP2K Hutter recipe |
| Nonlocal correlation | DRSLL / DF1 post-SCF \(E_c^{\mathrm{nl}}\) | existing `GpwDft::vdw_df` / `vdw_df_drsll` |

Full functional energy:

\[
E_{\mathrm{optB88\text{-}vdW}}
  = E_{\mathrm{optB88+x}} + E_{\mathrm{PW92}} + E_c^{\mathrm{nl}}(\mathrm{DRSLL})
  + E_{\mathrm{KS\ other}}
\]

Nonlocal term is **post-SCF** on the final density (same product depth as
`gpw-vdw-df-ice`); SCF potential / forces remain `gpw-vdw-scf-potential`.

## Public API (no Config bags, no factories)

```cpp
cpu::gpw::GpwDft gpw;
gpw.xc = cpu::gpw::XcKind::OptB88VdW; // semi-local optB88-x + PW92-c
gpw.vdw_df = true;                    // nonlocal DRSLL
gpw.build(atoms);
const auto& res = gpw.run(atoms);
// res.e_disp  [Ha]  — E_c^{nl}
// res.e_total [Ha]  — full KS + nonlocal
// Frame kEnergy [eV] via GPWEvaluator
```

`XCType::OptB88VdW` is the shared pointwise engine for
`evaluate_xc` / GPW grid XC.

## Acceptance

1. Unit: `evaluate_xc(OptB88VdW)` matches libxc exchange+PW92 composition
   (or B88@opt params + existing PW92) to **1e−8** on representative points;
   differs from `BLYP` and `PBE` at finite gradient.
2. Unit/grid: `xc_grid_optb88_vdw` constant-density limit matches
   `evaluate_xc` (σ=0).
3. Regression `reg_gpw-optb88-vdw` (water dimer, same geometry as
   `gpw-vdw-df-ice`):
   - SCF converges; finite `e_total`, `e_disp`.
   - **Dual gate A:** \(|E_{\mathrm{disp}}^{\mathrm{optB88}} - E_{\mathrm{disp}}^{\mathrm{PBE+vdW\text{-}DF}}|\) may be small (density-driven); primary
     pin is absolute \(E_{\mathrm{disp}}\) and/or
     \(|E_{\mathrm{tot}}^{\mathrm{optB88}} - E_{\mathrm{tot}}^{\mathrm{PBE}}| > 10^{-4}\) Ha
     (semi-local activation).
   - **Dual gate B:** offline CP2K pin when available; else ATV first-green
     freeze + libxc exchange pin (same policy as `gpw-blyp-water` /
     `gpw-vdw-df-ice`).
4. Perf: **perf_exempt** (E_disp / XC activation deliverable; same family as
   other nonlocal vdW regs).

## Out of scope

- Self-consistent nonlocal vdW potential (`gpw-vdw-scf-potential`)
- UKS product regression (RKS Γ only this slug; UKS may accept the XcKind)
- VV10 stock, optPBE-vdW, optB86b-vdW
- Plane-wave path, molecular Gaussian-only path regressions
