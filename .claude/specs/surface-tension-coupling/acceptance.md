# acceptance — surface-tension-coupling

oracle: lammps
oracle_version: 30 Mar 2026 - Development (patch_30Mar2026-1074-g7f680de296)

## Build (binding)

MD-only cmake: `ATV_BUILD_MD=ON ATV_BUILD_KSPACE=ON ATV_BUILD_TESTS=ON
ATV_BUILD_REGRESSIONS=ON`; HF/PW/GPW/CUDA/Torch **OFF**.

## Criteria

| id | criterion | type | tol / note |
| -- | --------- | ---- | ---------- |
| A1 | Public `cpu::md::SurfaceTensionNPT` (or `NoseHooverNPT` + `set_surface_tension` / `set_couple_xy`) is an `Integrator`; no public `*Config` / `make_*` | api | compile + unit |
| A2 | With couple-xy, after MD steps \(L_x = L_y\) within \(10^{-9}\) relative (or absolute \(10^{-10}\) Å) | science | unit + regression |
| A3 | Offline pin: LAMMPS `fix npt … couple xy` with \(P_\perp \neq P_z\) expands compressed LJ FCC; pin documents \(L_x=L_y\) and \(L_z \neq L_x\) after 200 steps | pin | REFERENCE.md |
| A4 | Regression dual-gate: (i) \(L_x=L_y\); (ii) \(L_z\) moves differently from lateral under nonzero \(\gamma\) (or fixed \(P_\perp \neq P_n\)) and volume stays finite | science | e2e ≤ 60 s |
| A5 | Instantaneous \(\gamma_\mathrm{inst}\) formula matches definition from stress (finite; documented) | science | unit or e2e |
| A6 | Fair wall vs same LAMMPS twin; `ratio = T_atv/T_lammps ≤ 10` or `perf_exempt` with reason | perf | REFERENCE.md |
| A7 | No live third-party at ctest; pin path under tools/metalwalls_oracle and/or REFERENCE.md | process | offline only |

## Dual-gate (regression)

1. **Structural:** couple-xy lockstep \(L_x = L_y\).
2. **Surface-tension path:** anisotropic response (\(L_z\) vs \(L_x\)) consistent with \(P_\perp \neq P_n\) / nonzero \(\gamma\).

## Done when

All A1–A7 verified on `wf/surface-tension-coupling`; regression PASS; commit local only.
