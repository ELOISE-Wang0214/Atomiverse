# Spec — `gpw-hirshfeld-kpoints`

## Goal

Enable **Hirshfeld volume ratios** \(V/V_{\mathrm{free}}\) from the **multi-k
RKS** GPW valence density and density-scale free-atom **TS / MBD@rsSCS**
parameters. Closes the roadmap §8 residual after Γ-only Hirshfeld ✓
(`gpw-mbd-hirshfeld`).

## Method

Same Hirshfeld weight / effective volume as `gpw-mbd-hirshfeld`
(`internal/hirshfeld_volumes.h`, exponential free-atom model). On the multi-k
path the molecular density is the **k-summed** real density

\[
\rho(\mathbf r)=\sum_{\mathbf k} w_{\mathbf k}\,\rho_{\mathbf k}(\mathbf r)
\]

collocated from the final Bloch \(P_k\) (same grid as the SCF Poisson/XC).
Volume ratios then feed the existing `ts_pbe` / `mbd_rsscs_pbe` `vol_ratio`
argument (isolated or periodic MBD via `mbd_k_mesh` unchanged).

## Public API

```cpp
GpwDft gpw;
gpw.mbd = true;            // or ts
gpw.hirshfeld = true;
gpw.k_mesh[0] = gpw.k_mesh[1] = gpw.k_mesh[2] = 2; // multi-k RKS
// optional: gpw.mbd_k_mesh = {2,2,2} for periodic Ewald MBD
```

- Requires `ts` **or** `mbd` (mutual exclusion with other disp flags unchanged)
- **Multi-k RKS** product path (Monkhorst–Pack or explicit `k_points`)
- Γ path remains the existing `fold_dispersion` (no behavior change)
- **UKS still throws** (follow-on; out of scope this slug)
- Ratios on `GpwResult::volume_ratios`; energy on `GpwResult::e_disp` [Ha]
- No public `*Config`; no `make_*` / `build_*` factories

## Out of scope

- UKS Hirshfeld (Γ or multi-k)
- Analytic \(dv_A/dR\) volume forces
- All-electron free-atom densities / Hirshfeld-I
- Live CP2K / live libMBD in ctest
- vdW-DF / XDM / VV10

## Offline oracle

CP2K 2026.1 stock has no TS/MBD `PAIR_POTENTIAL` (same as Γ Hirshfeld). Hard
pin is **ATV product** E_disp + volume ratios under fixed multi-k settings,
cross-checked vs free-atom multi-k MBD (volume_ratio = 1) and the dual
in-tree Hirshfeld kernel. Notes under
`tools/cp2k_oracle/inputs/gpw-hirshfeld-kpoints/`.

## Acceptance

See `acceptance.md`.
