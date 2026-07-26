# Spec — `gpw-hirshfeld-uks`

## Goal

Enable **Hirshfeld volume ratios** \(V/V_{\mathrm{free}}\) from the **UKS**
(Γ and multi-k) GPW total valence density \(\rho=\rho_\alpha+\rho_\beta\) and
density-scale free-atom **TS / MBD@rsSCS**. Closes the roadmap §8 residual
after RKS Γ + multi-k Hirshfeld ✓ (`gpw-mbd-hirshfeld`, `gpw-hirshfeld-kpoints`).

## Method

Same Hirshfeld weight / effective volume as `gpw-mbd-hirshfeld`
(`internal/hirshfeld_volumes.h`, exponential free-atom model). Hirshfeld is
**spin-independent**: the molecular density is the total electron density

\[
\rho(\mathbf r)=\rho_\alpha(\mathbf r)+\rho_\beta(\mathbf r)
\]

(Γ: collocate \(P_\alpha+P_\beta\); multi-k: \(\sum_k w_k(\rho_{k\alpha}+\rho_{k\beta})\)).
Volume ratios feed the existing `ts_pbe` / `mbd_rsscs_pbe` `vol_ratio`
argument (isolated or periodic MBD via `mbd_k_mesh` unchanged).

## Public API

```cpp
GpwDft gpw;
gpw.mbd = true;            // or ts
gpw.hirshfeld = true;
gpw.n_alpha = 7;           // UKS occupations
gpw.n_beta  = 5;
// optional multi-k:
// gpw.k_mesh[0] = gpw.k_mesh[1] = gpw.k_mesh[2] = 2;
```

- Requires `ts` **or** `mbd` (mutual exclusion with other disp flags unchanged)
- **Γ UKS** and **multi-k UKS** product paths
- RKS Γ / multi-k Hirshfeld paths unchanged
- Ratios on `GpwResult::volume_ratios`; energy on `GpwResult::e_disp` [Ha]
- No public `*Config`; no `make_*` / `build_*` factories

## Out of scope

- Analytic \(dv_A/dR\) volume forces
- All-electron free-atom densities / Hirshfeld-I
- Live CP2K / live libMBD in ctest
- vdW-DF / XDM / VV10 / D3/D4 Hirshfeld

## Offline oracle

CP2K 2026.1 stock has no TS/MBD `PAIR_POTENTIAL` (same as RKS Hirshfeld). Hard
pin is **ATV product** E_disp + volume ratios under fixed UKS settings,
cross-checked vs free-atom MBD (volume_ratio = 1) and the dual in-tree
Hirshfeld kernel. Notes under `tools/cp2k_oracle/inputs/gpw-hirshfeld-uks/`.

## Acceptance

See `acceptance.md`.
