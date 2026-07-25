# Acceptance — `gpw-kpoints-uks`

## Unit (`test_cpu_gpw_uks`)

| Gate | Check | Tol |
| ---- | ----- | --- |
| Build | `k_mesh=2,2,2` + UKS accepted | no throw |
| A1 Γ-limit | explicit k=Γ k-UKS vs real Γ-UKS (H₂ Na=2 Nb=0) | ≤ 1e-6 Ha |
| A2 closed-shell | multi-k UKS Na=Nb=1 vs multi-k RKS n_occ=1 | ≤ 5e-6 Ha |
| A3 time-reversal | MP 2×2×2 full 8 k vs TR 4 k (UKS) | ≤ 5e-6 Ha |

## Regression (`reg_gpw-kpoints-uks`)

| Gate | Check |
| ---- | ----- |
| Converge | O₂ periodic MP 2×2×2 UKS-PBE converges |
| Activation | \|E_MP − E_Γ\| ≥ 1e-4 Ha (k-sampling moves energy) |
| Magnetization | fixed occupations ⇒ M = Nα − Nβ = 2 |
| Absolute | e_tot vs REFERENCE pin (CP2K or frozen first-green) |
| No live oracle | ctest does not spawn CP2K |

## Perf

`T_atv / T_cp2k ≤ 10` when a fair twin wall is available; else `perf_exempt`
with reason.
