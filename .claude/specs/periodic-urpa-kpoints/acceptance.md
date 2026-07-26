# Acceptance — `periodic-urpa-kpoints`

| # | Check | Met |
|---|--------|-----|
| 1 | Spec under `.claude/specs/periodic-urpa-kpoints/` | yes |
| 2 | Offline CP2K input `tools/cp2k_oracle/inputs/periodic-urpa-kpoints.inp` | yes |
| 3 | Kernel `direct_rpa_uhf_energy_kpoints`; public `GpwDft::direct_rpa` multi-k UKS no throw | yes |
| 4 | Unit: \(N_k=1\) α=β multi-k URPA ≡ real plasmon / Γ URPA | yes |
| 5 | Regression dual-gate vs RKS multi-k \(E_c\) band; hard multi-k gates | yes (bit-match −0.0736273682 Ha) |
| 6 | Scoped ctest + `run_reg_periodic-urpa-kpoints` green; wall ≤ 60 s | yes (~2.7 s) |
| 7 | `perf_exempt` documented (no fair CP2K plasmon twin) | yes |
| 8 | Commit on `wf/periodic-urpa-kpoints`; build wiped | pending commit |

`acceptance_met=true` when 1–8 pass after commit+wipe.
