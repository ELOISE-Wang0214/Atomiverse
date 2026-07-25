# Acceptance — periodic-rpa-kpoints

| Gate | Criterion | How |
|------|-----------|-----|
| A1 | multi-k `direct_rpa` runs (no nk>1 throw) | `reg_periodic-rpa-kpoints` |
| A2 | `e_corr < 0`, finite; `n_kpts==2`; `ri`; Ω>0; e_total identity | regression asserts |
| A3 | N_k=1 multi-k ≡ real plasmon | `test_cpu_gpw_rpa` KpointsNk1RecoversRealPlasmon |
| A4 | first-green \(E_c\) band (tol 1e-4 Ha setup) | REFERENCE.md + reg pin |
| A5 | wall ≤ 60 s | regression |
| A6 | perf ≤10× or perf_exempt | no CP2K on pin host → exempt |

**acceptance_met** when A1–A5 green on worker; A6 exempt documented.
