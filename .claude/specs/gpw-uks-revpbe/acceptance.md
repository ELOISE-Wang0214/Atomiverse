# Acceptance — gpw-uks-revpbe

| Check | Criterion |
| --- | --- |
| Unit | `ctest -R '^(test_cpu_gpw_uks)$'` passes (RevPBE spin collapse, build, toy SCF) |
| Regression | `run_reg_gpw-uks-revpbe` PASS; \|E − E_CP2K\| ≤ 1e-3 Ha; \|S² − pin\| ≤ 0.15 |
| Dual-gate | \|E_revPBE − E_PBE\| > 1e-4 Ha on same O₂ MT twin |
| Frame | `frame.meta` energy [eV] matches `last_energy_hartree()` × conversion |
| Perf | T_atv / T_oracle ≤ 10.0, or `perf_exempt` with reason |
| API | No public `*Config` / `make_*` / `build_*`; ctor args + members only |
| Oracle | Offline pin only (`tools/cp2k_oracle/inputs/gpw-uks-revpbe/`); never live CP2K in ctest |

`acceptance_met = true` only when unit + regression green and (within_10x ∨ perf_exempt).
