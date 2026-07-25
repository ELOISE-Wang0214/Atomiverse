# Acceptance — `gpw-uks-pbe`

| # | Criterion | How checked |
|---|-----------|-------------|
| 1 | Spec exists under `.claude/specs/gpw-uks-pbe/` | path |
| 2 | `xc_grid_gga_spin` / PBE spin path unit tests green | `ctest -R test_cpu_gpw_uks` + `test_cpu_gpw_xc_grid_gga` |
| 3 | `GpwDft`/`GPWEvaluator` accept `XcKind::PBE` + UKS occupations | unit + regression |
| 4 | Offline CP2K input under `tools/cp2k_oracle/inputs/gpw-uks-pbe/` | path |
| 5 | `regressions/gpw/gpw-uks-pbe/{REFERENCE.md,reg_*.cpp}` vs pin | `run_reg_gpw-uks-pbe` |
| 6 | \(T_\mathrm{atv}/T_\mathrm{oracle}\le 10\) or `perf_exempt` in REFERENCE | REFERENCE |
| 7 | No public `*Config` / `make_*` / factories | review |

**green ⇔ acceptance_met=true** when 1–6 hold.
