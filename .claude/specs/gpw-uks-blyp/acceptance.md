# Acceptance — `gpw-uks-blyp`

| # | Criterion | How checked |
|---|-----------|-------------|
| 1 | Spec exists under `.claude/specs/gpw-uks-blyp/` | path |
| 2 | `xc_grid_blyp_spin` unit tests green (collapse + differs from PBE) | `ctest -R 'test_cpu_gpw_(uks|xc_grid_gga)'` |
| 3 | `GpwDft`/`GPWEvaluator` accept `XcKind::BLYP` + UKS occupations | unit + regression |
| 4 | Offline CP2K input under `tools/cp2k_oracle/inputs/gpw-uks-blyp/` | path |
| 5 | `regressions/gpw/gpw-uks-blyp/{REFERENCE.md,reg_*.cpp}` vs pin | `run_reg_gpw-uks-blyp` |
| 6 | Dual-gate \|E_BLYP−E_PBE\| > 1e-4 Ha in regression | reg binary |
| 7 | \(T_\mathrm{atv}/T_\mathrm{oracle}\le 10\) or `perf_exempt` in REFERENCE | REFERENCE |
| 8 | No public `*Config` / `make_*` / factories | review |

**green ⇔ acceptance_met=true** when 1–7 hold (8 is hard fail).
