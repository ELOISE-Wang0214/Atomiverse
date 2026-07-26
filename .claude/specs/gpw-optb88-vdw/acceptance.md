# Acceptance — `gpw-optb88-vdw`

## Hard gates

| # | Gate | Pass criterion |
| - | ---- | -------------- |
| 1 | Spec present | `.claude/specs/gpw-optb88-vdw/{spec,acceptance}.md` |
| 2 | API | `XcKind::OptB88VdW` + `XCType::OptB88VdW`; `vdw_df` for DRSLL; no `*Config` / `make_*` |
| 3 | Unit XC | `test_xc_functional` OptB88VdW cases green |
| 4 | Unit grid | `test_cpu_gpw_xc_grid_gga` OptB88 cases green |
| 5 | Regression | `reg_gpw-optb88-vdw` PASS vs REFERENCE.md |
| 6 | Offline pin | `tools/cp2k_oracle/inputs/gpw-optb88-vdw/` + REFERENCE citations |
| 7 | Perf | `perf_exempt=true` (E_disp / XC combo pin) |

## Soft / deferred

- Fair CP2K QS wall ratio when module available on pin host.
- UKS / multi-k / analytic force with OptB88VdW.
