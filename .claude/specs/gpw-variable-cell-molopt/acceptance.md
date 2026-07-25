# Acceptance — `gpw-variable-cell-molopt`

## Gates

| # | Check | Pass criterion |
| - | ----- | -------------- |
| 1 | Spec present | `.claude/specs/gpw-variable-cell-molopt/{spec,acceptance}.md` |
| 2 | Offline pin | `REFERENCE.md` cites CP2K 2026.1 a\* + E\* from CELL_OPT |
| 3 | Regression | `regressions/gpw/gpw-variable-cell-molopt/` + CMake `atv_add_regression` |
| 4 | Dual gate | `reg_gpw-variable-cell-molopt` PASS: \|Δa\|≤0.02 Å, \|ΔE\|≤1e-3 Ha, residual stress, b/c fixed |
| 5 | Budget | regression wall ≤ 60 s |
| 6 | Perf | within_10x vs T_cp2k CELL_OPT **or** perf_exempt + reason |
| 7 | Arch | no public Config bags / factories; Frame I/O |
| 8 | Commit | on `wf/gpw-variable-cell-molopt` only |

## Done when

`acceptance_met = true`, dual gate green, `status = green`, head committed.
