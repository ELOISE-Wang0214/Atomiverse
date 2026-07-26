# Acceptance — `second-gen-cpmd-uks-ot`

| ID | Gate | How |
| -- | ---- | --- |
| A1 | `GPWEvaluator::build` accepts UKS + `second_gen_cp` + `use_ot` | unit |
| A2 | `GpwDft::run` accepts UKS OT + `guess_P`/`guess_P_beta` | unit |
| A3 | Production MD force steps use partial OT (iters ≤ budget) | e2e |
| A4 | Relative NVE drift ≤ 0.10 on short O₂ traj | e2e |
| A5 | ⟨S²⟩ stays in [1.5, 2.5]; O–O in soft window | e2e |
| A6 | Offline oracle input + REFERENCE.md | docs |
| A7 | Dual gate: within_10x **or** `perf_exempt` | REFERENCE |

`acceptance_met=true` when A1–A7 hold and scoped `ctest` + `run_reg_*` green.
