# Acceptance — `second-gen-cpmd-uks`

| ID | Check | How |
| -- | ---- | --- |
| A1 | UKS + `second_gen_cp` builds; mutual exclusives throw | `test_cpu_gpw_second_gen_cpmd` |
| A2 | Warm production SCF iters ≤ `second_gen_cp_max_iter` | regression |
| A3 | RKS SGCP path unchanged (still builds) | unit |
| A4 | O₂ triplet NVE drift envelope + ⟨S²⟩ + finite coords | `reg_second-gen-cpmd-uks` |
| A5 | perf_exempt or within_10x | REFERENCE.md |

Green ⇒ `acceptance_met=true` on worker schema; delete this spec dir on merge to `dev`.
