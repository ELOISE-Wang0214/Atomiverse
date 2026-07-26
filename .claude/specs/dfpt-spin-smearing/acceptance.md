# Acceptance — dfpt-spin-smearing

| # | Gate | Met? |
| - | ---- | ---- |
| 1 | Spec exists under `.claude/specs/dfpt-spin-smearing/` | yes |
| 2 | Γ DFPT accepts spin_polarized + smearing (no spin+smear throw) | yes |
| 3 | Offline pin + REFERENCE.md + oracle note | yes |
| 4 | `reg_dfpt-spin-smearing` PASS (physics + FD pin) | yes |
| 5 | Unit tests updated (`test_cpu_pw_dfpt`) | yes |
| 6 | Dual gate / scoped ctest green | yes (rebased on `dev` @ `75655ede`; `test_cpu_pw_dfpt` PASS; `run_reg_dfpt-spin-smearing` PASS ~35 s wall) |
| 7 | Commit on `wf/dfpt-spin-smearing` | yes |

`acceptance_met=true` when 2–6 green and commit landed.
