# Acceptance — `dfpt-incommensurate-si`

| # | Gate | How checked |
| - | ---- | ----------- |
| 1 | Incommensurate \(q\) uses direct Sternheimer (no throw at commensurate gate) | unit + regression |
| 2 | Commensurate / Γ paths unchanged | existing regs not re-run by worker; no code path rewrite of SC branch |
| 3 | `reg_dfpt-incommensurate-si` PASS ≤ 60 s | `run_reg_dfpt-incommensurate-si` |
| 4 | Offline pins in REFERENCE.md | file present |
| 5 | Dual gate: correctness + perf (within_10x or perf_exempt) | REFERENCE + worker schema |
| 6 | No public `*Config` / factories | review |

`acceptance_met=true` when all rows hold on `wf/dfpt-incommensurate-si`.
