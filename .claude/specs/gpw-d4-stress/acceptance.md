# Acceptance — `gpw-d4-stress`

| # | Criterion | Gate |
| - | --------- | ---- |
| 1 | Analytic stress with `d4=true` does not throw (Γ) | unit + reg |
| 2 | Dual-gate max\|\σ_an−σ_fd\| ≤ tol on product twin | reg |
| 3 | E_disp matches dftd4 3.7.0 pin (water dimer periodic) | reg ≤ 1e-8 Ha |
| 4 | Off-diag product σ ≈ 0 | reg ≤ 1e-8 |
| 5 | Unit cheap H₂+D4 analytic vs FD finite + residual | unit |
| 6 | `REFERENCE.md` cites offline pin (dftd4 / dual-gate) | file |
| 7 | Scoped build+test green; perf_exempt OK | run |

`acceptance_met=true` when all rows pass on `wf/gpw-d4-stress`.
