# Acceptance — `pw-variable-cell-si`

| # | Gate | Criterion |
| - | ---- | --------- |
| 1 | Spec + pin | `REFERENCE.md` cites method, twin, offline pin values, tols |
| 2 | Relax converges | `StressConverged` (or residual within stol) within budget |
| 3 | Lattice | \(\lvert a_{\mathrm{ATV}} - a^\*\rvert \le\) a-tol (Å) |
| 4 | Residual stress | \(\max_d\lvert\sigma_{dd}\rvert \le 1.2\times\) stol |
| 5 | Energy finite | SCF converged; Frame energy [eV] finite |
| 6 | Wall | Regression ≤ 60 s |
| 7 | Perf dual-gate | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) **or** `perf_exempt` with reason (no fair NC-UPF PW twin is acceptable, matching `pw-analytic-stress-si`) |
| 8 | Arch | No public `*Config`; no `make_*`/`build_*`; Frame I/O only |
| 9 | Scoped build/test | `--target atomiverse test_cpu_minimize_* reg_pw-variable-cell-si` (+ related unit tests) green |

`acceptance_met=true` requires gates 1–6, 8–9 and gate 7 (within_10x **or** perf_exempt).
