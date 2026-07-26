# Acceptance — `gpw-uks-rpbe`

| # | Criterion | Gate |
|---|-----------|------|
| 1 | `XcKind::RPBE` UKS build accepted | unit `test_cpu_gpw_uks` |
| 2 | `xc_grid_rpbe_spin` unpolarized collapse | unit |
| 3 | O₂ triplet SCF converges | regression |
| 4 | \(|E - E_{\mathrm{CP2K}}| \le\) REFERENCE tol | regression |
| 5 | \(|E_{\mathrm{RPBE}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha | dual-gate |
| 6 | \(\langle S^2\rangle\) within REFERENCE tol | regression |
| 7 | Frame `kEnergy` [eV] matches internal Ha | regression |
| 8 | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or perf_exempt | REFERENCE |

`acceptance_met=true` only when all rows that apply are green on
`wf/gpw-uks-rpbe` after scoped build+test.
