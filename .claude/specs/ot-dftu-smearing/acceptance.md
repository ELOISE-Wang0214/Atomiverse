# Acceptance — `ot-dftu-smearing`

| ID | Gate | Tol / criterion |
| -- | ---- | --------------- |
| A0 | `build` accepts multi-k RKS `use_ot` + `plus_u` + FD/MP smear | no throw |
| A1 | OT+U+FD SCF converges | `converged == true` |
| A2 | \(\lvert F_{\mathrm{OT}} - F_{\mathrm{DIIS}}\rvert\) (same U, smear, k) | \(10^{-4}\) Ha |
| A3 | \(\lvert E_U\rvert\) non-trivial | \(> 10^{-8}\) Ha; OT vs DIIS \(E_U\) within \(10^{-3}\) Ha |
| A4 | \(N_{\mathrm{elec}}\) conserved under smear | \(\lvert N_e - N_{\mathrm{target}}\rvert < 10^{-5}\) |
| A5 | Activation \(\lvert F(U)-F(0)\rvert\) (OT, same smear) | \(\ge 10^{-4}\) Ha |
| A6 | Frozen first-green free-energy pin (REFERENCE.md) | \(10^{-3}\) Ha |
| A7 | Multi-k UKS OT + U + smear still throws at build (out of scope) | throws |

## Dual gate

| Item | Criterion |
| ---- | --------- |
| Correctness | A0–A7 green in `reg_ot-dftu-smearing` |
| Perf | \(r = T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\), **or** `perf_exempt` because CP2K has no OT+k twin (DIIS twin only; same pattern as `ot-dftu-kpoints`) |

`acceptance_met = true` when A0–A7 pass and perf gate records either
`within_10x` or `perf_exempt` with reason in REFERENCE.md.
