# Acceptance — periodic-gw-finite-q

| # | Gate | Pass criterion |
| - | ---- | -------------- |
| 1 | Multi-k SCF | Converged; `C_k.size()==n_kpts`; closed MP mesh |
| 2 | Finite-q product | `g0w0` does not throw for \(N_k>1\); `ri==true`; `n_kpts` set |
| 3 | Newton | `converged==true` (all \(k\), all MO QP solves) |
| 4 | Gaps | Min-direct QP and KS gaps finite and **strictly positive** |
| 5 | Ordering | On report spectrum: \(\varepsilon_{\mathrm{HOMO}}^{\mathrm{QP}} < \varepsilon_{\mathrm{LUMO}}^{\mathrm{QP}}\) |
| 6 | Occupations | `n_occ==8`, `n_virt==2` on the pin deck |
| 7 | Unit Γ recovery | Multi-k kernel at \(N_k=1\) ≡ real spectral QP (`test_cpu_gpw_g0w0`) |
| 8 | Absolute band | First-green `gap_ev` pin in REFERENCE.md (setup tol; replace with fair CP2K when available) |
| 9 | Wall | ≤ 60 s e2e |
| 10 | Perf | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) **or** `perf_exempt` + reason |

`acceptance_met=true` only when gates 1–9 green and (10 within_10x or perf_exempt).
