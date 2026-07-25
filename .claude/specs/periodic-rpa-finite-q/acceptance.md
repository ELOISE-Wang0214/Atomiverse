# Acceptance — periodic-rpa-finite-q

| # | Gate | Pass criterion |
| - | ---- | -------------- |
| 1 | Multi-k SCF | Converged; `C_k.size()==n_kpts`; closed MP mesh |
| 2 | Finite-q product | `direct_rpa` does not throw for \(N_k>1\); `ri==true`; `n_kpts` set |
| 3 | Correlation | `e_corr` finite and **strictly negative**; all \(\Omega>0\) |
| 4 | Identity | `e_total == e_scf + e_corr` (1e-12) |
| 5 | Spectrum size | `omega_ha.size() == N_q * N_k * n_occ * n_virt` with \(N_q=N_k\) (q-mesh = k-mesh) |
| 6 | Unit Γ recovery | Multi-k kernel at \(N_k=1\) ≡ real `direct_rpa_rhf_energy` (`test_cpu_gpw_rpa`) |
| 7 | Absolute band | First-green \(E_c\) pin in REFERENCE.md (setup tol; replace with fair CP2K when available) |
| 8 | Wall | ≤ 60 s e2e |
| 9 | Perf | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) **or** `perf_exempt` + reason |

`acceptance_met=true` only when gates 1–8 green and (9 within_10x or perf_exempt).
