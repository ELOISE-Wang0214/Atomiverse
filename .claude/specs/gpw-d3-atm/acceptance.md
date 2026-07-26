# Acceptance — gpw-d3-atm

| # | Check | Met |
| - | ----- | --- |
| 1 | Offline CP2K pin recorded (E_disp isolated water dimer + ATM) | yes |
| 2 | `dftd3_pbe(..., calculate_c9=true)` matches pin ≤ 1e-9 Ha | yes |
| 3 | Analytic forces vs FD ≤ 1e-7 Ha/Bohr with ATM | yes |
| 4 | `reg_gpw-d3-atm` green (`GpwDft` d3+d3_atm) | yes |
| 5 | Default `d3_atm=false` preserves two-body D3 pin | yes |
| 6 | Scoped build+test only; commit on `wf/gpw-d3-atm` | yes |

**Dual gate (2026-07-26):** `test_cpu_gpw_dftd3` PASS (1.94 s);
`run_reg_gpw-d3-atm` PASS — E_disp = −8.14141236471528e−4 Ha vs pin
−8.1414118193e−4 Ha (tol 1e−9). WALL_ATV ≈ 3.2 s.

**perf_exempt:** true — deliverable is E_disp pin; D3 ATM is \(O(N^3)\) on 6 atoms (µs); CP2K wall is full QS (~10.8 s).
