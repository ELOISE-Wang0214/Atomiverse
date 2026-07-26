# Acceptance — gpw-d3-atm

| # | Check | Met |
| - | ----- | --- |
| 1 | Offline CP2K pin recorded (E_disp isolated water dimer + ATM) | |
| 2 | `dftd3_pbe(..., calculate_c9=true)` matches pin ≤ 1e-9 Ha | |
| 3 | Analytic forces vs FD ≤ 1e-7 Ha/Bohr with ATM | |
| 4 | `reg_gpw-d3-atm` green (`GpwDft` d3+d3_atm) | |
| 5 | Default `d3_atm=false` preserves two-body D3 pin | |
| 6 | Scoped build+test only; commit on `wf/gpw-d3-atm` | |

**perf_exempt:** true — deliverable is E_disp pin; D3 ATM is \(O(N^3)\) on 6 atoms (µs); CP2K wall is full QS.
