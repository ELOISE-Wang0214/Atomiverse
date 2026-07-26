# Acceptance — `gpw-mbd-kpoint-forces`

| # | Criterion | Gate |
| - | --------- | ---- |
| A1 | Multi-k force assembly folds MBD when `mbd=true` | product path; no skip |
| A2 | Isolated free-atom MBD E_disp water dimer (libMBD pin) | \|ΔE\| ≤ 1e-8 Ha (reg) |
| A3 | Product MBD forces vs independent central FD | max\|ΔF\| ≤ 1e-6 Ha/Bohr |
| A4 | Translation invariance of MBD forces | \|Σ F_MBD\| ≤ 1e-8 Ha/Bohr |
| A5 | Multi-k wiring: \(F^{k+\mathrm{MBD}}-F^{k}\approx F_{\mathrm{MBD}}\) | max\|Δ\| ≤ 1e-8 Ha/Bohr |
| A6 | Regression `reg_gpw-mbd-kpoint-forces` green | exit 0, ≤ 60 s |
| A7 | Architecture | no public Config / factory |
| A8 | Perf | `within_10x` **or** `perf_exempt=true` |

`acceptance_met = true` iff A1–A7 hold and A8 is within_10x or perf_exempt.
