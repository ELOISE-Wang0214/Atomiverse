# Acceptance — `gpw-mbd-forces`

| # | Criterion | Gate |
| - | --------- | ---- |
| A1 | Isolated free-atom MBD E_disp water dimer matches libMBD pin | \|ΔE\| ≤ 1e-12 Ha (unit) / 1e-8 Ha (reg) |
| A2 | Product MBD forces match independent central FD of same energy | max\|ΔF\| ≤ 1e-6 Ha/Bohr |
| A3 | Net force (translation) | \|Σ F\| ≤ 1e-8 Ha/Bohr |
| A4 | `gpw_forces_*` / evaluator with `mbd=true` includes MBD forces | product path; no energy-only skip |
| A5 | Regression `reg_gpw-mbd-forces` green | exit 0, ≤ 60 s |
| A6 | Architecture | no public Config / factory |
| A7 | Perf | `within_10x` **or** `perf_exempt=true` (geometry-only pin) |

`acceptance_met = true` iff A1–A6 hold and A7 is within_10x or perf_exempt.
