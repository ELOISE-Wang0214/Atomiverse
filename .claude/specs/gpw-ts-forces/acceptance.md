# Acceptance — `gpw-ts-forces`

| # | Criterion | Gate |
| - | --------- | ---- |
| A1 | Isolated free-atom TS E_disp water dimer matches libMBD pin | \|ΔE\| ≤ 1e-12 Ha (unit) |
| A2 | Periodic pure-RS E_disp (rcut 50 Å, 12 Å cell) matches pin | \|ΔE\| ≤ 1e-12 Ha (unit) / 1e-8 Ha (reg) |
| A3 | Product TS forces match independent central FD of same energy | max\|ΔF\| ≤ 1e-6 Ha/Bohr |
| A4 | Net force (translation) | \|Σ F\| ≤ 1e-8 Ha/Bohr |
| A5 | `gpw_forces_*` / evaluator with `ts=true` includes TS forces | \(F(\mathrm{ts})-F(\mathrm{plain})\) matches standalone; no energy-only skip |
| A6 | Regression `reg_gpw-ts-forces` green | exit 0, ≤ 60 s |
| A7 | Architecture | no public Config / factory |
| A8 | Perf | `within_10x` **or** `perf_exempt=true` (geometry-only pin; CP2K has no stock TS) |

`acceptance_met = true` iff A1–A7 hold and A8 is within_10x or perf_exempt.
