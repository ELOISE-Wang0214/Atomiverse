# Acceptance — gpw-dftu-kpoints

| # | Criterion | Gate |
|---|-----------|------|
| 1 | `build()` allows k-mesh/k_points + `plus_u` | unit |
| 2 | Γ-limit: complex k={Γ} + U ≡ real Γ-UKS + U within 1e-6 Ha | unit/reg |
| 3 | MP multi-k + U SCF converges; \|E(U)−E(0)\| ≥ 1e-4 Ha | reg |
| 4 | Offline CP2K pin in `regressions/gpw/gpw-dftu-kpoints/REFERENCE.md` | docs |
| 5 | Dual gate: green ctest + within_10x or perf_exempt | reg |
| 6 | Multi-k DFT+U stress still throws (out of scope) | unit |

`acceptance_met=true` only when 1–6 hold on Release scoped build.
