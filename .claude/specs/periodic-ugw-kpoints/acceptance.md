# Acceptance — periodic-ugw-kpoints

| Gate | Criterion | Met? |
|------|-----------|------|
| API | multi-k UKS `g0w0` no throw; `n_kpts`, `ri`, Newton `converged` | **yes** |
| Unit | Nk=1 UHF α=β ≡ RHF spectral QP 1e-8 | **yes** (`UhfKpointsNk1DualGatesRhf`) |
| Dual-gate | \|gap_UGW − 7.529834\| ≤ 0.5 eV vs multi-k RKS first-green | **yes** (7.698404 eV) |
| Hard | finite positive min-direct QP/KS gaps; HOMO_QP < LUMO_QP α | **yes** |
| Oracle | offline input twin present; no live CP2K in ctest | **yes** |
| Perf | perf_exempt documented | **yes** |
| Wall | e2e ≤ 60 s | **yes** (~16 s) |

`acceptance_met = true` (2026-07-26).
