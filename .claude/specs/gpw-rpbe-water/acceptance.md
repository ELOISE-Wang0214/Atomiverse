# Acceptance — `gpw-rpbe-water`

| # | Criterion | Gate |
| - | --------- | ---- |
| 1 | `XCType::RPBE` pure GGA metadata | unit `test_xc_functional` |
| 2 | σ=0 matches PBE; finite-σ differs | unit `test_xc_functional` + `test_cpu_gpw_xc_grid_gga` |
| 3 | `XcKind::RPBE` SCF on H₂O MT converges | `reg_gpw-rpbe-water` |
| 4 | Dual-gate \|E_RPBE−E_PBE\| > 1e-4 Ha | `reg_gpw-rpbe-water` |
| 5 | Soft absolute pin / offline CP2K twin documented | `REFERENCE.md` |
| 6 | No live oracle; public eV on Frame | code review / reg |
| 7 | Perf: ≤10× CP2K or `perf_exempt` | `REFERENCE.md` |

**green** when 1–6 pass and (7 within_10x **or** perf_exempt).
