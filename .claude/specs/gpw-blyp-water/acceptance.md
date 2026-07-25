# Acceptance — `gpw-blyp-water`

## Dual gate (must all pass)

| # | Check | Criterion |
| - | ----- | --------- |
| 1 | Unit | `test_cpu_gpw_xc_grid_gga` covers BLYP LDA-limit (constant ρ) and/or Gaussian energy finite + differs from PBE |
| 2 | E2E energy | `reg_gpw-blyp-water` converges SCF; finite \(E_{\mathrm{BLYP}}\) [Ha] |
| 3 | Dual-gate | Same geometry: \(\lvert E_{\mathrm{BLYP}} - E_{\mathrm{PBE}}\rvert > 10^{-4}\) Ha |
| 4 | Absolute | \(\lvert E_{\mathrm{BLYP}} - E_{\mathrm{ref}}\rvert \le\) tol in `REFERENCE.md` (CP2K offline pin preferred; first-green ATV soft band until CP2K available) |
| 5 | Oracle hygiene | No live CP2K spawn in ctest; twin input only under `tools/cp2k_oracle/` |
| 6 | API | Public path uses `XcKind::BLYP` member/ctor — no `*Config` / `make_*` |
| 7 | Perf | `within_10x` **or** `perf_exempt` documented in REFERENCE |

## Done when

- Spec acceptance rows green on isolation branch `wf/gpw-blyp-water`
- Scoped build: `--target atomiverse test_cpu_gpw_xc_grid_gga reg_gpw-blyp-water`
- Commit on branch; `acceptance_met=true` in worker schema
