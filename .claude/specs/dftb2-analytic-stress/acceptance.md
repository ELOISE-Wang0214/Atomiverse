# Acceptance — `dftb2-analytic-stress`

## Correctness

- [x] `DFTB2` with `with_stress=true` on **restricted MioSketch** writes intensive
      Voigt `kStress` [eV/Å³] via **closed-form** pair contraction (no multi-SCC
      strain loop on the product path).
- [x] Analytic \(\sigma\) matches multi-SCC **affine strain-FD** of the same
      total energy on off-eq **H2** (diag components) to ≤ 5×10⁻⁴ eV/Å³
      (unit test `AnalyticStressMatchesStrainFdH2`).
- [x] Regression on off-eq **H2O**: finite energy; finite stress; analytic vs
      strain-FD on σ_xx/σ_yy within 5×10⁻³ eV/Å³; `REFERENCE.md` states
      `oracle: dftb+` and `oracle_version`.
- [x] No public `*Config` / `make_*` / `build_*`. SE-only cmake flags only.
- [x] Open-shell / ThreeOb use multi-SCC strain-FD stress (not claimed analytic);
      restricted MioSketch is the closed-form product path.

## Build / test (scoped)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_dftb2 reg_dftb2-analytic-stress
ctest --test-dir build -R '^(test_cpu_dftb2)$' --output-on-failure
cmake --build build --target run_reg_dftb2-analytic-stress
```

**Result (2026-07-27):** `test_cpu_dftb2` PASS; `run_reg_dftb2-analytic-stress` PASS
(H2O σ_xx |Δ|≈2.0e-4, σ_yy |Δ|≈1.3e-4 eV/Å³).

## Perf

- `perf_exempt=true` — no DFTB+ binary on pin host for fair wall ratio.
