# Acceptance — `dftb2-analytic-forces-threeob`

## Correctness

- [x] `DFTB2` with `sk_set=ThreeOb` and `with_forces=true` does **not** re-run
      full SCC finite differences (no per-coordinate re-SCC).
- [x] Analytic \(F=-\nabla E\) matches central re-SCC FD of the same total
      energy on off-eq organic H2O within **1.0 eV/Å** on gated components
      (unit test + regression; residuals ~0.84 / 0.24 eV/Å on F_O_x / F_H1_y).
- [x] Net force \(\lvert\sum_i F_i\rvert\) ≤ 1×10⁻⁶ eV/Å per Cartesian.
- [x] Energy finite on ThreeOb H2O force path; CH4 energy path unchanged.
- [x] No public `*Config` / factories.
- [x] `REFERENCE.md` states `oracle: dftb+` and `oracle_version` (+ literature).

## Build / test (scoped SE-only)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_dftb2_3ob \
      reg_dftb2-analytic-forces-threeob
ctest --test-dir build -R '^(test_cpu_dftb2_3ob)$' --output-on-failure
cmake --build build --target run_reg_dftb2-analytic-forces-threeob
```

**Result (2026-07-27):** `test_cpu_dftb2_3ob` PASS; `reg_dftb2-analytic-forces-threeob` PASS.

## Perf

- `perf_exempt=true` — no DFTB+ fair wall on pin host; analytic path replaces
  \(6N\) re-SCC force evaluations with one post-SCC contraction.
