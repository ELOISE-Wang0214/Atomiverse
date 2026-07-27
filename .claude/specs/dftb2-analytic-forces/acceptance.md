# Acceptance — `dftb2-analytic-forces`

## Correctness

- [x] `DFTB2` with `with_forces=true` does **not** re-run full SCC finite differences
      (no per-coordinate re-SCC; closed-form SK/Pulay/γ/repulsive contraction).
- [x] Analytic \(F=-\nabla E\) matches central re-SCC FD of the same total energy
      on off-eq **H2** to **≤ 5×10⁻³ eV/Å** (unit test).
- [x] On off-eq **H2O**, analytic vs re-SCC FD within **1.0 eV/Å** on gated
      components (charge-response residual on multi-element sketch SK; documented).
- [x] Net force \(\lvert\sum_i F_i\rvert\) ≤ 1×10⁻⁶ eV/Å (translation invariance).
- [x] Energy path finite bound H2O; no re-introduction of public `*Config` / factories.
- [x] `REFERENCE.md` states `oracle: dftb+` (formula / module boundary) and
      `oracle_version`, plus literature citations.

## Build / test (scoped)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_dftb2 reg_dftb2-analytic-forces
ctest --test-dir build -R '^(test_cpu_dftb2)$' --output-on-failure
cmake --build build --target run_reg_dftb2-analytic-forces
```

## Perf

- `perf_exempt=true` — no DFTB+ binary on pin host for fair wall.
