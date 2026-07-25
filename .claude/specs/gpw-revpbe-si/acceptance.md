# Acceptance — `gpw-revpbe-si`

## Correctness

1. `XCType::RevPBE` / `XcKind::RevPBE` compile and are pure GGA
   (`is_gga`, `hybrid_exchange_fraction == 0`).
2. At \(\sigma = 0\), `evaluate_xc(RevPBE)` matches `evaluate_xc(PBE)`
   (same PW92 LDA base; \(\kappa\) only enters \(F_x(s)\)).
3. At finite \(\sigma\), \(|\varepsilon_{\mathrm{revPBE}}-\varepsilon_{\mathrm{PBE}}| > 10^{-8}\).
4. Regression dual-gate: SCF converges; finite \(E_{\mathrm{revPBE}}\);
   \(|E_{\mathrm{revPBE}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha on Si₂ MT twin.
5. Soft absolute band vs first-green ATV (or offline CP2K \(E_{\mathrm{tot}}\)
   when binary available); tol \(10^{-3}\) Ha until twin replaces pin.
6. No live oracle spawn in ctest / regression binary.

## Performance

- Target \(r = T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) when CP2K wall
  twin is available.
- If no CP2K on pin host: `perf_exempt = true` (document reason), same as
  `gpw-pbesol-si` / `gpw-blyp-water`.

## Dual gate (scoped)

```bash
cmake --build build -j$(nproc) --target atomiverse \
  test_xc_functional test_cpu_gpw_xc_grid_gga reg_gpw-revpbe-si
ctest --test-dir build -R '^(test_xc_functional|test_cpu_gpw_xc_grid_gga)$' --output-on-failure
cmake --build build -j$(nproc) --target run_reg_gpw-revpbe-si
```

## Exit

`acceptance_met=true` when unit + regression green and pin/REFERENCE written.
Commit on `wf/gpw-revpbe-si`; wipe `build*`.
