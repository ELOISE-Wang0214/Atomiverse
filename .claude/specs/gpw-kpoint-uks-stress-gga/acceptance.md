# Acceptance — gpw-kpoint-uks-stress-gga

## Correctness

1. Multi-k UKS `gpw_stress` with `XcKind::{RevPBE,RPBE,BLYP}` does **not**
   throw for the high-spin H₂ 1×1×2 fixture.
2. Dual-gate: for each of RevPBE / RPBE / BLYP,
   `max|σ_an − σ_fd|` (diagonal) ≤ **1×10⁻⁴** Ha/Bohr³.
3. Off-diagonal product stress ≤ **1×10⁻⁸** (orthorhombic).
4. Activation: |E(RevPBE) − E(PBE)| ≥ **1×10⁻⁴** Ha on the same geometry
   (GGA not silently PADE/PBE).
5. Soft energy pins (first-green ATV) within **1×10⁻³** Ha after fill-in.
6. Unit test `UksKpointGgaAnalyticMatchesFdDiagonal` green.

## Performance

- Target \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) when a same-basis
  CP2K NUMERICAL twin wall is available.
- **`perf_exempt`** allowed if no CP2K module / no same-basis twin on pin host
  (same policy as `gpw-kpoint-uks-stress`).

## Build / test (scoped)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_gpw_stress \
      reg_gpw-kpoint-uks-stress-gga
ctest --test-dir build -R '^(test_cpu_gpw_stress)$' --output-on-failure
cmake --build build --target run_reg_gpw-kpoint-uks-stress-gga
```

## Done when

`acceptance_met=true`, dual-gate green, commit on `wf/gpw-kpoint-uks-stress-gga`.
