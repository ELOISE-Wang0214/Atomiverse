# Acceptance — `gpw-vdw-scf-potential`

## Dual gate

1. **Correctness (hard)**
   - Unit: synthetic-grid FD of \(E_c^{\mathrm{nl}}\) vs \(\int V_{\mathrm{nl}}\,\delta\rho\)
     for DRSLL and rVV10 (`test_cpu_gpw_vdw_df`, `test_cpu_gpw_rvv10`).
   - Regression `reg_gpw-vdw-scf-potential`:
     - Γ PBE + rVV10 water dimer (same geometry as ice pins): SC \(E_{\mathrm{disp}}\)
       within **1×10⁻⁴ Ha** of offline CP2K 2026.1 pin.
     - Analytic forces vs central FD (\(h=5\times10^{-4}\) Bohr) max |ΔF|
       ≤ **5×10⁻⁴ Ha/Bohr** on O atoms (Γ single-grid, tol loose for float32
       kernel + FD noise).
2. **Perf**
   - `perf_exempt = true`: SC RPS is \(O(N_{\mathrm{SCF}}\times n_\alpha)\) FFTs;
     primary deliverable is potential/force correctness (same policy as ice
     E_disp pins). Optional fair wall noted in REFERENCE if available.

## Out of scope this slug

- Analytic stress with nonlocal virial.
- Multigrid product.
- optB88-vdW / stock VV10.
