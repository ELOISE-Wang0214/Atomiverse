# Acceptance — `periodic-urpa-si`

- [x] Spec reviewed (this file + `spec.md`)
- [x] `GpwDft::direct_rpa` runs on Γ UKS; multi-k UKS / OT / unconverged throw
- [x] Unit: `ctest -R 'test_cpu_gpw_rpa'` green
- [x] Regression: `cmake --build build --target reg_periodic-urpa-si run_reg_periodic-urpa-si`
- [x] Dual-gate \(|E_c^{\mathrm{URPA}}-E_c^{\mathrm{RKS\ first-green}}|\le 2\times10^{-3}\) Ha
      (got −0.2658538251 vs −0.2658539776; Δ ≈ 1.5e-7 Ha)
- [x] Offline CP2K input twin present; never live-oracle in ctest
- [x] `perf_exempt` (CP2K frequency RI-RPA ≠ plasmon; wall ~18 s ≤ 60 s)
- [x] Commit on `wf/periodic-urpa-si`; build* wiped
