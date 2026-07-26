# Acceptance — periodic-ugw-si

- [x] Spec written under `.claude/specs/periodic-ugw-si/`
- [x] `direct_rpa_uhf` in `mp2_kernels` — E_c dual-gates RHF on α=β tensors
- [x] `GpwDft::g0w0` Γ UKS path (spin-conserving dRPA + spectral Σ + MT Σ_x + spin v_xc)
- [x] Multi-k UKS / OT / unconverged guards; RKS path unchanged
- [x] Unit `test_cpu_gpw_g0w0` green
- [x] Regression `reg_periodic-ugw-si` green — dual-gate QP gap 5.719564 eV vs RKS 5.719545 eV
- [x] Offline twin `tools/cp2k_oracle/inputs/periodic-ugw-si.inp`
- [x] `perf_exempt` (no fair CP2K spectral twin); wall ~14 s @ 8 OMP
- [x] No public Config / factories

**acceptance_met = true** (2026-07-26)
