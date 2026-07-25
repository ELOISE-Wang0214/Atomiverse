---
spec: gpw-uks-stress
created: 2026-07-25
criteria:
  - id: ac-001
    summary: Γ-UKS analytic stress no longer throws
    type: code
    pass_when: |
      gpw_stress / energy_fixed_P accept is_uks() at Γ; multi-k UKS still
      throws with a clear message (out of scope).
    status: met
  - id: ac-002
    summary: Unit dual-gate UKS analytic vs FD
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_stress.cpp includes a cheap Γ-UKS fixture
      (H₂ high-spin or O₂-scale) with max|σ_an−σ_fd|_diag below unit tol.
    status: met
  - id: ac-003
    summary: Regression e2e dual-gate + energy pin
    type: runtime
    pass_when: |
      regressions/gpw/gpw-uks-stress/reg_gpw-uks-stress builds under
      ATV_BUILD_GPW, never launches CP2K, dual-gate max diag |σ_an−σ_fd|
      ≤ REFERENCE tol, |E−E_cp2k| ≤ energy tol, off-diag ~0.
    status: met
  - id: ac-004
    summary: Offline pin documented
    type: docs
    pass_when: |
      REFERENCE.md cites CP2K version (or documents pin host without binary),
      twin input path under tools/cp2k_oracle/inputs/gpw-uks-stress/, and
      hard-coded pins used by the regression.
    status: met
  - id: ac-005
    summary: API iron laws
    type: code
    pass_when: |
      No public *Config / make_* / build_* for this feature.
    status: met
  - id: ac-006
    summary: Perf gate
    type: runtime
    pass_when: |
      T_atv/T_oracle ≤ 10 when fair CP2K wall exists; else perf_exempt with
      reason (no CP2K binary on pin host is acceptable).
    status: met
out_of_scope:
  - k-point UKS stress
  - triclinic shear stress
  - DFT+U / D4 / hybrid / meta-GGA stress
  - CUDA
  - live oracle in ctest
---

# Acceptance — `gpw-uks-stress`

**Done** means product Γ-UKS analytic fixed-P cell stress is available via
`gpw_stress` / `GPWEvaluator::compute_stress`, dual-gated against multi-SCF
strain-FD on the O₂ UKS twin, with offline energy pin documented.

## Green condition

`acceptance_met=true` when AC-001–AC-005 pass and (AC-006 within_10x **or**
`perf_exempt`).
