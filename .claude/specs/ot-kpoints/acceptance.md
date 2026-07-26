---
spec: ot-kpoints
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Multi-k RKS + use_ot builds and converges
    type: runtime
    pass_when: |
      GpwDft with use_ot=true and k_mesh or explicit k_points (closed-shell
      RKS, no smear, no DFT+U) no longer throws at build; SCF converges on
      the regression system.
    status: met
  - id: ac-002
    summary: OT vs DIIS multi-k energy parity
    type: runtime
    pass_when: |
      On the same geometry/grids/k-mesh, |E_OT − E_DIIS| ≤ tol in unit test
      and regression (≤ 1e-4 Ha class on compact systems).
    status: met
  - id: ac-003
    summary: Explicit-Γ OT matches real Γ OT
    type: runtime
    pass_when: |
      k_points = {(0,0,0)} + use_ot recovers the real Γ OT energy within
      1e-6 Ha on a compact H₂ toy (identity gate).
    status: met
  - id: ac-004
    summary: Fail-fast on unsupported combos
    type: runtime
    pass_when: |
      Multi-k UKS + use_ot, multi-k OT + smearing, multi-k OT + plus_u still
      throw at build/run with clear messages (no silent fallback).
    status: met
  - id: ac-005
    summary: Regression + offline pin contract
    type: runtime
    pass_when: |
      regressions/gpw/ot-kpoints builds under ATV_BUILD_GPW, never launches
      CP2K, asserts gates in REFERENCE.md; perf_exempt if no fair CP2K wall
      on pin host, else within_10x.
    status: met
  - id: ac-006
    summary: API iron laws
    type: code
    pass_when: |
      No new public *Config / make_* / build_* for this feature; parameters
      remain GpwDft members.
    status: met
out_of_scope:
  - multi-k UKS OT
  - OT + DFT+U
  - OT + metallic smearing
  - CUDA
  - live CP2K in ctest
---

# Acceptance — `ot-kpoints`

**Done** means multi-k closed-shell GPW can use `use_ot` productively, with
OT/DIIS fixed-point parity and documented offline pin strategy, without
reopening banned public Config/factory surface.
