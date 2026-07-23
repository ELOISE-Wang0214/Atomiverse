---
spec: concrete-nddo-method
created: 2026-07-23
criteria:
  - id: ac-001
    summary: NDDOSolver not a public bag API
    type: code
    pass_when: |
      Public user path is AM1/RM1 (evaluator or method members) + Frame;
      NDDOSolver is internal or tombstoned; no NDDOSystem public return.
    status: pending
  - id: ac-002
    summary: API surface unit test
    type: runtime
    pass_when: |
      tests/cpu/semiempirical/test_nddo_api_surface passes.
    status: pending
  - id: ac-003
    summary: Existing AM1 regressions still green
    type: runtime
    pass_when: |
      regressions/semiempirical/am1_water_force (and related) exit 0 unchanged gates.
    status: pending
out_of_scope:
  - New semiempirical parameterizations
  - DFTB product tables
  - CUDA
---

# Acceptance — concrete-nddo-method

Done means NDDO public API matches concrete-method architecture while AM1 regressions stay green.
