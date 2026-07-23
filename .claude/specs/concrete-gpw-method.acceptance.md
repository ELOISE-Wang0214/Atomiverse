---
spec: concrete-gpw-method
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public GPW API is concrete method with members
    type: code
    pass_when: |
      User path is Frame + method members; public GpwConfig/GPWSolver bag removed
      or tombstoned; no factories.
    status: pending
  - id: ac-002
    summary: Unit tests for GPW method API
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_gpw_method passes.
    status: pending
  - id: ac-003
    summary: GPW regressions pass on method API
    type: runtime
    pass_when: |
      Existing GPW regressions exit 0 with prior hard-coded tolerances on new API.
    status: pending
out_of_scope:
  - New GPW physics
  - PW/NDDO migrations
  - CUDA
---

# Acceptance — concrete-gpw-method

Done means GPW user API matches concrete-method architecture without changing validated numbers.
