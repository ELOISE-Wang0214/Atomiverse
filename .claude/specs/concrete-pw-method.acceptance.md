---
spec: concrete-pw-method
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public PW API is concrete method with members
    type: code
    pass_when: |
      User path is Frame + method members (no public PWConfig bag / PWSolver
      as morphing solver); public factories absent.
    status: pending
  - id: ac-002
    summary: Unit tests for PW method API
    type: runtime
    pass_when: |
      tests/cpu/planewave/test_pw_method passes.
    status: pending
  - id: ac-003
    summary: Existing PW regression retargeted to method API
    type: runtime
    pass_when: |
      Plane-wave regression exits 0 with prior hard-coded energy tolerance on new API.
    status: pending
out_of_scope:
  - Driver-phase rewrite of PW SCF
  - GPW/NDDO concrete migrations
  - CUDA
---

# Acceptance — concrete-pw-method

Done means plane-wave DFT user API matches the concrete-method rule without changing validated energies.
