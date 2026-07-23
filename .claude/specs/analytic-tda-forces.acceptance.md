---
spec: analytic-tda-forces
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Analytic CIS/TDA forces on RHF public API
    type: code
    pass_when: |
      RHF exposes analytic TDA nuclear forces [eV/Å] without public IntegralPack;
      FD path remains available as oracle or explicit option.
    status: pending
  - id: ac-002
    summary: Unit tests match analytic vs FD
    type: runtime
    pass_when: |
      tests/cpu/scf/test_analytic_tda_forces passes comparing analytic to FD.
    status: pending
  - id: ac-003
    summary: H2 analytic TDA force regression
    type: runtime
    pass_when: |
      regressions/scf/h2_analytic_tda_force exits 0 with hard-coded tolerance gate.
    status: pending
out_of_scope:
  - TDDFT forces
  - NAC vectors
  - CUDA
---

# Acceptance — analytic-tda-forces

Done means closed-shell RHF TDA nuclear forces are analytic, agree with FD on a textbook system, and ship a regressions/ gate.
