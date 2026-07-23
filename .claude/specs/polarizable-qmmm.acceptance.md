---
spec: polarizable-qmmm
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public polarizable QM/MM Potential
    type: code
    pass_when: |
      cpu::QmmmPolarizable (or equivalent) Potential runs induced-dipole SCC
      with member parameters; no factory Config bag.
    status: pending
  - id: ac-002
    summary: Unit tests for polarizable embedding
    type: runtime
    pass_when: |
      tests/cpu/multiscale/test_qmmm_polarizable passes.
    status: pending
  - id: ac-003
    summary: Polarizable QM/MM regression
    type: runtime
    pass_when: |
      regressions/multiscale/qmmm_pol_h2o exits 0 with hard-coded energy gate.
    status: pending
out_of_scope:
  - Drude dynamics
  - Link atoms
  - CUDA
---

# Acceptance — polarizable-qmmm

Done means polarizable MM embedding exists on CPU as a public Potential with a textbook regression.
