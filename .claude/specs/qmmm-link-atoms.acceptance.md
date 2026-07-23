---
spec: qmmm-link-atoms
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Link-atom boundary on public QM/MM Potential
    type: code
    pass_when: |
      Public QM/MM path accepts covalent cut pairs and places H caps with
      force projection; no Config bag factory.
    status: pending
  - id: ac-002
    summary: Unit tests for link placement
    type: runtime
    pass_when: |
      tests/cpu/multiscale/test_link_atom passes.
    status: pending
  - id: ac-003
    summary: Link-atom QM/MM regression
    type: runtime
    pass_when: |
      regressions/multiscale/qmmm_link_ch3oh exits 0 with hard-coded gate.
    status: pending
out_of_scope:
  - Polarizable MM
  - ONIOM
  - CUDA
---

# Acceptance — qmmm-link-atoms

Done means covalent QM/MM cuts are supported via H link atoms with a regressions/ textbook gate.
