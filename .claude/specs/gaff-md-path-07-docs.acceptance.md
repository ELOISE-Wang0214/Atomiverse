---
spec: gaff-md-path-07-docs
created: 2026-07-21
criteria:
  - id: ac-001
    summary: force-fields roadmap states Amber special scaling and GAFF regression
    type: docs
    pass_when: |
      docs/roadmap/force-fields.md mentions special_bonds scales including 0.5,
      CoulLong real-space specials, and points at gaff_peo_fragment (or equivalent)
      regression; does not claim multi-type LJ is absent if code has type tables.
    status: pending
  - id: ac-002
    summary: Tutorial subsection documents how to enable Amber weights
    type: docs
    pass_when: |
      ch10 or ch11 contains a short subsection with set_special_bonds(0,0,0.5)
      and load_lammps_data full/charge note plus link to the regression path.
    status: pending
  - id: ac-003
    summary: Remaining gaps are listed honestly
    type: docs
    pass_when: |
      Docs list at least: multi-term fourier component, full .ff text parser,
      CUDA Amber parity, paper-length protocol — as not done or future work.
    status: pending
out_of_scope:
  - New engine features
  - Paper figure regeneration
---

# Acceptance — gaff-md-path-07-docs

Done means roadmap and tutorial text match the GAFF path delivered in 01–06 and
do not over-claim production polymer parity.
