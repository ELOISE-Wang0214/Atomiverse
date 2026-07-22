---
spec: backend-mirror-03-paths
created: 2026-07-22
criteria:
  - id: ac-001
    summary: old short domain/infra dir names removed
    type: code
    pass_when: |
      Under include/atomiverse/cpu and cuda, directories named md, mc, minimizer, vib,
      pw, observer, recorder, comm, async, kspace do not exist.
    status: pending
  - id: ac-002
    summary: full-word directories present on cpu
    type: code
    pass_when: |
      include/atomiverse/cpu contains dynamics, montecarlo, minimize, communication,
      asynchronous, planewave, observe, record, vibration (as applicable to prior
      content), longrange, potentials.
    status: pending
  - id: ac-003
    summary: cuda mirrors implemented domain paths
    type: code
    pass_when: |
      For every domain/infra directory that exists on both backends after the rename,
      the relative path spelling matches (e.g. both dynamics/, both longrange/).
    status: pending
  - id: ac-004
    summary: tests mirror new paths; old test dirs gone
    type: code
    pass_when: |
      tests/cpu/md, tests/cpu/mc, tests/cpu/minimizer, tests/cpu/comm (etc. old names)
      do not exist; corresponding tests live under new full-word paths and are in ctest.
    status: pending
  - id: ac-005
    summary: path_words regression uses new includes only
    type: runtime
    pass_when: |
      regressions/layout/path_words/ builds with includes under
      atomiverse/cpu/dynamics and atomiverse/cpu/communication only (new paths).
    status: pending
  - id: ac-006
    summary: full suite green after renames
    type: runtime
    pass_when: |
      Full check + ctest pass with zero references to deleted directory path segments
      in compile commands for product targets.
    status: pending
out_of_scope:
  - Driver class rename; Backend dissolve; docs-only harness rewrite
---

# Acceptance — backend-mirror-03-paths

Done means full-word directory names on disk for both backends, tests rewritten, no old path shims.
