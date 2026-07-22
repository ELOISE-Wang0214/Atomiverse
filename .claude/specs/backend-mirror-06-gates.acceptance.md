---
spec: backend-mirror-06-gates
created: 2026-07-22
criteria:
  - id: ac-001
    summary: architecture scrubs exist and are registered
    type: code
    pass_when: |
      tests/architecture/ contains scrubs covering no facade, no top-level shared
      math, mirror paths, banned short dir names, driver vocab, no runtime/loop/shared
      dirs; each is registered in ctest.
    status: pending
  - id: ac-002
    summary: scrubs pass on the tree after 01–05
    type: runtime
    pass_when: |
      All new architecture scrub tests pass under ctest.
    status: pending
  - id: ac-003
    summary: no old-layout test directories remain
    type: code
    pass_when: |
      tests/cpu and tests/cuda have no directories named md, mc, minimizer, comm,
      async, vib, pw, observer, recorder, kspace, facade.
    status: pending
  - id: ac-004
    summary: architecture_gates regression documents gates
    type: runtime
    pass_when: |
      regressions/layout/architecture_gates/ exists and encodes the expected gate
      names as hard-coded literals (buildable or data-checked per repo regression style).
    status: pending
  - id: ac-005
    summary: full suite green with gates on
    type: runtime
    pass_when: |
      Full check + ctest pass including architecture scrubs.
    status: pending
out_of_scope:
  - Docs prose rewrite; implementing new CUDA physics domains
---

# Acceptance — backend-mirror-06-gates

Done means layout decisions are enforced by ctest scrubs and tests live only on the new tree.
