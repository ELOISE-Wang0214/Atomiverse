---
spec: backend-mirror-02-absorb
created: 2026-07-22
criteria:
  - id: ac-001
    summary: top-level longrange and potentials dirs gone
    type: code
    pass_when: |
      include/atomiverse/longrange/ and include/atomiverse/potentials/ do not exist;
      src/longrange/ does not exist.
    status: pending
  - id: ac-002
    summary: cuda kspace renamed to longrange
    type: code
    pass_when: |
      include/atomiverse/cuda/longrange/ exists; include/atomiverse/cuda/kspace/ does not;
      product code has no include path containing cuda/kspace.
    status: pending
  - id: ac-003
    summary: pair_kernel lives under cuda/potentials
    type: code
    pass_when: |
      pair kernel symbols used by CUDA pair styles are declared under
      include/atomiverse/cuda/potentials/ (not top-level potentials/).
    status: pending
  - id: ac-004
    summary: longrange unit tests under tests/cpu/longrange
    type: code
    pass_when: |
      tests/longrange/ does not exist; former FMM/operator unit tests live under
      tests/cpu/longrange/ and are registered with ctest.
    status: pending
  - id: ac-005
    summary: backend_owned_longrange regression builds
    type: runtime
    pass_when: |
      regressions/layout/backend_owned_longrange/ builds using only
      atomiverse/cpu/longrange/ public headers (no top-level longrange include).
    status: pending
  - id: ac-006
    summary: full suite green after absorb
    type: runtime
    pass_when: |
      Full project check + ctest pass after absorb/rename.
    status: pending
out_of_scope:
  - Domain directory English renames; Driver rename; Backend dissolve
---

# Acceptance — backend-mirror-02-absorb

Done means every operator lives under a backend tree, top-level shared math is gone, cuda uses `longrange/` not `kspace/`, tests rewritten to backend paths.
