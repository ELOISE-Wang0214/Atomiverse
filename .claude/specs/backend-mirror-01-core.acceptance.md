---
spec: backend-mirror-01-core
created: 2026-07-22
criteria:
  - id: ac-001
    summary: facade/ directory removed from include tree
    type: code
    pass_when: |
      Path include/atomiverse/facade/ does not exist.
    status: pending
  - id: ac-002
    summary: core holds Plan/Stage/Role/Verdict contracts
    type: code
    pass_when: |
      include/atomiverse/core/plan.h, stage.h, role.h (or equivalent agreed names in the
      spec Design), and verdict/iteration status headers exist and define Plan, slot_of,
      RoleFlag/ComponentRole, and Verdict.
    status: pending
  - id: ac-003
    summary: zero atomiverse/facade includes in product code
    type: code
    pass_when: |
      rg 'atomiverse/facade' over include/, src/, tests/, examples/, regressions/,
      python/, app/ returns no matches.
    status: pending
  - id: ac-004
    summary: tests/facade removed; tests/core owns scrub
    type: code
    pass_when: |
      tests/facade/ does not exist; tests/core/stage_typing_scrub.sh exists and is
      registered in ctest.
    status: pending
  - id: ac-005
    summary: core_contract_includes regression builds
    type: runtime
    pass_when: |
      regressions/layout/core_contract_includes/ builds as an out-of-library consumer
      including atomiverse/core plan+stage headers only (no facade path).
    status: pending
  - id: ac-006
    summary: full suite green after move
    type: runtime
    pass_when: |
      cmake build + ctest --test-dir build --output-on-failure passes for the
      configured CPU test set used by the project check line.
    status: pending
out_of_scope:
  - Domain renames, Backend dissolve, shared math absorb, Driver class rename
---

# Acceptance — backend-mirror-01-core

Done means: contracts live only under `core/`, `facade/` is gone, tests and includes were rewritten (not shimmed), and the suite is green.
