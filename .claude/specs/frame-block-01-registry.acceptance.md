---
slug: frame-block-01-registry
criteria:
  - id: ac-001
    summary: The registry is THE single field enumeration
    type: code
    pass_when: |
      grow, set_capacity, copy_row, pack_exchange, and unpack_exchange contain
      NO per-field hand-written lists — each iterates the registry. A grep gate
      (same style as test_stage_typing_scrub) fails the build if a per-atom
      field name is enumerated outside core/frame_schema.h in those functions.
      Every cpu::Atom member has exactly one schema row (a unit test
      cross-checks count and names both ways).
    status: pending

  - id: ac-002
    summary: molrs-aligned mapping API with compile-time typed keys
    type: code
    pass_when: |
      frame.has("atoms","x")-style presence checks and frame.span(kX) typed
      access work; keys reuse the molrs exchange schema names 1:1; an extension
      column can be inserted at runtime, is registry-owned, survives
      grow/set_capacity, and rides copy_row/pack_exchange when its semantic
      tags say so. Bookkeeping (nlocal/nghost/nmax, tag map, SystemScalars,
      Box) stays typed attrs — NOT reachable via string keys.
    status: pending

  - id: ac-003
    summary: The 13.7x cliff is machine-guarded
    type: code
    pass_when: |
      Debug builds count per-launch key resolutions and step_once asserts the
      count stays below the configured bound; a deliberately-violating test
      component trips the assertion. Release builds compile the counter out
      (zero overhead, verified by inspecting that the assert path is under
      NDEBUG or an ATV debug macro).
    status: pending

  - id: ac-004
    summary: Behavior is bitwise-preserved
    type: code
    pass_when: |
      tests/cpu/test_materialize.cpp passes UNCHANGED (golden round-trip
      contract). The existing equivalence pins hold: AtomGrow /
      AtomCopyRow / transaction-rollback / migration tests all green.
      Full ctest suite passes.
    status: pending

  - id: ac-005
    summary: Performance gate — the LJ ⚡ kernel does not regress
    type: performance
    pass_when: |
      benchmarks/md/lj/performance before/after comparison shows <= 2% per-step
      regression on the benchmarked LJ case (the tree's only ⚡ kernel, most
      sensitive sentinel). Run at acceptance time, not in ctest.
    status: pending

  - id: ac-006
    summary: Memory gate — memory_bytes parity
    type: performance
    pass_when: |
      Atom::memory_bytes() ("Atoms" category) for the same LJ case is equal
      before/after up to the registry's constant overhead (< 16 KB).
    status: pending
---

# Acceptance — frame-block-01-registry

Step 1 of the signed Frame verdict (variant C route): registry + ColKey +
mapping API, members untouched. Gates from the verdict travel with every
frame-block spec.
