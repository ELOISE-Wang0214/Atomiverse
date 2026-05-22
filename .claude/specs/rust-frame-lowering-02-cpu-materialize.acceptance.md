---
slug: rust-frame-lowering-02-cpu-materialize
criteria:
  - id: ac-001
    summary: cpu::Frame gains a charge field and a shared SystemScalars
    type: code
    pass_when: |
      cpu::Frame::Atoms has a `charge` xtensor field with Doxygen units, and
      cpu::Frame embeds a SystemScalars sub-struct defined in
      include/atomiverse/core/system_scalars.h.
    status: pending
  - id: ac-002
    summary: molrs.Frame <-> cpu::Frame materialization round-trips
    type: runtime
    pass_when: |
      tests/cpu/test_materialize.cpp materializes a molrs.Frame into a
      cpu::Frame and writes it back, asserting positions / velocities / forces /
      charge / topology / box / system scalars are element-wise equal within
      dtype-cast tolerance; the test passes under ctest.
    status: pending
  - id: ac-003
    summary: missing-field behaviour is loud for required, silent-zero for optional
    type: runtime
    pass_when: |
      tests/cpu/test_materialize.cpp asserts cpu::materialize throws naming the
      field when atoms.x or atoms.type is absent, and leaves zeros (no throw)
      when velocity or charge is absent; the assertions pass under ctest.
    status: pending
  - id: ac-004
    summary: a cpu::md::MD run on a materialized cpu::Frame yields finite energy
    type: runtime
    pass_when: |
      A cpu::md::MD run started from a cpu::Frame produced by cpu::materialize
      runs several steps and leaves a finite frame energy; asserted in
      tests/cpu/test_materialize.cpp and passing under ctest.
    status: pending
---

# Acceptance criteria

- **ac-001** is the structural deliverable — the `charge` field closes the
  Capability/field gap, and `SystemScalars` is the cross-backend shared shape
  `-04` will reuse. Verified by reading the headers, hence `code`.
- **ac-002** is the core function of this spec: data survives the
  materialize/write-back round trip.
- **ac-003** pins the boundary contract — required fields fail loud, optional
  fields default to zero — so downstream drivers can rely on it.
- **ac-004** proves a materialized frame is actually runnable, the practical
  end of the CPU ingress path.
