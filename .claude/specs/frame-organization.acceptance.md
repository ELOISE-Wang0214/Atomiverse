---
slug: frame-organization
criteria:
  - id: ac-001
    summary: The decision is made on probe evidence, not speculation
    type: code
    pass_when: |
      All four probes (manybody-potentials, md-constraints, path-driver-neb,
      mc-gcmc-npt) are at least code-complete, and their .claude/NOTES.md
      architectural verdicts are summarized into one table of "holes the real
      workloads punched in Frame" — specifically: where EAM's per-atom rho lives,
      where the constraint's x_unc/v_unc shadow arrays live, whether the MC
      Frame-patch forced sub-Frame semantics, and whether ReplicaExpansion's
      expand/collapse forced Frame composability. No decision is recorded before
      this evidence exists.
    status: pending

  - id: ac-002
    summary: A signed-off verdict among (A) backend-only, (B) backend x domain, (C) composable blocks
    type: code
    pass_when: |
      A written verdict picks one option with its evidence, and the user signs off.
      Option (B) backend x domain must be explicitly addressed, since it multiplies
      the conversion boundaries from 3 (materialize/write_back, to_cuda/to_cpu,
      torch::transfer) to 3xD and breaks cross-domain workflows like AIMD (an MD
      driver bound to an SCF evaluator would need two Frame types at once).
    status: pending

  - id: ac-003
    summary: If (C) is chosen — blocks align with molrs, and the 3 boundaries survive
    type: code
    pass_when: |
      Conditional on the verdict being (C): the block model mirrors molrs's Block
      columns (the exchange schema already IS blocks on the other side of the FFI),
      the hand-rolled has_bonded_topology_/ensure_bonded_storage() lazy-allocation in
      cpu/frame/atom.h is replaced by a real optional block, the number of conversion
      boundaries stays at 3, and tests/cpu/test_materialize.cpp's pinned round-trip
      contract still passes unchanged.
    status: pending
---

# Acceptance — frame-organization

A decision spec. No implementation code lands until the verdict is signed off.
Deliberately sequenced AFTER the four probes: they are the evidence.
