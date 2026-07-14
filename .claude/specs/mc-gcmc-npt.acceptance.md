---
slug: mc-gcmc-npt
criteria:
  - id: ac-001
    summary: Variable-N is LIFTED from comm/, not reinvented in mc/
    type: code
    pass_when: |
      The row-copy primitive is promoted to a public cpu::Atom::copy_row(dst, src),
      and src/cpu/comm/atom_migration.cpp's private static move_row (atom_migration.cpp:134,
      previously unreachable behind #ifdef ATV_HAS_MPI) is DELETED and rewired to call it —
      single definition, no duplication. GCMC's delete/rollback payload uses the EXISTING
      public, MPI-free comm::pack_exchange / unpack_exchange / max_exchange_doubles
      (atom_pack.h:112,119,127). A grep shows mc/ contains no hand-rolled atom-row
      serializer and no hand-rolled array growth.
    status: verified
    last_checked: 2026-07-14

  - id: ac-002
    summary: Transaction is Frame-shaped — the bare old_xyz/new_xyz arrays are GONE
    type: code
    pass_when: |
      Per the user directive "我们的输入输出尽可能是Frame，所以不应该有 old_xyz *_xyz
      这种东西": mc::Transaction no longer carries the seven parallel scalar vectors
      (moved/old_x/old_y/old_z/new_x/new_y/new_z, transaction.h:32-38). It carries a
      Frame PATCH — the touched atoms as a small cpu::Frame, plus the Box — and a verb
      (Displace | Insert | Delete | Box). A grep over mc/ finds no raw coordinate array.
      Rollback splices the patch Frame back and restores the frame EXACTLY for each verb:
      a rejected insertion restores nlocal with no orphan row; a rejected deletion
      resurrects the atom with every field intact (x/v/f/mass/charge/type/tag/image +
      bonded topology + special lists — which a bare xyz array could never have held);
      a rejected volume move restores BOTH the box cell and the coordinates. Unit tests
      compare the full frame field-by-field against a pre-move copy and assert bitwise
      equality.
    status: verified
    last_checked: 2026-07-14

  - id: ac-003
    summary: muVT and NPT acceptance rules as OutcomePolicy implementations
    type: code
    pass_when: |
      [2026-07-14: the seam itself was already delivered by the driver-handler refactor —
      SamplingDriver and its hardcoded metropolis_accept_ are deleted; OutcomePolicy
      (cpu/driver/outcome_policy.h) is the swappable slot and cpu::Metropolis its NVT
      implementation.] What remains: a muVT policy (N, V, mu, Lambda — read off
      LoopState.frame) and an NPT policy (P*dV - N kB T ln(V'/V)) each implement
      OutcomePolicy and are attached via set_outcome() by the Gcmc / NptMc presets.
      Both draw from the run's shared RNG stream (ctx.mc.rng), like Metropolis does.
      The existing NVT MC results are unchanged (the NIST LJ validation still passes).
    status: verified
    last_checked: 2026-07-14

  - id: ac-004
    summary: Tags stay unique across delete + insert
    type: code
    pass_when: |
      A global next-free-tag allocator exists. After a delete-compaction followed by an
      insertion, no two atoms share a tag. A unit test that performs 1000 random
      insert/delete cycles asserts tag uniqueness throughout. (Without this, the
      identity-tag seeding in Atom::grow / set_capacity at atom.cpp:71-78 and :151-154
      silently produces duplicate tags.)
    status: verified
    last_checked: 2026-07-14

  - id: ac-005
    summary: The anisotropic-box stale-grid bug is fixed
    type: code
    pass_when: |
      CellList::grid_stale() compares all three box edges, not only length(0)
      (cell_list.cpp:131 — box_at_last_build_ is a single scalar today). A unit test that
      changes ONLY the y or z edge triggers a grid rebuild; before the fix it did not,
      silently dropping pairs — the exact failure mode documented at cell_list.cpp:253-266.
    status: verified
    last_checked: 2026-07-14

  - id: ac-006
    summary: Multiple movers have correct semantics via a selector
    type: code
    pass_when: |
      Registering TranslateMove and InsertDeleteMove side by side does NOT fire both every
      sweep. A MoveSelector draws which move type is attempted per sweep with configured
      probabilities, and the per-type attempt/acceptance counts are reported. A unit test
      asserts the drawn move-type frequencies match the configured probabilities.
    status: verified
    last_checked: 2026-07-14

  - id: ac-007
    summary: Regression — GCMC LJ adsorption isotherm vs literature
    type: scientific
    pass_when: |
      regressions/mc/gcmc_lj/ runs muVT Monte Carlo on the Lennard-Jones fluid through the
      public API and reproduces a published <N>/density at a given (mu, V, T) state point,
      with REFERENCE.md citing the source (e.g. the NIST SRSW LJ muVT reference data, or
      Frenkel & Smit's tabulated LJ equation-of-state point). Detailed balance is
      demonstrated. Runs <= 60 s.
    status: pending

  - id: ac-008
    summary: NPT-MC reproduces the equation of state
    type: scientific
    pass_when: |
      NPT Monte Carlo on the LJ fluid at a state point reproduces the published mean
      density for that (P, T), cited in REFERENCE.md. The volume-move acceptance rule
      includes the P*dV and N kB T ln(V'/V) terms.
    status: pending

  - id: ac-009
    summary: The architectural verdict is written down
    type: code
    pass_when: |
      .claude/NOTES.md records the headline finding — that the Frame's "fixed-N assumption"
      was a myth: Atom::grow/set_capacity, comm::pack_exchange/unpack_exchange, the
      migration insert/delete paths, and NeighborList's size_stale trigger already
      implemented variable N; the fixed-N assumption lived ONLY in Frame's constructor and
      in mc::Transaction. Also records that GCMC correctly needs NO new Driver category —
      confirmed twice over: the driver-handler refactor deleted categories entirely, and
      the acceptance rule became an OutcomePolicy — and that trial energy is O(N)
      full-recompute by deliberate MC contract (mc.h run() docs), with energy_of_rows()
      deferred.
    status: verified
    last_checked: 2026-07-14
---

# Acceptance — mc-gcmc-npt

Probe #4. The lesson is reuse: `comm/` already solved variable-N; lift it, don't rebuild it.
