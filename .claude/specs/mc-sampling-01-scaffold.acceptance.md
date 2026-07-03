---
slug: mc-sampling-01-scaffold
criteria:
  - id: ac-001
    summary: sampling::Phase enum has correct values, order, and COUNT sentinel
    type: code
    pass_when: |
      include/atomiverse/facade/stage.h defines namespace ATV::sampling with
      enum class Phase { propose, energy, end_of_step, COUNT } in that order
      and no `accept` value; a doc comment states accept/reject is driver
      outer-loop control, not a phase.
    status: pending
  - id: ac-002
    summary: phase_names(sampling::Phase) matches the enum order and arity
    type: runtime
    pass_when: |
      test_mc_scaffold asserts phase_names(sampling::Phase{}) == {"propose",
      "energy","end_of_step"} and static_cast<int>(sampling::Phase::COUNT)==3.
    status: pending
  - id: ac-003
    summary: ComponentRole::Mover exists on the structural axis
    type: code
    pass_when: |
      include/atomiverse/facade/component_role.h declares ComponentRole::Mover
      after Barostat with the "Proposes a trial state change ..." doc comment.
    status: pending
  - id: ac-004
    summary: cpu::mc::Transaction is a documented pure-data struct with no factory
    type: code
    pass_when: |
      include/atomiverse/cpu/mc/transaction.h defines struct Transaction in
      namespace ATV::cpu::mc holding LocalInt indices + StorageReal old/new
      x,y,z, has Doxygen @brief with units, a validating constructor, and no
      make_/resolve_/factory free function; grep of the repo shows Transaction
      is never a member of cpu::Frame.
    status: pending
  - id: ac-005
    summary: Transaction round-trips move data and rejects mismatched lengths
    type: runtime
    pass_when: |
      test_mc_scaffold constructs a Transaction with N moved atoms, reads back
      identical indices and old/new positions, and asserts the constructor
      throws on mismatched index/position array lengths.
    status: pending
  - id: ac-006
    summary: cpu::Context carries proposal + rng handles defaulting to null
    type: code
    pass_when: |
      include/atomiverse/cpu/component/context.h forward-declares mc::Transaction
      and adds `mc::Transaction* proposal = nullptr;` and
      `std::mt19937_64* rng = nullptr;` with pair_list-style doc comments; a
      default cpu::Context has both null (asserted in test_mc_scaffold).
    status: pending
  - id: ac-007
    summary: SamplingDriver base compiles with bootstrap-seeded Metropolis run() and seam
    type: code
    pass_when: |
      include/atomiverse/cpu/driver/sampling_driver.h defines
      cpu::SamplingDriver : public CPUDriver<sampling::Phase> owning
      proposal_/rng_/e_cur_; run() seeds e_cur_ from an initial full-Frame
      energy before the sweep loop and applies inline min(1,exp(-dE/(kB*T)))
      with commit/rollback; a `// ponytail:` seam comment marks the swappable
      acceptance-policy point; the header builds under -DATV_BUILD_TESTS=ON.
    status: pending
  - id: ac-008
    summary: validate_components() throws without a Mover and without an Energy writer
    type: runtime
    pass_when: |
      test_mc_scaffold shows SamplingDriver::validate_components() throws when
      no ComponentRole::Mover is registered, and throws when a stub Mover is
      present but no component advertises Capability::Energy.
    status: pending
  - id: ac-009
    summary: full check + test suite is green with the dormant scaffold
    type: runtime
    pass_when: |
      cmake --build then ctest (with -DATV_BUILD_TESTS=ON) passes including the
      newly registered test_mc_scaffold target.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** pin the `sampling::Phase` skeleton: static shape (values, order, no `accept`) plus the runtime `phase_names` arity/order check.
- **ac-003** pins the `Mover` structural role.
- **ac-004 / ac-005** pin `Transaction` as documented pure data with a validating (non-factory) constructor, Frame-purity preserved (grep), and a runtime round-trip + validation check.
- **ac-006** pins the two new `Context` handles and their null defaults.
- **ac-007** pins the `SamplingDriver` base: bootstrap-seed invariant, inline Metropolis commit/rollback, the `// ponytail:` seam, and that the dormant header compiles.
- **ac-008** pins the two-branch `validate_components()` gate using test-only stubs.
- **ac-009** is the overall build + fast-gate `ctest` green bar.
