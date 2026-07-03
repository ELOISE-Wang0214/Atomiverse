---
slug: mc-sampling-02-move
criteria:
  - id: ac-001
    summary: cpu::mc::MC and TranslateMove exist with correct tags, no factories
    type: code
    pass_when: |
      include/atomiverse/cpu/mc/mc.h defines ATV::cpu::mc::MC subclassing
      SamplingDriver with a set_mover() using set_singular_default_ and a
      validating MC{T,seed} constructor; include/atomiverse/cpu/mc/translate_move.h
      defines TranslateMove with role()==ComponentRole::Mover and
      capabilities()==Capability::None and a validating δ constructor; no
      free factory function (make_*/resolve_*) appears in either file.
    status: pending
  - id: ac-002
    summary: ZeroEnergy head-bound at sampling energy phase order -999
    type: code
    pass_when: |
      The MC constructor binds cpu::ZeroEnergy at sampling::Phase::energy with
      order -999 as the head of that phase, mirroring MD's pre_force -999 bind.
    status: pending
  - id: ac-003
    summary: Bootstrap e_cur_ equals fresh full-Frame energy to 1e-12
    type: runtime
    pass_when: |
      test_cpu_mc_move asserts |e_cur_ - fresh_full_energy(initial_frame)| <= 1e-12
      after build()/first evaluation, and the test passes.
    status: pending
  - id: ac-004
    summary: Rejected move restores positions exactly and leaves e_cur_ unchanged
    type: runtime
    pass_when: |
      After a forced rejection, the test asserts every frame position matches the
      pre-move value to <= 1e-15 (restored from Transaction.old_pos) and e_cur_ is
      unchanged; the test passes.
    status: pending
  - id: ac-005
    summary: Accepted move updates e_cur_ to full recompute to 1e-12
    type: runtime
    pass_when: |
      After a forced acceptance, the test asserts |e_cur_ - fresh_full_energy(new_state)|
      <= 1e-12; the test passes.
    status: pending
  - id: ac-006
    summary: Metropolis criterion is exact for downhill and T to 0
    type: runtime
    pass_when: |
      The test asserts a ΔE <= 0 move is accepted unconditionally (probability 1)
      and that at T -> 0 only non-uphill moves are accepted; the test passes.
    status: pending
  - id: ac-007
    summary: Same seed reproduces sequence and trajectory; different seed differs
    type: runtime
    pass_when: |
      The test asserts two MC runs with identical seed produce identical
      accept/reject sequences and identical final positions, and that a different
      seed produces a different sequence; the test passes.
    status: pending
  - id: ac-008
    summary: Acceptance ratio decreases monotonically with δ and mid-δ is in band
    type: runtime
    pass_when: |
      The test asserts accepted-fraction is monotonically non-increasing across an
      increasing δ sweep and that a mid-range δ falls within a loose ~20-60% band
      (trend + band only, no point value); the test passes.
    status: pending
  - id: ac-009
    summary: build() throws on missing Mover or missing Energy evaluator
    type: runtime
    pass_when: |
      The test asserts cpu::mc::MC::build() throws when the Mover is omitted, and
      throws when no Capability::Energy evaluator is registered (via
      SamplingDriver::validate_components); a fully-assembled MC builds without
      throwing; the test passes.
    status: pending
  - id: ac-010
    summary: MC sources and test registered in CMake and full suite passes
    type: runtime
    pass_when: |
      src/cpu/mc/translate_move.cpp is in the atomiverse source list and
      test_cpu_mc_move is registered in tests/CMakeLists.txt; a clean configure +
      build succeeds and ctest passes with test_cpu_mc_move green.
    status: pending
---

# Acceptance criteria

- ac-001 / ac-002 (code): structural — the new driver and Mover exist with the correct role/capability tags, validating constructors, no factory functions, and the ZeroEnergy head-bind that guarantees a zeroed energy accumulator per trial.
- ac-003 through ac-008 (runtime): the core-correctness unit tests — energy bootstrap, exact rollback, energy caching on accept, Metropolis exactness, determinism, and the acceptance-ratio trend/band sanity check.
- ac-009 (runtime): build-time validation surfaced through the concrete driver.
- ac-010 (runtime): CMake wiring plus a green full build + test suite.
