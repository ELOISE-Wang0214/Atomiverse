---
slug: pe-ff-migrate-dreiding-ua
criteria:
  - id: ac-001
    summary: OPLS-like constants gone from the codebase; Dreiding-UA constants in place
    type: code
    pass_when: |
      `grep -RnE "kAngleTheta0\s*=\s*112\.7|kAngleK_kcal\s*=\s*58\.35|2\.0058|-2\.5128|0\.2710|6\.2900"
      benchmarks/ tests/` returns zero matches (the OPLS-like
      literals are completely gone from the codebase outside
      `SOURCES.md` historical context).
      `grep -RnE "kAngleTheta0\s*=\s*109\.5|kAngleK_kcal\s*=\s*60|1\.73|-4\.49|0\.776|6\.99"
      benchmarks/ tests/` returns ≥ 30 matches (each of the 15
      affected files contributes the three updated constants).
    status: pending

  - id: ac-002
    summary: New pin test asserts Dreiding-UA values
    type: code
    pass_when: |
      `tests/cpu/md/test_pe_ff_constants.cpp` exists and is
      registered in `tests/CMakeLists.txt`. `ctest -R
      test_pe_ff_constants` reports 1/1 passing. Its assertions
      cover all three migrated constants by direct comparison
      against the literal Dreiding-UA values.
    status: pending

  - id: ac-003
    summary: Full ctest suite stays green with no PE-test regression
    type: runtime
    pass_when: |
      `cmake --build build -j && ctest --test-dir build
      --output-on-failure` exits 0. Total passing count ≥ 90
      (89 post-`pe-potential-abstraction` baseline + the new
      `test_pe_ff_constants` pin). No PE test regresses
      (`test_pe_*`, `test_mpi_pe_*`, `test_lj_exclusions`,
      `test_nose_hoover_t_ramp_n2` all PASS).
    status: pending

  - id: ac-004
    summary: SOURCES.md FF discrepancy marked resolved
    type: docs
    pass_when: |
      `benchmarks/md/reference/polyethylene/SOURCES.md`
      "Atomiverse FF discrepancy" subsection contains a
      "resolved 2026-05-19" annotation (or later date) and the
      commit SHA that closed this spec. The original mismatch
      table is preserved as a historical footnote.
    status: pending

  - id: ac-005
    summary: msu-pe-validation Phase D status flipped to ready
    type: docs
    pass_when: |
      `.claude/specs/msu-pe-validation.md` status snapshot row
      for Phase D no longer contains "blocked on FF
      migration"; new text references this spec's commit SHA
      and reads `🟡 ready to run` (or equivalent wording
      indicating the FF blocker is cleared).
    status: pending

  - id: ac-006
    summary: Single-chain Atomiverse-vs-LAMMPS energy decomposition agrees to 1e-6 eV per term
    type: scientific
    evaluator_hint: human
    pass_when: |
      On the same all-trans C20H42 united-atom chain geometry,
      Atomiverse's `BondHarmonic + AngleHarmonic +
      DihedralMultiHarmonic` evaluator and LAMMPS running
      `lammps_msu_paper.in` report `E_bond`, `E_angle`,
      `E_dihedral` agreeing to 1e-6 eV absolute per term.
      This is the fundamental "Atomiverse now runs the MSU
      paper's FF" check. NOT mandated by this spec to be
      automated — it's the natural follow-up that
      `msu-pe-validation.md` Phase A.2 covers. Closed by
      manual attestation (`/mol:close --manual`) after the
      operator runs the comparison externally.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** are static + runtime checks `/mol:impl`
  can flip directly. Together they prove the constants are correctly
  swapped and nothing breaks.
- **ac-004 / ac-005** are documentation updates downstream consumers
  need to see (audit-trail + parent-spec status sync). Per `/mol:impl`
  rule they stay `pending` for human attestation; observed-met values
  are recorded on each criterion.
- **ac-006** is the **scientific** load-bearing gate — the closing
  proof that Atomiverse and MSU now agree on the physics. It's not
  automated by this spec because (a) the constant swap itself is
  obviously correct under any test that exercises the kernels, and (b)
  the cross-engine numerical comparison belongs to `msu-pe-validation`
  Phase A.2 / D. Closed by `/mol:close --manual` after the operator
  verifies.
