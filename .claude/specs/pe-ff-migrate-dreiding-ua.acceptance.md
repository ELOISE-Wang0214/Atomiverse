---
slug: pe-ff-migrate-dreiding-ua
criteria:
  - id: ac-001
    summary: All FF constant declarations in tree match Dreiding-UA values
    type: code
    pass_when: |
      `grep -RnE "kAngleTheta0\s*=\s*112\.7|kAngleK_kcal\s*=\s*58\.35|2\.0058|-2\.5128|0\.2710|6\.2900"
      benchmarks/ tests/` returns zero matches (the OPLS-like
      literals are completely gone from the codebase, outside
      `SOURCES.md` historical footnotes and existing run-output
      CSVs). `grep -RnE "kAngleTheta0\s*=\s*109\.5|kAngleK_kcal\s*=\s*60|1\.73|-4\.49|0\.776|6\.99"
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
    summary: Full ctest cpu+mpi suite stays green
    type: runtime
    pass_when: |
      `cmake --build build -j && ctest --test-dir build
      --output-on-failure` exits 0; total passing count is
      ≥ pre-migration baseline (87 — recorded by
      `parallel-md-finish-06-close-out` on 2026-05-19) plus
      the new pin test, so ≥ 88. No PE test regresses
      (`test_pe_*`, `test_mpi_pe_*`, `test_lj_exclusions`,
      `test_nose_hoover_t_ramp_n2` all PASS).
    status: pending

  - id: ac-004
    summary: pe_chain.mol regenerated for new half-apex (35.25°)
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/pe_chain.mol`
      contains a position table where atom 2's `y` coordinate
      is within ±1e-4 of `1.248645` (= 1.529 · cos(35.25°))
      and atom 2's `z` coordinate is within ±1e-4 of
      `0.882455` (= 1.529 · sin(35.25°)). The docstring at
      the top references `θ0 = 109.5°` (not `112.7°`).
    status: pending

  - id: ac-005
    summary: Small-scale lammps.in updated and re-runs cleanly
    type: runtime
    pass_when: |
      `grep -E "^\s*angle_coeff" benchmarks/md/reference/polyethylene/lammps.in`
      shows `1 60 109.5` (not `58.35 112.7`); same for
      `dihedral_coeff` showing the new K1..K4 set.
      `cd benchmarks/md/reference/polyethylene &&
      /opt/homebrew/bin/lmp_mpi -in lammps.in -log
      /tmp/lammps_smoke.log` exits 0 in < 30 s.
    status: pending

  - id: ac-006
    summary: Re-committed thermo + stress CSVs still pass parallel-md-finish-02 shape checks
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/thermo_reference.csv`
      and `stress_strain_reference.csv` were updated by the
      `lmp_mpi` re-run in ac-005. Each file has ≥ 11 lines
      with the documented header; final T in `thermo_reference.csv`
      ∈ [80, 150] K; final strain in `stress_strain_reference.csv`
      ≥ 0.3.
    status: pending

  - id: ac-007
    summary: SOURCES.md marks FF discrepancy resolved
    type: docs
    pass_when: |
      `benchmarks/md/reference/polyethylene/SOURCES.md`
      "Known discrepancies" subsection contains a "resolved"
      annotation dated 2026-05-19 or later, with a link / SHA
      pointer to the commit that closed this spec. The
      original discrepancy table is preserved as a historical
      footnote.
    status: pending

  - id: ac-008
    summary: msu-pe-validation Phase D status flipped to ready
    type: docs
    pass_when: |
      `.claude/specs/msu-pe-validation.md` status snapshot row
      for Phase D no longer contains "blocked on FF
      migration"; new text references the Dreiding-UA migration
      commit SHA and reads `🟡 ready to run` (or equivalent
      wording indicating the spec is no longer FF-blocked).
    status: pending

  - id: ac-009
    summary: Single-chain Atomiverse-vs-LAMMPS energy decomposition agrees to 1e-6 eV per term
    type: scientific
    evaluator_hint: human
    pass_when: |
      On the same all-trans C20H42 united-atom chain geometry,
      Atomiverse's `HarmonicBond + HarmonicAngle +
      MultiHarmonicDihedral` evaluator and LAMMPS running
      `lammps_msu_paper.in` report `E_bond`, `E_angle`,
      `E_dihedral` agreeing to 1e-6 eV absolute per term.
      This is the fundamental "we run the MSU paper's FF"
      check; requires a small `tools/single_chain_ff_check.py`
      (or equivalent) which is NOT mandated by this spec but
      is the natural follow-up. Closed by manual attestation
      (`/mol:close --manual`) after the operator has run the
      comparison externally.
    status: pending

  - id: ac-010
    summary: Potential-abstraction audit note committed
    type: docs
    pass_when: |
      `.claude/notes/potential-abstraction-audit.md` exists.
      Its body answers the six audit questions in the spec's
      "Architecture audit task" section, with one row per
      question per backend (CPU / CUDA / torch). Each row
      records: the current state (one paragraph), the gap
      severity (low / medium / high), and one concrete
      refactor proposal. At the bottom, the note names a
      proposed follow-up spec slug (e.g.
      `potential-base-class-refactor`) and lists the smallest
      concrete refactor that would close the most severe
      gap. The note references this spec's commit SHA.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-004 / ac-005 / ac-006 / ac-007 / ac-008** are static checks /mol:impl can verify directly: file existence, content regex match, exit codes. These flip to `verified` automatically at /mol:impl close-out.
- **ac-003** is a runtime gate: the full `ctest` suite. The pre-migration baseline of 87 passing tests is the documented floor; the new pin test (ac-002) brings the floor to 88.
- **ac-009** is the **scientific** load-bearing acceptance: it's the test that closes the loop on "Atomiverse now runs the MSU paper's FF". It cannot be verified by /mol:impl alone — it requires running both engines on a known geometry and comparing per-term energies. Either a small `tools/single_chain_ff_check.py` is added (in a follow-up spec or commit) and run, OR the operator manually does the comparison and closes via `/mol:close --manual`. The criterion stays `pending` until then.
