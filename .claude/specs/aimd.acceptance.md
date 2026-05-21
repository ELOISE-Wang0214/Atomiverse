---
slug: aimd
criteria:
  - id: ac-001
    summary: MD::set_qm registers an SCFEvaluator as the force source
    type: runtime
    pass_when: |
      tests/cpu/md/test_aimd_set_qm.cpp builds an MD via set_qm(SCFEvaluator{...})
      + set_integrator, calls build() and run() on an H2 STO-3G fixture, and the
      resulting frame.energy is finite and inside the bonded H2 well (< -25 eV);
      the test passes under ctest.
    status: pending
  - id: ac-002
    summary: set_pair and set_qm are mutually exclusive
    type: runtime
    pass_when: |
      tests/cpu/md/test_aimd_set_qm.cpp asserts that set_pair-then-set_qm and
      set_qm-then-set_pair each throw std::invalid_argument, and that a second
      set_qm with a different component throws; all assertions pass under ctest.
    status: pending
  - id: ac-003
    summary: ASPCPredictor::predict uses the Kolafa B_j coefficients
    type: runtime
    pass_when: |
      tests/cpu/scf/test_aspc_predictor.cpp computes ASPCPredictor coefficients
      for nvec = 2..6 and asserts each matches the verified table
      (nvec=3 -> {2.5,-2,0.5}, etc.) to 1e-12 and sums to 1; the test passes
      under ctest.
    status: pending
  - id: ac-004
    summary: NVE AIMD via set_qm conserves total energy
    type: runtime
    pass_when: |
      An H2 STO-3G NVE AIMD run through set_qm (dt = 0.5 fs, 200 steps) keeps the
      E_total = energy + kinetic peak-to-peak envelope below 0.10 eV and the
      energy inside the bonded H2 well, asserted in
      tests/cpu/md/test_aimd_set_qm.cpp and passing under ctest.
    status: pending
  - id: ac-005
    summary: ASPC warm-start cuts mean SCF iterations per step
    type: runtime
    pass_when: |
      In tests/cpu/scf/test_aimd_conservation.cpp, an AIMD trajectory with
      SCFEvaluator::predictor() enabled records a strictly lower mean SCF
      iterations/step (measured past predictor fill) than the same trajectory
      with cold-start; the assertion passes under ctest.
    status: pending
  - id: ac-006
    summary: ASPC docstring no longer misattributes ASPC to Car-Parrinello
    type: code
    pass_when: |
      grep of include/atomiverse/cpu/scf/aspc_predictor.h finds no
      "Car-Parrinello" / "Kuhne" / "Kühne" string and no "Lagrange" in any
      ASPCPredictor doc comment; the header describes ASPC as a
      Born-Oppenheimer SCF initial-guess accelerator.
    status: pending
  - id: ac-007
    summary: Full test suite stays green
    type: runtime
    pass_when: |
      `cmake --build build` succeeds and `ctest --test-dir build` reports all
      tests passing, including the two new test executables.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** verify the `set_qm` ergonomic surface and the pair-vs-QM
  force-source mutual exclusion — the core "first-class composition" deliverable.
- **ac-003** is the binding correctness bar for the ASPC fix: the coefficient
  table is the externally verifiable signature of faithful Kolafa ASPC (vs the
  current Lagrange coefficients).
- **ac-004** is the defining BOMD physics validation — energy conservation over a
  short NVE trajectory, bound chosen consistent with the existing
  `test_aimd_conservation.cpp` envelope test and literature drift figures. It is
  a deterministic ctest assertion (peak-to-peak envelope < 0.10 eV), hence
  `runtime`.
- **ac-005** confirms the warm-start actually accelerates the SCF, the practical
  payoff of a correct cross-step predictor.
- **ac-006** is the documentation correctness fix (ASPC is a BO initial-guess
  accelerator, not second-gen Car-Parrinello); verified by grep, hence `code`.
- **ac-007** guards against regressions across the full suite, including the
  existing ASPC-vs-fresh-single-point test which the corrected `B_j`
  coefficients must not break.
