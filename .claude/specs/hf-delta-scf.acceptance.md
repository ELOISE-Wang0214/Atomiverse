---
slug: hf-delta-scf
criteria:
  - id: ac-001
    summary: Anchor cation UHF energies and <S2> match pinned PySCF (default gate)
    type: scientific
    pass_when: |
      test_cpu_scf_delta_scf passes in the default ctest run: for He+, H2+,
      H2O+ at STO-3G (neutral geometry, charge=+1, n_alpha=n_beta+1) the
      converged scf_solve energy differs from the pinned PySCF 2.13.1 value in
      tests/cpu/scf/delta_scf_reference.h by <= 1e-7 Ha and
      (spin_contamination + 0.75) differs from the pinned <S2> by <= 1e-4.
    status: pending
  - id: ac-002
    summary: Full 5x2 sweep matches pinned PySCF (E2E tier)
    type: scientific
    pass_when: |
      val_scf_delta_scf_sweep (benchmarks/qc/rhf/validation,
      -DATV_BUILD_BENCHMARKS=ON) passes: all 10
      (system, basis) pairs from delta_scf_reference.h satisfy the 1e-7 Ha
      energy and 1e-4 <S2> tolerances.
    status: pending
  - id: ac-003
    summary: Variational bound dSCF <= -eps_HOMO asserted for every pair
    type: code
    pass_when: |
      Both test targets assert I_dSCF <= -eps_HOMO programmatically (anchors in
      test_delta_scf.cpp; all 10 pairs in test_delta_scf_sweep.cpp), and the
      assertions pass.
    status: pending
  - id: ac-004
    summary: He+/H2+ STO-3G equality with Koopmans to 1e-9 Ha
    type: code
    pass_when: |
      test_cpu_scf_delta_scf asserts |I_dSCF - (-eps_HOMO)| <= 1e-9 Ha for He
      and H2 at STO-3G, and the assertion passes.
    status: pending
  - id: ac-005
    summary: Warm-start and cold-start converge to the same cation solution
    type: code
    pass_when: |
      test_cpu_scf_delta_scf shows, for the anchor systems, that the
      SCFInitialGuess-seeded solve (P_alpha = P_beta = P_neutral/2) and the
      default Hcore-guess solve agree to <= 1e-8 Ha, guarding against
      excited/hole-shifted solutions.
    status: pending
  - id: ac-006
    summary: Example 12_delta_scf builds and prints the comparison table
    type: runtime
    pass_when: |
      With -DATV_BUILD_RHF=ON, target 12_delta_scf compiles and, when run,
      prints one row per (system, basis) pair containing -eps_HOMO [eV],
      I_dSCF [eV], the relaxation difference, and <S2>, using kHartree2eV from
      core/units.h (no inline conversion constants).
    status: pending
  - id: ac-007
    summary: Tutorial dSCF section present and numbers match pinned references
    type: docs
    pass_when: |
      docs/tutorial/ch03-open-shell.md contains a dSCF section (English) with
      the two-line variational argument, spin bookkeeping (n_alpha = n_beta+1,
      M=2), the SCFInitialGuess warm-start pattern, spin-contamination and
      excited-solution caveats, the He/H2 equality teaching point, the
      qualitative vertical-vs-experiment comparison (H2O 12.6 eV, CH4 13.6 eV),
      and a cross-reference to the ch02 Koopmans section; every numeric value
      quoted appears verbatim in tests/cpu/scf/delta_scf_reference.h.
    status: pending
  - id: ac-008
    summary: Default ctest gate respected; sweep only under E2E flag
    type: runtime
    pass_when: |
      Default configure runs test_cpu_scf_delta_scf but
      reports test_cpu_scf_delta_scf_sweep as Disabled (listed in
      ATV_E2E_TESTS in tests/CMakeLists.txt), and the full suite
      `ctest --test-dir build --output-on-failure` passes in both
      configurations.
    status: pending
---

# Acceptance criteria

ac-001/ac-002 are the scientific regression against the pinned PySCF 2.13.1
values embedded in the spec (Domain basis table); ac-003/ac-004/ac-005 are the
physics-invariant code assertions (variational bound, exact-in-basis equality,
excited-solution guard); ac-006 is the runnable teaching deliverable; ac-007
ties the tutorial prose to engine-produced numbers; ac-008 protects the ~16 s
default ctest gate.
