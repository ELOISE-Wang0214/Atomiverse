---
slug: hf-robust-02-canonical-orthogonalization
criteria:
  - id: ac-001
    summary: canonical_orthogonalization unit-correct on synthetic near-singular S
    type: code
    pass_when: |
      tests/cpu/linalg/test_canonical_orthogonalization.cpp passes in the default
      ctest gate: for a synthetic S with one ~1e-12 eigen-direction, X has shape
      n×(n−1), max|XᵀSX − I| < 1e-12, n_dropped == 1, s_min matches the injected
      eigenvalue; tau=0 keeps all n columns; a tau dropping every column throws
      std::runtime_error.
    status: pending
  - id: ac-002
    summary: Healthy-basis no-regression — dsygvd path and all existing values unmoved
    type: code
    pass_when: |
      With default lindep_threshold=1e-8: every pre-existing test under
      tests/cpu/scf/ passes without modification; clean H2/6-31G RHF energy in
      test_scf_lindep.cpp equals the pre-change value to 1e-12 Ha; with
      lindep_threshold=0 the clean-system result is bitwise identical to the
      current dsygvd output.
    status: pending
  - id: ac-003
    summary: Near-degenerate H2+H+ system converges with m=4 and matches PySCF
    type: scientific
    pass_when: |
      For H2 + H⁺ / 6-31G (H at 0,0,0; 0,0,0.7414; 0,0,1e-5 Å; charge=+1),
      RHF with lindep_threshold=1e-8 converges, drops exactly 2 overlap
      eigenvectors (n_mo == 4 of n_ao == 6), and |E − E_ref| < 1e-6 Ha where
      E_ref is a freshly regenerated PySCF 2.13.1 (cart=True) value for this
      exact geometry/charge, recorded with its provenance in
      tests/cpu/scf/test_scf_lindep.cpp.
    status: pending
  - id: ac-004
    summary: Tutorial section exists and its numbers are regression-pinned
    type: docs
    pass_when: |
      docs/tutorial/ch02-rhf.md contains a section covering Roothaan–Hall,
      Löwdin S^{−1/2} vs canonical orthogonalization, and linear-dependence
      thresholds (noting PySCF ≥ 2.13 default 1e-6 vs Atomiverse 1e-8); every
      numeric value quoted there (s_min values, dropped-column count, both
      energies) appears as a pinned constant checked by test_scf_lindep.cpp.
    status: pending
  - id: ac-005
    summary: New public symbols carry Doxygen with units
    type: docs
    pass_when: |
      canonical_orthogonalization, CanonicalOrth, EigenSolver::step_orthogonal,
      and SCFConfig::lindep_threshold each have @brief plus @param/@return
      (units stated; dimensionless where applicable), consistent with the
      repository Doxygen convention.
    status: pending
  - id: ac-006
    summary: Full build and default ctest gate pass with new tests inside the gate
    type: runtime
    pass_when: |
      cmake --build build -j succeeds warning-clean for the touched files and
      ctest --test-dir build --output-on-failure passes with both new tests
      registered in the default (non-E2E) gate, each completing in < 1 s.
    status: pending
---

# Acceptance criteria

- **ac-001/ac-002** are `code`-typed: verifiable by `/mol:impl` running the new unit test and the untouched existing SCF suite.
- **ac-003** is `scientific`: it binds the implementer to regenerate the PySCF 2.13.1 reference for the H2 + H⁺ system (the −1.1267347374 Ha number in the domain basis was produced with a zero-charge ghost center and is **not** directly reusable; only m = 4 carries over).
- **ac-004** enforces the tutorial-regression coupling: docs numbers and test constants must be the same literals.
- **ac-006** protects the ~16 s default ctest budget; neither new test may land in the E2E-disabled tail of `tests/CMakeLists.txt`.

All file paths verified against the repo (branch `dev`); `tests/cpu/linalg/` is a new directory. Key placement confirmations: `X_` lives at `include/atomiverse/cpu/scf/scf_solver.h:99`, the Löwdin build at `src/cpu/scf/scf_solver.cpp:91`, the dsygvd solve at lines 99 and 199–205, and a tutorial-regression precedent exists (`test_cpu_gpw_h2o_tutorial`, `tests/CMakeLists.txt:164`).
