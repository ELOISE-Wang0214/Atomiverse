---
slug: scf-review-03-posthf-dedup
criteria:
  - id: ac-001
    summary: EOM intermediates defined only in eom_ccsd_common; no free build_imds
    type: code
    pass_when: |
      grep over src/cpu/scf shows make_tau and `struct EomImds` DEFINED only in
      src/cpu/scf/internal/eom_ccsd_common.cpp (declared in
      include/atomiverse/cpu/scf/internal/eom_ccsd_common.h); the intermediates
      build is the computing constructor EomImds(const CCSDSolution&) and the
      free function `build_imds` no longer exists ANYWHERE under src/ or
      include/ (this clears its 3 pre-existing no_factory_scrub hits). ccsd.cpp,
      eom_ee_ccsd.cpp, eom_ea_ccsd.cpp, eom_ip_ccsd.cpp contain no local
      definition of these symbols and each #include the new header.
    status: verified
    last_checked: 2026-07-24
  - id: ac-002
    summary: CI engine defined only in ci_engine; slater_determinant untouched
    type: code
    pass_when: |
      choose, enumerate_strings, occ_list, jw_single_sign, diagonal_element,
      single_element, same_spin_double_element, opposite_spin_double_element,
      ci_element, dense_ci_solve, make_rdm12, `struct CiSolution` are DEFINED
      only in src/cpu/scf/internal/ci_engine.cpp (declared in
      include/atomiverse/cpu/scf/internal/ci_engine.h). fci.cpp and casscf.cpp
      contain none and #include ci_engine.h. include/atomiverse/cpu/scf/
      slater_determinant.h is unchanged (git diff empty).
    status: verified
    last_checked: 2026-07-24
  - id: ac-003
    summary: dedup grep-gate scrub registered and green
    type: test
    pass_when: |
      ctest test test_posthf_dedup_scrub (bash tests/architecture/
      posthf_dedup_scrub.sh <repo-root>) exits 0, and fails non-zero if any
      hoisted symbol is defined in any src/cpu/scf/*.cpp outside the two
      internal engine files.
    status: verified
    last_checked: 2026-07-24
  - id: ac-004
    summary: canonical Slater-Condon element pinned to hand value
    type: test
    pass_when: |
      tests/cpu/scf/test_ci_engine.cpp: for a hand-built 2-orbital/2-electron
      case with literal h,g, ci_element returns the analytically computed
      diagonal, single-excitation (with Jordan-Wigner sign), and
      opposite-spin double-excitation matrix elements within 1e-12 of the
      hard-coded expected literals.
    status: verified
    last_checked: 2026-07-24
  - id: ac-005
    summary: make_tau pinned to hand value
    type: test
    pass_when: |
      tests/cpu/scf/test_eom_ccsd_common.cpp: make_tau on a hand case
      reproduces tau = t2 + 0.5*fac*P(ij)P(ab) t1 t1 element-wise within
      1e-12 of hard-coded expected literals (fac=1.0 and fac=0.5).
    status: verified
    last_checked: 2026-07-24
  - id: ac-006
    summary: regression proves FCI == CASCI == literature via unified engine
    type: runtime
    pass_when: |
      regressions/scf/ci-engine-equiv-h2-sto3g/reg_ci-engine-equiv-h2-sto3g
      exits 0: RHF::fci and RHF::casci(2,2) on H2/STO-3G both equal the
      embedded PySCF 2.14.0 full-CI literal -1.137283834488502 Ha (recorded
      in REFERENCE.md with tool/version/command/date) to 1e-6 Ha, AND
      |E_fci - E_casci| < 1e-12; runs in < 60 s; no live third-party call.
    status: verified
    last_checked: 2026-07-24
  - id: ac-007
    summary: existing EOM/CCSD/CI regressions unchanged
    type: runtime
    pass_when: |
      reg_eom-ee-ccsd, reg_eom-ea-ccsd, reg_eom-ip-ccsd, reg_eom-cc2-h2,
      reg_eom-sf-ccsd-o2, reg_ccsd-h2-sto3g, reg_ccsd-t-h2o-sto3g,
      reg_cc2-h2-sto3g, reg_ccsd-grad-h2o-sto3g, reg_casci-h2-sto3g,
      reg_casscf-h2-sto3g, reg_nevpt2-h2-sto3g all exit 0 with their
      pre-existing pinned literals and tolerances (embedded goldens
      unchanged; git diff of each reg_*.cpp shows no numeric edit).
    status: verified
    last_checked: 2026-07-24
  - id: ac-008
    summary: full check + test suite green
    type: runtime
    pass_when: |
      HF-enabled configure (-DATV_BUILD_HF=ON -DATV_BUILD_TESTS=ON) builds,
      ctest --output-on-failure is all green, and
      cmake --build build --target regressions completes with every case
      exiting 0.
    status: verified
    last_checked: 2026-07-24
---

# Acceptance criteria

- **ac-001 / ac-002 (code):** the two hoists physically relocate the duplicated symbols into the two new `internal/` engine files; the source-of-duplication files retain no local copy. ac-002 also enforces the librarian rule that `slater_determinant.h` is not extended.
- **ac-003 (test):** a grep-gate scrub (mirroring `tests/architecture/no_factory_scrub.sh`) makes re-duplication a build failure.
- **ac-004 (test):** pins the *canonical* fci Slater-Condon evaluator against a hand-computable 2-electron matrix element — this is the discovered-bug detector: if the old casscf fork disagreed, this test states the correct value.
- **ac-005 (test):** pins the single hoisted `make_tau` against its analytic definition.
- **ac-006 (runtime):** the spec's regression example — FCI and CASCI(full space) now route the same engine and must both reproduce the cited full-CI literal, proving the merge introduced no numeric change and both entry points agree.
- **ac-007 (runtime):** behavior-preserving guarantee for every downstream method (EOM-EE/EA/IP/CC2/SF, CCSD/(T)/CC2, CASCI/CASSCF/NEVPT2).
- **ac-008 (runtime):** whole-suite green, including unit tiers and the regressions target.
