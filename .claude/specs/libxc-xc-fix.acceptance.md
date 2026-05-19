---
slug: libxc-xc-fix
criteria:
  - id: ac-001
    summary: libxc integrates into CMake under ATV_USE_LIBXC
    type: code
    pass_when: |
      Configuring with -DATV_USE_LIBXC=ON succeeds, defines ATV_HAS_LIBXC,
      and links Libxc::xc into the atomiverse target; configuring with
      -DATV_USE_LIBXC=OFF builds cleanly with no libxc references.
    status: pending

  - id: ac-002
    summary: Point-level XC parity native vs libxc to 1e-12
    type: scientific
    pass_when: |
      `ctest --test-dir build -R test_xc_vs_libxc --output-on-failure`
      passes; every (exc, vrho_alpha, vrho_beta, vsigma_aa, vsigma_ab,
      vsigma_bb) tuple at the documented test points (LSDA + PBE; un-,
      partially, fully polarized; low and high density) agrees to
      |native - libxc| < 1e-12.
    status: pending

  - id: ac-003
    summary: Lebedev angular grid integrates spherical harmonics exactly
    type: scientific
    pass_when: |
      `ctest -R test_grid` passes; for grids {6, 110, 302}, Σ w_i |Y_l^m|^2
      equals 1.0 to within 1e-14 for every (l, m) up to the grid's stated
      exactness degree.
    status: pending

  - id: ac-004
    summary: Mura-Knowles radial grid integrates Gaussian moments
    type: scientific
    pass_when: |
      `ctest -R test_grid` passes; numerical vs analytic
      ∫ r^{2+n} exp(-α r²) dr agrees to relative error < 1e-10 for
      n in {0,2,4,6}, α in {0.1, 1.0, 10.0}.
    status: pending

  - id: ac-005
    summary: Bragg-Slater table covers Z=1..36 and throws otherwise
    type: code
    pass_when: |
      `bragg_slater_radius(Z)` returns the published value for every
      Z in 1..36 (He=0.31, not 0.25); for Z outside the table it throws
      std::invalid_argument. No silent 1.0 Å fallback path remains in
      `src/cpu/scf/grid.cpp`.
    status: pending

  - id: ac-006
    summary: UKS-LSDA total energy matches PySCF after the fix
    type: scientific
    pass_when: |
      `ctest -R test_uks_solve` passes with: |E(H) - E_PySCF| < 1e-6 Ha,
      |E(He)| < 1e-5 Ha, |E(Li)| < 1e-4 Ha, |E(H2,eq)| < 1e-5 Ha,
      and every point of the H2 PES scan within 1e-4 Ha.
    status: pending

  - id: ac-007
    summary: Native XC path agrees with libxc path end-to-end
    type: scientific
    pass_when: |
      Running `test_uks_solve` twice — once with `evaluate_xc_libxc`
      forced and once with native `evaluate_xc` — yields total energies
      that agree to within 1e-12 Ha on every test system.
    status: pending

  - id: ac-008
    summary: Full test suite passes with libxc on and off
    type: runtime
    pass_when: |
      `ctest --test-dir build --output-on-failure` reports 100% pass
      under both `-DATV_USE_LIBXC=ON` and `-DATV_USE_LIBXC=OFF`
      configurations.
    status: pending
---

# Acceptance criteria

These criteria are the binding bar for `libxc-xc-fix`. The scientific
gates (ac-002, ac-003, ac-004, ac-006, ac-007) require domain
validation by the `scientist` evaluator and cannot be discharged by
a green build alone. ac-001, ac-005, ac-008 are code/runtime gates
discharged by `/mol:impl` directly.
