---
slug: ks-xc-grid-response-gradient
criteria:
  - id: ac-001
    summary: becke_weight_derivative matches finite-difference of partition weights
    type: code
    pass_when: |
      tests/cpu/scf/test_ks_grid_response_gradient.cpp builds and the
      per-point case asserting becke_weight_derivative agrees with a central
      difference of becke_partition_weights_for_test (h=1e-4 Bohr) to
      <= 1e-6 passes.
    status: pending
  - id: ac-002
    summary: per-point dw_g/dR_A sums to zero over atoms (rigid-shift identity)
    type: code
    pass_when: |
      the test_ks_grid_response_gradient.cpp case asserting
      sum_A dw_g/dR_A == 0 per grid point and axis to <= 1e-10 on a
      multi-atom geometry passes.
    status: pending
  - id: ac-003
    summary: net translational XC force vanishes without projection
    type: runtime
    pass_when: |
      with SCFConfig::xc_grad_project_translation = false, the
      test_ks_grid_response_gradient.cpp case asserting
      ||sum_A dE_xc/dR_A|| per axis on water PBE/def2-SVP is <= 1e-8 Ha/Bohr
      passes under ctest.
    status: pending
  - id: ac-004
    summary: analytic KS gradient matches moving-grid central-difference FD
    type: runtime
    pass_when: |
      the test_ks_grid_response_gradient.cpp moving-grid FD case (grid
      regenerated per displaced geometry, h=1e-3 Bohr, SCF tol 1e-9) asserts
      per-component agreement <= 1e-5 Ha/Bohr for water PBE/def2-SVP and
      PBE/6-31G* (RKS and UKS), and passes under ctest.
    status: pending
  - id: ac-005
    summary: translational-invariance projection is gated behind a debug flag
    type: code
    pass_when: |
      grep of src/cpu/scf/gradient.cpp shows the projection block formerly at
      lines 447-469 is conditioned on SCFConfig::xc_grad_project_translation,
      and SCFConfig declares that flag with default false in
      include/atomiverse/cpu/scf/scf_solve.h.
    status: pending
  - id: ac-006
    summary: existing dshell-gradient tests and full ctest suite stay green
    type: runtime
    pass_when: |
      ctest --test-dir build --output-on-failure passes with zero failures,
      including all cases in test_ks_dshell_gradient.cpp.
    status: pending
---

# Acceptance criteria

- **ac-001** — Unit-level fidelity of the ported weight derivative against a numerical
  oracle built from the existing test accessor `becke_partition_weights_for_test`. Catches
  sign and chain-rule errors in the `grids_response_becke` port before any SCF runs.
- **ac-002** — The rigid-shift identity `sum_A dw_g/dR_A = 0` must hold per grid point; it is
  a structural property independent of any reference data and isolates scatter/normalisation
  bugs in `dpbecke`.
- **ac-003** — The core first-principles diagnostic: with the projection shortcut disabled,
  the net translational XC force must vanish on its own. No external reference data needed.
- **ac-004** — End-to-end physics validation against a moving-grid finite-difference oracle,
  covering LDA and GGA, RKS and UKS, two basis sets. The `* rho` convention factor and the
  prior-attempt garbage-error class are caught here.
- **ac-005** — Confirms the projection shortcut is removed from the default path and only
  reachable via an explicit debug flag, so the production gradient is first-principles.
- **ac-006** — No-regression gate: the predecessor fixed-grid tests and the whole suite
  remain green.
