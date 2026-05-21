---
slug: phase4-vibrations-02-frequencies
criteria:
  - id: ac-001
    summary: nist constants and IUPAC-2013 masses present, no hard-coded 5140 literal
    type: code
    pass_when: |
      include/atomiverse/cpu/vib/nist_constants.h defines au2hz from the
      CODATA-2014 two-factor form (HBAR, E_MASS, BOHR_SI, ATOMIC_MASS,
      LIGHT_SPEED_SI); grep across include/atomiverse/cpu/vib finds no
      hard-coded 5140.x wavenumber literal; atomic_masses.h contains the
      values 1.008, 12.011, 14.007, 15.999 for H/C/N/O.
    status: pending
  - id: ac-002
    summary: MassWeighting divides Hessian by sqrt(m_i m_j) correctly
    type: runtime
    pass_when: |
      test_cpu_vib_mass_weighting passes: a hand-built 6x6 Hessian on a
      2-atom system with known masses produces H^mw_{iα,jβ} equal to
      H_{iα,jβ}/sqrt(m_i m_j) to machine precision, and mass_amu throws
      on an unsupported atomic number.
    status: pending
  - id: ac-003
    summary: TR projection yields exact 3N-6 / 3N-5 genuine mode count
    type: runtime
    pass_when: |
      test_cpu_vib_tr_projection passes: at an optimized H2O geometry the
      projector annihilates all 6 TR vectors and diagonalization yields
      exactly 3 finite modes; for a linear molecule (CO2 or HCN)
      is_linear is true and exactly 3N-5 modes are returned.
    status: pending
  - id: ac-004
    summary: imaginary modes reported as negative wavenumbers and counted
    type: runtime
    pass_when: |
      test_cpu_vib_tr_projection (or test_cpu_vib_frequencies_vs_pyscf)
      passes: a non-stationary or saddle-point geometry yields exactly
      one mode with a negative wavenumber and HessianResult.n_imaginary
      equals 1.
    status: pending
  - id: ac-005
    summary: frequencies match PySCF harmonic_analysis within 5 cm-1
    type: runtime
    pass_when: |
      test_cpu_vib_frequencies_vs_pyscf passes: the cpu::vib::Hessian
      driver on H2O, NH3, CH4 (RHF or PBE / 6-31G) produces frequencies
      whose max deviation over the projected 3N-6 / 3N-5 modes is
      <= 5 cm-1 versus the checked-in PySCF analytic-Hessian references
      in tests/cpu/vib/reference/pyscf_frequencies.json.
    status: pending
  - id: ac-006
    summary: full ctest suite stays green
    type: runtime
    pass_when: |
      ctest --test-dir build --output-on-failure passes with zero
      failures after the new vib components and default_bindings
      extension land.
    status: pending
---

# Acceptance criteria

- **ac-001 (code).** Static check that the wavenumber conversion is derived from `nist`-mirrored CODATA-2014 constants, never a magic `5140.x` literal, and that the IUPAC-2013 mass table values are present. Guards against rounding drift (~0.1 cm⁻¹) and isotope-averaging mistakes (~0.7 cm⁻¹).
- **ac-002 (runtime).** Unit-level verification of the `MassWeighting` transform and the throwing mass lookup on an isolated hand-built Hessian.
- **ac-003 (runtime).** Verifies the translation/rotation projector and `_get_rotor_type` decision: exactly 3N−6 modes for nonlinear H₂O, 3N−5 for a linear molecule.
- **ac-004 (runtime).** Verifies imaginary-mode handling — negative-wavenumber convention and `n_imaginary` diagnostic — on a non-stationary or saddle geometry.
- **ac-005 (runtime, domain validation).** The binding physical-correctness bar: numerical-Hessian frequencies within 5 cm⁻¹ of PySCF analytic-Hessian ground truth on H₂O, NH₃, CH₄. Wall-time is a non-binding note; benzene/formaldehyde optional.
- **ac-006 (runtime).** Regression gate — no existing test broken by the new components or the `default_bindings()` extension.
