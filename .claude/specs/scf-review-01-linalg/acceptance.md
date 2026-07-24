---
slug: scf-review-01-linalg
criteria:
  - id: ac-001
    summary: Shared linalg primitives replace all in-file duplicates
    type: code
    pass_when: |
      grep across src/cpu/scf finds zero local `sym_power` definitions
      (rpa.cpp/g0w0.cpp/bse.cpp/tdhf.cpp all call linalg::sym_power), zero
      local `dgeev_` extern/wrapper blocks (eom_ea/eom_ip/eom_ee/tdhf/adc all
      route linalg::eig_nonsymmetric), and zero local PCG loops in
      cphf_pack.cpp/ucphf_pack.cpp; sym_power + eig_nonsymmetric are declared in
      include/atomiverse/cpu/linalg/eigensolve.h, solve_pcg in
      include/atomiverse/cpu/linalg/linear_solve.h, and OrbitalHessian::direct_rpa
      in internal/orbital_hessian.h.
    status: verified
    last_checked: 2026-07-24
  - id: ac-002
    summary: New-primitive unit tests build and hit hard-coded goldens
    type: code
    pass_when: |
      test_cpu_linalg_eigensolve and test_cpu_linalg_linear_solve build and pass
      under ctest; sym_power(S,-0.5) reproduces the hard-coded S^{-1/2} to <=1e-12,
      eig_nonsymmetric returns {2,3} for the triangular fixture, solve_pcg matches
      the dense solve to <=1e-10 on the SPD fixture and throws on max_iter exceeded.
    status: verified
    last_checked: 2026-07-24
  - id: ac-003
    summary: Existing response regressions still reproduce cited values
    type: runtime
    pass_when: |
      after the port, the existing scf regressions that route through the hoisted
      code (direct-rpa-h2-sto3g, g0w0-h2-sto3g, ch2o_tdhf, adc2-h2co-sto3g,
      eom-ee-ccsd, eom-ea-ccsd, eom-ip-ccsd, static-polariz-h2o) all build and
      reproduce their cited values unchanged (exit 0).
    status: verified
    last_checked: 2026-07-24
  - id: ac-004
    summary: New linalg-hoist regression reproduces cited dRPA E_c and CPHF alpha
    type: runtime
    pass_when: |
      regressions/scf/scf-review-01-linalg/reg_scf-review-01-linalg builds, runs
      <=60 s, exits 0, and reproduces (a) the cited dRPA correlation energy
      (Furche PRB 64,195120 methodology; same literature value gated by
      direct-rpa-h2-sto3g) to <=1e-12 Ha with kOmega2Tol=1e-10 preserved, and
      (b) the cited CPHF static isotropic polarizability to <=1e-6 a.u.; all
      reference values are hard-coded literals with a citation/date comment and
      no live third-party oracle is invoked.
    status: verified
    last_checked: 2026-07-24
  - id: ac-005
    summary: No new architecture-scrub violations from the hoisted symbols
    type: code
    pass_when: |
      tests/architecture/no_factory_scrub.sh and
      tests/cpu/compatibility_surface_scrub.sh report NO violation attributable
      to the new symbols — sym_power, eig_nonsymmetric, and solve_pcg are
      RealMatrix/RealVector free math with no make_/build_/*_config surface
      introduced. (Both scrubs are red on dev with PRE-EXISTING violations —
      build_imds x3 / adc build_m*_ x3 / build_se_part owned by specs 03 and 05,
      gpw molrs::Element sites owned by 05 — their count must not grow here;
      full green is 03/05's gate, not this spec's.)
    status: verified
    last_checked: 2026-07-24
  - id: ac-006
    summary: New symbols carry Doxygen with units
    type: docs
    pass_when: |
      sym_power and eig_nonsymmetric (cpu/linalg/eigensolve.h), solve_pcg
      (cpu/linalg/linear_solve.h), and OrbitalHessian::direct_rpa each carry
      @brief + @param (with units: Hartree for A/B/M, dimensionless for S) + @return.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-005 (code):** the hoist is structurally complete — one
  definition per primitive in `cpu/linalg`, all in-file copies deleted, new
  unit tests green on hard-coded goldens, and the factory / compatibility scrubs
  confirm the new free functions are the sanctioned no-owner math exception, not
  factories. Verified by `/mol:impl` at delivery.
- **ac-003 (runtime):** behavior-preservation gate — every existing response
  regression that consumes the hoisted code still reproduces its cited value,
  proving the port did not move any number.
- **ac-004 (runtime, the required regression):** the permanent numeric-invariant
  gate for this sub-spec. dRPA E_c is held to the scientist's ≤1e-12 Ha bar
  (bitwise is the target, ≤1e-12 the verifiable floor; resolves open question
  (c)), and the CPHF PCG path is gated by a cited static polarizability. Cited
  literals only — no oracle at ctest time.
- **ac-006 (docs):** public/internal Doxygen with units per repo convention.
