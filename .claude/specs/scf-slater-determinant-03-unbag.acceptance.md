---
slug: scf-slater-determinant-03-unbag
criteria:
  - id: ac-001
    summary: scf_result.h deleted; SCFResult survives only as scf::toy
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/scf_result.h does not exist; `rg -l 'SCFResult'
      include src tests examples regressions` returns only
      include/atomiverse/cpu/scf/solve.h and src/cpu/scf/solve.cpp (toy namespace).
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: Method surface unbagged — determinant()/energy()/iterations(), void run()
    type: code
    pass_when: |
      mean_field_method.h declares `const SlaterDeterminant& determinant() const`,
      `double energy() const` (Hartree), `int iterations() const`; both run()
      overloads return void; MeanFieldOutputs::result is
      shared_ptr<const SlaterDeterminant>; `rg 'last_result|on_evaluate'` over
      include/atomiverse/cpu/scf, src/cpu/scf, src/cpu/recorder, tests/cpu/scf,
      tests/cpu/recorder, regressions/scf and the eight migrated examples
      returns zero hits.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: Pinned parity energies and eigenvalues unchanged through new surface
    type: runtime
    pass_when: |
      test_method_driver_parity passes with the five pinned Hartree literals
      (lines 154-158) and pinned eigenvalue vectors unmodified, reading
      energy()/determinant()/converged().
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: MeanFieldOutputs handle survives method move into a driver
    type: runtime
    pass_when: |
      test_density_cube passes: handle taken before add(std::move(scf)) still
      yields a valid determinant view after the run (n_writes > 0, integrated
      density matches n_electrons at existing tolerance).
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: Nested force-slot (AIMD) path behavior preserved
    type: runtime
    pass_when: |
      test_method_driver_wire and test_driver_nesting pass: frame.meta kEnergy
      [eV] matches energy()*kHartreeInEv at existing tolerances; inner SCF
      convergence does not stop the outer loop.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: Spin-contamination invariants preserved on SlaterDeterminant
    type: scientific
    pass_when: |
      UHF H/STO-3G and UKS He/STO-3G report |spin_contamination| < kPhysTol via
      determinant(); test_delta_scf doublet invariants (⟨S²⟩ = 0.75 exact cases,
      sub-Koopmans bound) pass at existing tolerances.
    status: pending
  - id: ac-007
    summary: H2O RHF/STO-3G literature regression via Frame meta
    type: scientific
    pass_when: |
      reg_h2o_rhf_sto3g prints PASS: energy read from frame.meta.get(kEnergy)
      converted to Hartree within 1e-3 Ha of -74.961754 Ha; convergence read
      from rhf.converged(); no SCFResult in the TU.
    status: pending
  - id: ac-008
    summary: MP2/FCI route through determinant + e_scf; guards relocated
    type: runtime
    pass_when: |
      test_mp2, test_fci and parity Mp2AndFciStillRouteThroughDrivenMethod pass;
      calling mp2()/fci() on a non-converged reference throws
      std::invalid_argument (asserted by unit test).
    status: verified
    last_checked: 2026-07-20
  - id: ac-009
    summary: Internal packs consume const SlaterDeterminant&
    type: code
    pass_when: |
      gradient.h, integral_derivatives.h, mp2_pack.h, fci_pack.h contain no
      SCFResult token and no scf_result.h include; ScfGradient::compute and
      accumulate_two_electron_gradient take const SlaterDeterminant&; build
      succeeds with ATV_BUILD_HF=ON.
    status: verified
    last_checked: 2026-07-20
  - id: ac-010
    summary: OrbitalInfoObserver uses det.homo()/lumo(), values unchanged
    type: runtime
    pass_when: |
      test_orbital_info_observer passes with identical HOMO/LUMO/gap eV values;
      orbital_info_observer.cpp no longer indexes eigenvalues_alpha manually
      (scan deleted, homo()/lumo() called).
    status: verified
    last_checked: 2026-07-20
  - id: ac-011
    summary: Full suite and examples green
    type: runtime
    pass_when: |
      cmake build with ATV_BUILD_TESTS=ON ATV_BUILD_HF=ON ATV_BUILD_EXAMPLES=ON
      compiles all migrated examples; ctest full suite passes with no new
      failures vs the pre-change baseline.
    status: verified
    last_checked: 2026-07-20
  - id: ac-012
    summary: Docs off the SCFResult user path
    type: docs
    pass_when: |
      tutorial ch03/ch04/ch06/ch07, developer/cpu-internals.md, api/index.md,
      roadmap/post-hf.md no longer present SCFResult/last_result() as the user
      path; methods.h and single_point.h doc comments show
      determinant()/energy()/frame.meta.
    status: pending
---

# Acceptance criteria

ac-001/002/009 are the debt-termination gates (bag gone, surface unbagged, internals typed). ac-003–005 pin behavior preservation of the three run entries (standalone, Evaluate component, nested force-slot). ac-006/007 are the domain-validation gates (spin invariants per Amos & Hall 1961; cited literature energy). ac-008/010 cover the two consumer relocations with semantic decisions (guard move, homo/lumo delegation). ac-011/012 close the loop on the whole tree and the documented user path.
