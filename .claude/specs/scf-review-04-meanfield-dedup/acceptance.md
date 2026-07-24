---
slug: scf-review-04-meanfield-dedup
criteria:
  - id: ac-001
    summary: FD Hessian hoisted to a protected CRTP template; RKS/UKS/ROHF delegate
    type: code
    evaluator_hint: "grep: central_difference_nuclear_hessian"
    pass_when: |
      MeanFieldMethod<Derived> declares protected central_difference_nuclear_hessian
      (mean_field_method.h) defined in mean_field_method.cpp; rks.cpp/uks.cpp/rohf.cpp
      hessian() bodies contain no ±h FD displacement loop and call the template;
      the tree builds with -DATV_BUILD_HF=ON.
    status: pending
  - id: ac-002
    summary: FD-template unit test reproduces a quadratic mock's analytic Hessian
    type: code
    evaluator_hint: "ctest -R test_fd_hessian_template"
    pass_when: |
      tests/cpu/scf/test_fd_hessian_template.cpp (analytic F=-K(R-R0) mock, no SCF)
      passes: recovered Hessian equals the hard-coded K in eV/Å² to 1e-8, geometry
      restored to r0, ordering col=3A+d verified.
    status: pending
  - id: ac-003
    summary: Hessian precondition guard deduped to a private CRTP helper
    type: code
    evaluator_hint: "grep: check_hessian_preconditions"
    pass_when: |
      MeanFieldMethod<Derived> has private check_hessian_preconditions(method_name, pack)
      throwing on non-convergence, spherical pack, and ECP pack; rhf.cpp and uhf.cpp
      hessian() call it after run(); the RI guard stays inline pre-run; builds green.
    status: pending
  - id: ac-004
    summary: GradientResult gains FromHartreeBohr ctor; make_* glue deleted
    type: code
    evaluator_hint: "grep: make_result|make_force_result"
    pass_when: |
      gradient_result.h defines FromHartreeBohr{} tag, converting ctor, default ctor,
      and member-wise eV ctor; make_result and make_force_result appear nowhere under
      src/cpu/scf; gradient.cpp / mp2_gradient.cpp / ccsd_gradient.cpp build against
      the tag ctor.
    status: pending
  - id: ac-005
    summary: Converting-ctor unit test matches hard-coded eV/Å + eV literals
    type: code
    evaluator_hint: "ctest -R test_gradient_result"
    pass_when: |
      tests/cpu/scf/test_gradient_result.cpp passes: GradientResult(FromHartreeBohr{},
      grad_ha_bohr, energy_ha) yields forces=-convert(..HartreePerBohr->eVPerAngstrom)
      and energy=convert(Hartree->eV) to hard-coded literals; default and brace
      member-wise construction still compile.
    status: pending
  - id: ac-006
    summary: Existing FD-Hessian and gradient regressions reproduce byte-identical values
    type: runtime
    evaluator_hint: "regressions target; byte-identical replay vs pinned PySCF 2.14.0"
    pass_when: |
      regressions/scf/{rks-hessian-h2o-sto3g,uks-hessian-o2-sto3g,rohf-hessian-o2,
      uhf-hessian-o2,mp2-grad-h2o-sto3g,ump2-grad-o2-sto3g,ccsd-grad-h2o-sto3g} all
      pass their existing gates AND emit numbers bit-identical to pre-refactor
      (e.g. RKS H_xx(0,0)=-9.632120712096 eV/Å²; MP2 O ∂E/∂z=-9.528798092138e-02
      Ha/Bohr; CCSD O ∂E/∂z=-1.106969560091e-01 Ha/Bohr).
    status: pending
  - id: ac-007
    summary: FD replay invariants (h, fresh-solver, ordering, threading) frozen
    type: code
    evaluator_hint: "grep: scf_dedup_scrub.sh"
    pass_when: |
      central_difference_nuclear_hessian freezes h=1e-4 Bohr, constructs a fresh
      Derived per ±h displacement via the caller functor, loops atom-major/axis-minor
      (col=3A+d) serially with no added displacement-level parallelism; the scrub
      tests/architecture/scf_dedup_scrub.sh (no FD loop in rks/uks/rohf.cpp; no
      make_result/make_force_result under src/cpu/scf) is registered and green.
    status: pending
  - id: ac-008
    summary: Doxygen with units on new template method, guard helper, and ctors
    type: docs
    pass_when: |
      central_difference_nuclear_hessian, check_hessian_preconditions, and the new
      GradientResult ctors each carry @brief/@param/@return with units per doc.style
      doxygen.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-003 / ac-004 / ac-007 (code):** the three dedups landed with the
  agreed shape — a protected CRTP FD template, a private precondition guard helper,
  and a tag-disambiguated `GradientResult` converting ctor — and the grep gate
  proves the copy-pasted bodies and the `make_*` factory helpers are gone.
- **ac-002 / ac-005 (code):** RED-first unit tests pin the extracted behavior in
  isolation (no SCF): the FD template against a quadratic mock with a known constant
  Hessian, and the converting ctor against hand-computed eV/Å + eV literals.
- **ac-006 (runtime):** the binding replay gate. Behavior-preserving means the seven
  pre-existing FD-Hessian and MP2/CCSD-gradient regressions reproduce their pinned
  PySCF 2.14.0 literals **bit-identically**, not merely within tolerance. No new
  `regressions/<slug>` is added by design (each touched feature already ships one).
- **ac-008 (docs):** public/protected new symbols documented with units.
