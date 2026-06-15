---
slug: minimizer-line-search
criteria:
  - id: ac-001
    summary: Strong-Wolfe search returns alpha satisfying N&W 3.7a and 3.7b on a quadratic
    type: scientific
    pass_when: |
      test_cpu_minimizer_line_search asserts, on the convex quadratic
      fixture, that the alpha returned by strong_wolfe() satisfies both
      f(x+αp) ≤ f(x)+c1·α·gᵀp (3.7a) and |g(x+αp)ᵀp| ≤ c2·|gᵀp| (3.7b)
      numerically, and the test passes.
    status: pending
  - id: ac-002
    summary: Armijo backtracking mode returns alpha satisfying N&W 3.6a
    type: scientific
    pass_when: |
      test_cpu_minimizer_line_search asserts the Armijo-mode alpha
      satisfies f(x+αp) ≤ f(x)+c1·α·gᵀp on the quadratic fixture and the
      test passes.
    status: pending
  - id: ac-003
    summary: Descent precheck resets p=-g when an uphill direction is injected
    type: code
    pass_when: |
      A test injects a direction with gᵀp ≥ 0; the optimizer/helper
      replaces it with p=-g (gᵀp < 0) before searching, asserted by
      test_cpu_minimizer_line_search, and the test passes.
    status: verified
    last_checked: 2026-06-15
  - id: ac-004
    summary: Every stored L-BFGS (s,y) pair satisfies yᵀs > 0
    type: scientific
    pass_when: |
      test_cpu_minimizer_lbfgs asserts yᵀs > eps·||s||·||y|| for every
      curvature pair stored across a full minimization of the quadratic;
      the test passes.
    status: pending
  - id: ac-005
    summary: reeval_forces handle re-runs the force phase at the current geometry
    type: code
    pass_when: |
      cpu::Context exposes std::function<void(bool)> reeval_forces;
      OptimizeDriver::setup_pass binds a non-null closure that, when
      invoked, updates frame.energy (and frame.force when need_grad) for
      the frame's current geometry, exercised by a line-search unit test.
    status: verified
    last_checked: 2026-06-15
  - id: ac-006
    summary: LBFGS reaches max-force < 0.01 eV/Å with strong-Wolfe enabled
    type: scientific
    pass_when: |
      A regression test optimizes the molecular fixture (H2O or aspirin
      at RHF/STO-3G) with default strong-Wolfe LBFGS and reaches
      max atomic force < 0.01 eV/Å within the iteration budget; the
      test passes.
    status: pending
  - id: ac-007
    summary: Strong-Wolfe path uses fewer force evaluations than fixed-step to the same minimum
    type: scientific
    pass_when: |
      The regression test records force-evaluation count and final
      energy for both the fixed-step and strong-Wolfe paths on the same
      fixture; strong-Wolfe reaches the tight threshold in fewer force
      evaluations (or fewer outer steps) and the final energies agree
      within tolerance; the test passes.
    status: pending
  - id: ac-008
    summary: CG strong-Wolfe (c2=0.1) with PR+ converges on the quadratic and a small molecule
    type: scientific
    pass_when: |
      test_cpu_minimizer_cg drives CG with default strong-Wolfe
      (c2=0.1) and PR+ direction; it converges to the analytic minimum
      of the quadratic and to max-force < 0.01 eV/Å on the small
      molecule; the test passes.
    status: pending
  - id: ac-009
    summary: Fixed-step fallback reproduces legacy behaviour
    type: code
    pass_when: |
      With LineSearchMode::Fixed selected, the optimizer reproduces the
      pre-change fixed-step trajectory (regression assertion in
      test_cpu_minimizer_lbfgs) and the test passes.
    status: verified
    last_checked: 2026-06-15
  - id: ac-010
    summary: Line search makes identical decisions and final energy across MPI ranks
    type: scientific
    pass_when: |
      test_mpi_minimize_ls_n1, _n2, _n4 each pass: the recorded
      accept/reject decision sequence and the final converged energy are
      identical (within tight tolerance) across the three rank counts.
    status: pending
  - id: ac-011
    summary: Example 05 converges at the tight threshold without the loosened max_dtol
    type: runtime
    pass_when: |
      examples/05_geom_opt.cpp no longer sets a loosened max_dtol, builds,
      and runs to convergence at the tight force threshold; the full
      check + test suite (ctest) passes.
    status: verified
    last_checked: 2026-06-15
  - id: ac-012
    summary: All new/changed public symbols carry Doxygen with units
    type: docs
    pass_when: |
      LineSearchMode, strong_wolfe/armijo_backtrack, the LBFGS/CG
      line_search setters, and Context::reeval_forces each have @brief and
      (where applicable) @param/@return with units per the project Doxygen
      convention.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 (scientific):** the binding numerical check that the search satisfies the cited N&W conditions — the mathematical heart of the spec.
- **ac-003 (code):** the descent-direction safeguard is observable in a unit test.
- **ac-004 (scientific):** positive curvature guarantee that keeps the BFGS update PD.
- **ac-005 (code):** the re-entrant force-eval mechanism exists and works — the design's structural keystone.
- **ac-006 / ac-007 / ac-008 (scientific):** intrinsic convergence/efficiency validation on real systems for both optimizers.
- **ac-009 (code):** the Fixed fallback is preserved.
- **ac-010 (scientific):** DD parity across n=1,2,4.
- **ac-011 (runtime):** the user-visible payoff — example converges tightly without loosened tolerances.
- **ac-012 (docs):** Doxygen coverage per project convention.
