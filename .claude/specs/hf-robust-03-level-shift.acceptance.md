---
slug: hf-robust-03-level-shift
criteria:
  - id: ac-001
    summary: Default level_shift=0.0 leaves all existing reference energies unchanged
    type: code
    pass_when: |
      SCFConfig::level_shift defaults to 0.0; the full existing test suite
      (test_cpu_scf_solver, test_cpu_scf_accelerators, UHF/UKS tests) passes
      with zero edits to existing expected values, and the new anchor assert
      E(H2O/6-31G RHF, cart) = -75.9839744727 Ha holds to 1e-6 Ha at λ=0.
    status: verified
    last_checked: 2026-07-02
  - id: ac-002
    summary: Unit test proves spectral effect — virtuals shift by λ, occupied unchanged
    type: code
    pass_when: |
      tests/cpu/scf/test_accelerators.cpp contains passing tests where, for a
      synthetic (F, S, D_occ) fixture, generalized eigenvalues of
      apply_level_shift(F,S,D_occ,λ) equal those of F with every virtual
      eigenvalue shifted by exactly λ (tol 1e-12) and every occupied one
      unchanged (tol 1e-12); and apply_level_shift is an element-wise identity
      when |λ| <= 1e-4.
    status: verified
    last_checked: 2026-07-02
  - id: ac-003
    summary: H2O/6-31G RHF energy invariant under λ=0.5 to 1e-9 Ha
    type: scientific
    pass_when: |
      test_cpu_scf_level_shift runs H2O/6-31G RHF (geometry O 0 0 0.1173;
      H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692 Å) with level_shift=0.5 and 0.0,
      both converge, and |E(0.5) - E(0.0)| < 1e-9 Ha.
    status: pending
  - id: ac-004
    summary: Reported orbital energies are shift-free after convergence
    type: scientific
    pass_when: |
      In the same H2O/6-31G pair of runs, max_i |ε_i(λ=0.5) - ε_i(λ=0)| over
      ALL orbitals (occupied and virtual) in result.eigenvalues_alpha is
      <= 1e-6 Ha, demonstrating the post-convergence λ=0 diagonalization.
    status: pending
  - id: ac-005
    summary: Hard case — level shift rescues a divergent stretched-H2O SCF
    type: scientific
    pass_when: |
      test_cpu_scf_level_shift encodes the tutorial hard case (stretched H2O,
      both O-H ≈ 1.7 Å, 6-31G RHF, init_guess=Hcore, use_diis=false,
      damping_factor=0, fixed max_iter): the λ=0 run reports converged==false
      while the λ=1.0 run reports converged==true with pinned iteration count,
      and its energy matches a DIIS-converged reference run of the same
      geometry to < 1e-8 Ha.
    status: pending
  - id: ac-006
    summary: DIIS extrapolates the unshifted Fock; shift applied after, per PySCF
    type: code
    pass_when: |
      In src/cpu/scf/scf_solver.cpp the DIIS error matrices are computed and
      diis_->step invoked strictly before apply_level_shift, which is invoked
      strictly before eigen_->step, with a structural comment documenting
      damping → DIIS → level shift → diagonalize and the density-damping vs
      PySCF-Fock-damping difference; grep finds no code path applying the
      shift before or inside DIIS.
    status: verified
    last_checked: 2026-07-02
  - id: ac-007
    summary: Tutorial section exists and every cited number is regression-pinned
    type: docs
    pass_when: |
      docs/tutorial/ch02-rhf.md contains an English section "SCF convergence
      is a numerical problem" covering DIIS, damping (noting the engine's
      density-damping convention), level shift (Saunders & Hillier, DOI
      10.1002/qua.560070407, post-convergence de-shift rule), and
      Jensen-style practical guidance; every energy, Δε, and iteration count
      quoted in the section appears as an assertion (same constants) in
      tests/cpu/scf/test_level_shift.cpp.
    status: pending
  - id: ac-008
    summary: Default ctest gate respected
    type: runtime
    pass_when: |
      ctest --test-dir build passes in full; tests added to the DEFAULT gate
      by this spec contribute < 2 s combined wall time (any solve-to-convergence
      check belongs in benchmarks/<domain>/validation, not in ctest).
    status: verified
    last_checked: 2026-07-02
---

# Acceptance criteria

- **ac-001 / ac-002 (code)** — no-regression at λ=0 plus the algebraic contract of the helper: the level shift is a pure spectral operation (virtuals +λ, occupied fixed) and a strict no-op below the PySCF skip threshold.
- **ac-003 / ac-004 (scientific)** — PySCF-2.13.1-verified parity target: converged energy invariant under λ (PySCF measured 1.1e-13 Ha; spec bar 1e-9), and the mandatory post-convergence λ=0 diagonalization keeps reported ε shift-free (PySCF measured 2.5e-7; spec bar 1e-6).
- **ac-005 (scientific)** — the tutorial's demonstration case is itself a binding regression: a deterministic setting where λ>0 visibly rescues divergence. Implementation may tune the exact stretch factor/λ during the number-generation run, but test and tutorial must share the same constants and the rescue property is non-negotiable.
- **ac-006 (code)** — the PySCF ordering (damping → DIIS on unshifted F → shift → diagonalize) is the parity-critical invariant; it is enforced structurally in the loop and documented, including the deliberate density-vs-Fock damping divergence kept out of scope.
- **ac-007 (docs)** — tutorial content in English per request; numbers produced by the new code, not hand-copied.
- **ac-008 (runtime)** — the ~16 s default ctest budget from `tests/` policy is preserved.
- **ac-008 measurement (2026-07-02)**: default-gate additions `test_cpu_scf_level_shift` 0.36 s + LevelShift unit tests inside `test_cpu_scf_accelerators` (0.05 s total binary) = < 0.5 s combined, well under the 2 s bar; no E2E-tail entry needed.
