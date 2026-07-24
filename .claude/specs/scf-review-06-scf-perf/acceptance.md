---
slug: scf-review-06-scf-perf
criteria:
  - id: ac-001
    summary: XC-basis cache path is bitwise-equal to the per-iteration rebuild
    type: code
    pass_when: |
      tests/cpu/scf/test_xc_grid_engine.cpp (TestXcGridEngine) builds Phi/dPhi
      via XCGridBasis::build then calls the precomputed-basis
      evaluate_xc_batched overload, and asserts Vxc_alpha, Vxc_beta, E_xc,
      n_alpha_integrated, n_beta_integrated are bit-for-bit equal to the
      original rebuild path for both an LDA and a PBE (GGA) case on the stored
      H2O grid fixture; a grep assertion confirms the dead `active[]` array is
      gone. Test compiles under -DATV_BUILD_HF=ON and passes.
    status: verified
    last_checked: 2026-07-24
  - id: ac-002
    summary: HamiltonianBuilder rebuilds the AO-on-grid cache only per geometry
    type: code
    pass_when: |
      With the mutable xc_basis_ cache, prepare() sets xc_basis_.built = false
      and step() (const) lazily builds on the first KS iteration keyed off the
      cached grid_ and is_gga(xc_type); two consecutive step() calls at a fixed
      geometry produce identical Vxc with the AO-on-grid eval executed exactly
      once (verified in TestXcGridEngine via a build-count probe or by asserting
      the second step performs no eval_basis pass). Existing RKS/UKS solver unit
      tests (tests/cpu/scf/test_uks_solver.cpp, test_scf_solver.cpp) still pass.
    status: verified
    last_checked: 2026-07-24
  - id: ac-003
    summary: Fused jk_rhf is algebraically exact vs the two-pass g_rhf + k
    type: code
    evaluator_hint: "gate: max_abs(F_fused - F_twopass) <= 1e-12 Ha"
    pass_when: |
      tests/cpu/scf/test_operators.cpp (TestShellBlockPairList) builds J and K
      from ShellBlockPairList::jk_rhf in one sweep and, separately, J via g_rhf
      (J = G + 0.5*K) and K via k() on a stored H2O/def2-SVP shell-block fixture
      containing d shells; for alpha_hf in {0.20 (B3LYP), 0.25 (PBE0)} it
      assembles F_fused = Hcore + J - 0.5*alpha_hf*K + Vxc and F_twopass and
      asserts per-element max_abs(F_fused - F_twopass) <= 1e-12 Ha. K is emitted
      unscaled (RSH-separable). Test passes.
    status: verified
    last_checked: 2026-07-24
  - id: ac-004
    summary: operators.cpp micro-opts preserve numerics and are thread-invariant
    type: code
    pass_when: |
      TestIntegralPack asserts IntegralPack::coulomb_exchange_fused produces
      bit-for-bit identical G for omp_set_num_threads(1) vs N after the K
      contraction gains an OpenMP reduction; and a grep/inspection confirms
      ShellBlockPairList::k and g_rhf no longer zero-init blk[] with `= {}` but
      zero exactly na*nb*nc*nd. Existing ERI unit tests (test_rys_eri.cpp,
      test_integrals.cpp) and TestShellBlockPairList remain green.
    status: verified
    last_checked: 2026-07-24
  - id: ac-005
    summary: thread_local ERI reuse removes per-quartet heap alloc, keeps values
    type: code
    pass_when: |
      tests/cpu/scf/test_rys_eri.cpp asserts eri_block values for a fixed
      (la,lb,lc,ld) block quartet are identical single-threaded vs inside a
      #pragma omp parallel region after ang_idx (rys_eri_block_canonical) and
      the non-canonical scratch (rys_eri_block) become static thread_local
      (resize-if-smaller, scratch prefix filled before accumulation); a grep
      assertion confirms no local `std::vector<AngIdx>` remains in the canonical
      kernel. Test passes.
    status: verified
    last_checked: 2026-07-24
  - id: ac-006
    summary: Fused-path B3LYP H2O regression reproduces the cited PySCF energy
    type: runtime
    evaluator_hint: "regression target: b3lyp-fused-jk-h2o"
    pass_when: |
      regressions/scf/b3lyp-fused-jk-h2o/reg_b3lyp_fused_jk_h2o.cpp runs a
      restricted B3LYP/def2-SVP H2O single point through the fused jk_rhf path
      and reproduces the hard-coded PySCF total energy in
      regressions/scf/b3lyp-fused-jk-h2o/REFERENCE.md (tool, version, command,
      date embedded as literals; no third-party call at runtime) within the
      stated tolerance (<= 1e-6 Ha), in <= 60 s. Built behind
      ATV_BUILD_REGRESSIONS and run via the regressions target.
    status: verified
    last_checked: 2026-07-24
  - id: ac-007
    summary: All DFT/SCF regression families still reproduce their cited values
    type: runtime
    pass_when: |
      The full unit ctest passes and the DFT/SCF regression families reproduce
      their cited values unchanged (<= 1e-10 Ha drift on hybrid SCF energies):
      regressions/scf/{h2_full_grid_alda, blyp-water-sto3g, full-lr-tddft-h2-pbe,
      uks-full-lr-o2-lsda, uks-tda-o2-lsda, h2-b3lyp-hybrid-tda,
      h2-b3lyp-rpa-tddft, b2plyp-h2o, d3bj-molecular-pbe, dft-d4-h2o-pbe,
      rks-hessian-h2o-sto3g, uks-hessian-o2-sto3g, sf-tddft-o2, spin-flip-tda-o2}.
    status: verified
    last_checked: 2026-07-24
  - id: ac-008
    summary: Advisory before/after SCF Fock wall-time on one regression case
    type: performance
    evaluator_hint: "advisory: reported in impl summary, agent-auto; not a ctest gate"
    pass_when: |
      The implementation summary reports a measured before/after wall-time for
      the SCF Fock build (or full single point) on ONE regression case
      (e.g. b3lyp-fused-jk-h2o or h2-b3lyp-hybrid-tda) showing the fused-JK +
      Phi-cache path is not slower and is expected faster; this is advisory
      evidence only — no numeric threshold gates the build or ctest.
    status: verified
    last_checked: 2026-07-24
    verified_by: agent-auto
    note: measured (git-stash A/B, 2 runs each) B3LYP/def2-SVP H2O single
      point — before 43.86/45.13 s, after 39.19/36.90 s (~14% faster on a
      24-basis fixture; the Phi-cache + one-sweep win scales with system size).
---

# Acceptance criteria

- **ac-001 / ac-002 (Item 15)** — Phi/dPhi cache correctness and per-geometry lifecycle. The cache must be a pure performance change: cached vs rebuild is bitwise-equal, and the AO-on-grid evaluation runs once per geometry, not once per SCF iteration.
- **ac-003 (Item 16)** — the binding scientist invariant. Fusing J and K into one sweep is algebraically exact; the ≤1e-12 per-element gate exists specifically to catch the classic one-sweep-JK bug where the 8-fold ERI permutation multiplicities unpack differently for J than for K. K is emitted unscaled so a later RSH erf-attenuated K can be added.
- **ac-004 / ac-005 (Item 20)** — the three micro-hotspots (block zero-init, serial stored-ERI K, per-quartet heap alloc) preserve numerics and stay thread-count-invariant.
- **ac-006 / ac-007** — end-to-end numeric guard: the new fused restricted-hybrid path reproduces a cited PySCF value, and every existing DFT/SCF regression family is unchanged.
- **ac-008** — advisory perf evidence only (benchmarks/ is removed; no timing in tests/): a measured before/after wall-time on one regression case, reported in the impl summary, never gating ctest.
