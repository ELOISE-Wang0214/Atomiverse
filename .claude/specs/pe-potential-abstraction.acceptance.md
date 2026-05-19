---
slug: pe-potential-abstraction
criteria:
  - id: ac-001
    summary: Design note committed before any source change lands
    type: docs
    pass_when: |
      `.claude/notes/potential-abstraction-design.md` exists
      and was committed BEFORE the first source-tree change in
      this spec's commit chain. The note covers all six points
      listed in the "Design pass" task (Form contract, Walker
      skeletons, ATV_HOST_DEVICE macro, SIMD opt-in fast path,
      multi-rank halo survival, name-resolution layout).
    status: pending

  - id: ac-002
    summary: Form headers exist and are header-only
    type: code
    pass_when: |
      `include/atomiverse/potentials/forms/{lj,harmonic_bond,harmonic_angle,multi_harmonic_dihedral}.h`
      all exist. Each declares exactly one `struct` with a
      `Params` member, a `static constexpr inline` (or
      `ATV_HOST_DEVICE`-annotated) `eval` function, and no
      `.cpp` companion. `grep -L "static.*constexpr.*inline"`
      across the four headers returns zero files.
    status: pending

  - id: ac-003
    summary: CPU Walker templates exist with the documented signatures
    type: code
    pass_when: |
      `include/atomiverse/cpu/potentials/{pair,bond,angle,dihedral}_walker.h`
      all exist, each declares a `template <class Form> class
      XWalker : public ComponentBase` with `build` / `launch` /
      `destroy` overrides, `role()` returns `Evaluator`,
      `capabilities()` returns at least `Energy | Force`. A
      `static_assert` inside the template verifies
      `Form::eval` is invocable with the expected signature.
    status: pending

  - id: ac-004
    summary: Existing CPU potentials collapsed to type aliases
    type: code
    pass_when: |
      `include/atomiverse/cpu/md/lj.h` declares
      `using LJ = PairWalker<potentials::LJForm>` (or
      equivalent) — no manual `class LJ` body remains.
      Same for `cpu::HarmonicBond`, `cpu::HarmonicAngle`,
      `cpu::MultiHarmonicDihedral`. `src/cpu/md/lj.cpp` is
      either deleted or reduced to a one-line stub. Same for
      the three bonded `.cpp` files.
    status: pending

  - id: ac-005
    summary: Zero-overhead — LJ::launch perf within ±2% of pre-refactor baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      `pe_profile_single --chains-per-axis 4 --n-atoms-per-chain 18 --steps 1000 --warmup 100`
      under the new code reports `LJ` row `total_ms` within
      ±2% of the pre-refactor 154.79 ms recorded in
      `benchmarks/md/reference/polyethylene/profile_1.csv` on
      the same host (Apple M4, Homebrew Open-MPI, OMP=10).
      If outside ±2%, the design note must record the
      inlining-failure root cause and the assembly-inspection
      evidence; mere regression is a failure.
    status: pending

  - id: ac-006
    summary: HarmonicBond::launch perf within ±2% of pre-refactor baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      Same as ac-005 but for the `HarmonicBond` row in
      `profile_1.csv`: pre-refactor 117.36 ms; new value within
      ±2%.
    status: pending

  - id: ac-007
    summary: Bit-identical / 5e-12 LJ NVE preserved across refactor
    type: scientific
    pass_when: |
      `ctest -R test_lj_nve --output-on-failure` passes; the
      total NVE energy at step 1000 differs from the
      pre-refactor reference by ≤ 5e-12 relative (the same
      tolerance parent `parallel-md-engine.acceptance.md` ac-012
      sets for the FMA-reorder envelope). The refactor moves
      code around but cannot change FP behaviour.
    status: pending

  - id: ac-008
    summary: All existing PE tests continue green at 1 and 8 ranks
    type: runtime
    pass_when: |
      `ctest --test-dir build --output-on-failure` exits 0.
      All `test_pe_*` (single-chain, multi-chain, NVE, NVT),
      `test_mpi_pe_*` (chains_multirank, chains_8rank,
      chain_crossing, chain_drifting, full_multirank,
      nvt_multirank, msu_pipeline_n8), and `test_lj_*`
      tests pass. Total test count ≥ (pre-Phase-5 baseline 87
      + new tests added by this spec). No PE / LJ test is
      skipped or removed in the refactor.
    status: pending

  - id: ac-009
    summary: "One new file" demonstration — Morse bond reaches NVE in ≤ 2 touched files
    type: code
    pass_when: |
      `tests/cpu/md/test_potential_walker_morse_smoke.cpp`
      exists and passes. The Morse implementation introduced
      to support it consists of exactly: (a) one new header
      `include/atomiverse/potentials/forms/morse.h` declaring
      `MorseForm` with `Params { D, alpha, r0 }` + `eval`,
      and (b) one `using MorseBond = BondWalker<MorseForm>`
      type alias (the test header is permitted to host this
      alias). `git log --name-only --diff-filter=A
      --pretty=format: HEAD~1..HEAD -- include/atomiverse/potentials/`
      shows ≤ 1 Form header added per "new FF form" commit.
      The design note records the touched-file count in a
      "Live extension demo" subsection.
    status: pending

  - id: ac-010
    summary: CUDA mirror compiles clean (no tests required)
    type: code
    pass_when: |
      `include/atomiverse/cuda/potentials/{pair,bond,angle,dihedral}_walker.cuh`
      all exist. CUDA classes (`cuda::HarmonicBond`,
      `cuda::HarmonicAngle`, `cuda::MultiHarmonicDihedral`,
      `cuda::LJ`) reduced to type aliases like the CPU side.
      `cmake --build build --target atomiverse_cuda` compiles
      without error or new warnings. No CUDA tests are
      added by this spec (per user direction). The compile
      gate is the only check.
    status: pending

  - id: ac-011
    summary: Audit findings consumed and cross-referenced
    type: docs
    pass_when: |
      The design note `.claude/notes/potential-abstraction-design.md`
      references `.claude/notes/potential-abstraction-audit.md`
      (produced by `pe-ff-migrate-dreiding-ua`) and explicitly
      addresses each "high"-severity gap in the audit by name.
      Gaps the design does not address are listed at the
      bottom of the design note with an explicit deferral
      rationale.
    status: pending
---

# Acceptance criteria

- **ac-001** locks in "design before code." The design note is the contract this refactor implements; without it the spec degenerates into ad-hoc rewriting.
- **ac-002 / ac-003 / ac-004** are structural file-existence + content-regex checks — auto-flippable by `/mol:impl`.
- **ac-005 / ac-006** are the **zero-overhead enforcement gates**. The user's stated constraint ("called many times, zero overhead") translates into measurable ±2% bounds against the pre-refactor `profile_1.csv` numbers committed by `parallel-md-finish-01`. These criteria are owed to `/mol:bench` (or a manual re-run on the dev box with `/mol:close --manual` if the bench harness isn't configured for this benchmark).
- **ac-007** is the scientific load-bearing gate: the refactor must not change FP behaviour. Bit-identicality at the 5e-12 envelope inherits from the parent FMA-reorder tolerance refinement.
- **ac-008** is the runtime no-regression floor for the full test suite.
- **ac-009** is the **live demonstration**: implementing Morse-bond in ≤ 2 new files is the visible proof the abstraction works as advertised. If this fails, the abstraction may be technically present but practically not useful.
- **ac-010** is the CUDA compile-only gate per user direction. No tests; just the compile.
- **ac-011** ensures this spec actually consumes the audit input rather than orbiting it.
