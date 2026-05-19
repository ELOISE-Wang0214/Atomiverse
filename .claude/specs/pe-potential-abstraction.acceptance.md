---
slug: pe-potential-abstraction
criteria:
  - id: ac-001
    summary: Design note committed before any source change lands
    type: docs
    pass_when: |
      `.claude/notes/potential-base-classes-design.md` exists
      and was committed BEFORE the first source-tree change in
      this spec's commit chain. The note covers all six points
      in the "Design pass" task (base-class headers, port
      sequence, ownership helper, `final` strategy, backwards-
      compat aliases, CUDA mirror layout) and explicitly
      addresses each "high"-severity finding from
      `.claude/notes/potential-abstraction-audit.md`.
    status: pending

  - id: ac-002
    summary: Four CPU base classes exist with LAMMPS-style shape
    type: code
    pass_when: |
      `include/atomiverse/cpu/potentials/{pair,bond,angle,dihedral}.h`
      all exist. Each declares an abstract class deriving
      from `cpu::ComponentBase`, with a `final` override of
      `launch()` that calls a pure-virtual `compute()`, a
      `role()` returning `Evaluator`, and a `capabilities()`
      default. `cpu::Pair` additionally exposes `cutoff()`.
      `cpu::Bond` / `cpu::Angle` / `cpu::Dihedral` expose a
      protected `select_owned_tuples(...)` helper for the
      multi-rank ownership filter.
    status: pending

  - id: ac-003
    summary: Four existing CPU potentials rewritten under LAMMPS names
    type: code
    pass_when: |
      All four exist and derive from the matching base:
      `cpu::PairLJCut : public cpu::Pair`,
      `cpu::BondHarmonic : public cpu::Bond`,
      `cpu::AngleHarmonic : public cpu::Angle`,
      `cpu::DihedralMultiHarmonic : public cpu::Dihedral`.
      The matching old `.cpp` files (`src/cpu/md/lj.cpp`,
      `src/cpu/potentials/harmonic_bond.cpp`, etc.) are
      deleted or reduced to a one-line include-only stub.
    status: pending

  - id: ac-004
    summary: Backwards-compat aliases preserve existing user-code API
    type: code
    pass_when: |
      `include/atomiverse/cpu/md/lj.h` contains
      `using LJ = PairLJCut;` (or equivalent type alias).
      Same for `using HarmonicBond = BondHarmonic;`,
      `using HarmonicAngle = AngleHarmonic;`,
      `using MultiHarmonicDihedral = DihedralMultiHarmonic;`.
      Existing user code referring to the old names
      (`sim.set_pair(cpu::md::LJ{...})`, `sim.add(cpu::HarmonicBond{...})`)
      compiles without change.
    status: pending

  - id: ac-005
    summary: Zero-overhead — PairLJCut::compute perf within ±2% of pre-refactor baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      `./build/benchmarks/pe_profile_single --chains-per-axis 4
      --n-atoms-per-chain 18 --steps 1000 --warmup 100 --out
      /tmp/profile_after.csv` on the same host produces a `LJ`
      row whose `total_ms` is within ±2% of the pre-refactor
      154.79 ms recorded in
      `benchmarks/md/reference/polyethylene/profile_1.csv`.
      If outside ±2%, the design note must record the root
      cause (missed `final`, missed inline, copy not elided)
      and the resolution before this criterion flips to
      verified.
    status: pending

  - id: ac-006
    summary: Zero-overhead — BondHarmonic::compute perf within ±2% of pre-refactor baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      Same as ac-005 for the `HarmonicBond` row in
      `profile_1.csv`: pre-refactor 117.36 ms; new value
      within ±2%.
    status: pending

  - id: ac-007
    summary: LJ NVE bit-identicality preserved at 5e-12
    type: scientific
    pass_when: |
      `ctest -R test_lj_nve --output-on-failure` passes. The
      NVE total energy at step 1000 differs from the
      pre-refactor reference by ≤ 5e-12 relative (parent
      `parallel-md-engine.acceptance.md` ac-012 FMA-tolerance
      envelope). The refactor is a name-and-structure change
      only; FP behaviour must not move.
    status: pending

  - id: ac-008
    summary: All existing PE + LJ tests continue green
    type: runtime
    pass_when: |
      `cmake --build build -j && ctest --test-dir build
      --output-on-failure` exits 0. All `test_pe_*`,
      `test_lj_*`, `test_mpi_pe_*`, `test_mpi_lj_*` tests
      pass at every rank count they're labelled for (1, 2,
      8). The four new regression tests (`test_pair_lj_cut`,
      `test_bond_harmonic`, `test_angle_harmonic`,
      `test_dihedral_multi_harmonic`) also pass.
    status: pending

  - id: ac-009
    summary: "One new class" demo — PairMorse reaches NVE with two new files
    type: code
    pass_when: |
      `tests/cpu/md/test_pair_morse_smoke.cpp` exists and
      passes. The Morse implementation introduced to support
      it consists of exactly: `include/atomiverse/cpu/potentials/pair_morse.h`
      and `src/cpu/potentials/pair_morse.cpp`. The design
      note records this file count in a "Live extension demo"
      subsection. No other file in the codebase is modified
      to enable the new style (no central registry, no
      Driver edit — the user just `add(cpu::PairMorse{...})`).
    status: pending

  - id: ac-010
    summary: CUDA mirror compiles clean (no tests)
    type: code
    pass_when: |
      `include/atomiverse/cuda/potentials/{pair,bond,angle,dihedral}.cuh`
      all exist. The four CUDA derived classes (`cuda::PairLJCut`,
      `cuda::BondHarmonic`, `cuda::AngleHarmonic`,
      `cuda::DihedralMultiHarmonic`) exist with the same
      shape as their CPU peers. Backwards-compat aliases
      (`using cuda::LJ = cuda::PairLJCut`, etc.) preserve
      the old names. `cmake --build build --target atomiverse_cuda`
      succeeds with no new warnings. No CUDA tests are added
      by this spec (per user direction); the compile gate
      is the only check.
    status: pending

  - id: ac-011
    summary: Audit findings consumed and cross-referenced
    type: docs
    pass_when: |
      `.claude/notes/potential-base-classes-design.md`
      references `.claude/notes/potential-abstraction-audit.md`
      (produced by `pe-ff-migrate-dreiding-ua`) and explicitly
      addresses each "high"-severity gap in the audit by name.
      Gaps the design does not address are listed at the
      bottom of the design note with an explicit deferral
      rationale.
    status: pending
---

# Acceptance criteria

- **ac-001** locks in "design before code" — the design note is the contract this rewrite implements.
- **ac-002 / ac-003 / ac-004** are structural file-existence + content-regex checks — auto-flippable by `/mol:impl`.
- **ac-005 / ac-006** are the **zero-overhead enforcement gates**. User's stated constraint translates into measurable ±2% bounds against the `profile_1.csv` numbers committed by `parallel-md-finish-01`. Owed to `/mol:bench` or `/mol:close --manual` after manual re-run.
- **ac-007** is the scientific load-bearing gate: refactor must not change FP behaviour. Inherits the 5e-12 FMA-reorder envelope.
- **ac-008** is the runtime no-regression floor.
- **ac-009** is the **live demonstration**: adding a new pair style in exactly two files (`.h` + `.cpp`) — visible proof the architecture works.
- **ac-010** is the CUDA compile-only gate per user direction. Compile clean = pass; no functional tests.
- **ac-011** ensures the audit input is consumed, not orbited.
