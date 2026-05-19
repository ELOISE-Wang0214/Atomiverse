---
slug: parallel-md-engine
criteria:
  - id: ac-001
    summary: NoseHoover exposes set_T_target runtime setter
    type: code
    pass_when: |
      include/atomiverse/cpu/md/nose_hoover.h declares
      `void set_T_target(double T)` as a public member of
      `cpu::md::NoseHoover`, and src/cpu/md/nose_hoover.cpp
      defines its body.

  - id: ac-002
    summary: Nose-Hoover ramp test exercises set_T_target without driver rebuild
    type: code
    pass_when: |
      tests/mpi/test_nose_hoover_t_ramp.cpp exists; calls
      sim.run() ≥ 5 times interleaved with set_T_target;
      asserts that &sim.backend().analysis_executor() (or
      equivalent backend identity) is unchanged across the
      sequence; ctest -R test_nose_hoover_t_ramp passes at
      n=2 ranks.

  - id: ac-003
    summary: polyethylene_tension.cpp composes the 5-stage protocol in one process
    type: code
    pass_when: |
      benchmarks/md/polyethylene_tension.cpp constructs
      exactly one cpu::md::MD instance and calls
      NoseHoover::set_T_target at least twice in distinct
      stage transitions (verifiable by grep for
      `set_T_target` returning ≥ 2 occurrences in distinct
      stage blocks).

  - id: ac-004
    summary: Multi-stage tension benchmark runs end-to-end at 8 ranks
    type: runtime
    pass_when: |
      `mpirun -n 8 ./build/benchmarks/md/polyethylene_tension`
      completes all 5 stages, exits 0, writes a non-empty
      thermo CSV and stress-strain CSV, with no MPI error
      and resident memory growth ≤ 5% over the run.

  - id: ac-005
    summary: LAMMPS reference inputs present and consistent
    type: code
    pass_when: |
      benchmarks/md/reference/polyethylene/lammps.in,
      thermo_reference.csv, stress_strain_reference.csv,
      and README.md all exist; README.md records a specific
      LAMMPS git tag (regex `lammps[- _]?(stable|patch)`),
      a hardware fingerprint section, and a reproduction
      command sequence.

  - id: ac-006
    summary: MSU launchers exist and accept RANKS arg
    type: code
    pass_when: |
      tools/run_msu_atomiverse.sh and tools/run_msu_lammps.sh
      exist, are executable (mode +x), pass shellcheck, and
      `bash <script> --help` (or no-arg invocation) prints
      a usage line referencing a RANKS positional argument.

  - id: ac-007
    summary: plot_scaling.py emits strong + weak scaling PNGs
    type: code
    pass_when: |
      tools/plot_scaling.py exists; running
      `python tools/plot_scaling.py <fixture.csv>
      --out-dir <tmp>` produces scaling_strong.png and
      scaling_weak.png with non-zero size in <tmp>.

  - id: ac-008
    summary: pe_profile_single benchmark exists with scoped timers
    type: code
    pass_when: |
      benchmarks/md/pe_profile_single.cpp exists and
      contains scoped-timer instrumentation for at least
      these symbols: HaloExchangeForward, HaloExchangeReverse,
      NeighborList, LJ, HarmonicBond, HarmonicAngle,
      MultiHarmonicDihedral, VelocityVerlet, ReduceEnergy,
      SnapshotPublisher.

  - id: ac-009
    summary: Per-component profile CSVs committed at 1 and 8 ranks
    type: code
    pass_when: |
      benchmarks/md/reference/polyethylene/profile_1.csv
      and profile_8.csv exist with header row
      `component,total_ms,pct_of_step,n_calls,mean_us` and
      ≥ 10 data rows each; profile_32.csv exists OR
      scaling_report.md documents "32-rank profile pending
      cluster access".

  - id: ac-010
    summary: plot_profile.py emits stacked-bar PNG from profile CSV
    type: code
    pass_when: |
      tools/plot_profile.py exists; running
      `python tools/plot_profile.py
      benchmarks/md/reference/polyethylene/profile_8.csv`
      produces profile_8.png in the same directory with
      non-zero size.

  - id: ac-011
    summary: LJ kernel SIMD speedup ≥ 1.4× over scalar baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      `./build/benchmarks/md/bench_md_step --force-only=lj`
      run on a 20-chain × 1000-monomer PE box reports a
      mean step time ≤ (1 / 1.4) × the recorded scalar-
      baseline number stored in
      benchmarks/md/reference/polyethylene/scaling_report.md;
      build was configured with -march=native.

  - id: ac-012
    summary: LJ NVE near-bit-identical at 5e-12 after SIMD path (FMA reorder tolerance)
    type: scientific
    pass_when: |
      `ctest -R test_lj_nve --output-on-failure` passes;
      total energy at step 1000 differs from the scalar-
      baseline run by ≤ 5e-12 relative. Rationale: the
      xsimd ket-loop replaces the scalar `a*b + c` chain
      with FMA-fused `fma(a,b,c)` accumulation, which is
      *more* accurate than the scalar path but reorders
      partial sums; bit-identicality at 1e-12 is therefore
      not achievable without disabling FMA (which would
      sacrifice 1.2× of the speedup). 5e-12 is the
      empirically observed maximum per-component force
      deviation on the 20-chain PE box; energy drift is
      still well below the 1e-4 NVE drift gate.

  - id: ac-013
    summary: Bonded ILP speedup ≥ 1.2× over baseline
    type: performance
    evaluator_hint: bench
    pass_when: |
      `./build/benchmarks/md/bench_md_step --force-only=bonded`
      reports mean step time ≤ (1 / 1.2) × the recorded
      pre-B.3 baseline number stored in scaling_report.md.

  - id: ac-014
    summary: Per-term bonded energies bit-identical at 1e-12 after ILP
    type: scientific
    pass_when: |
      `ctest -R "test_pe_single_chain|test_pe_small_box"`
      passes; E_bond, E_angle, E_dihedral at every reported
      step differ from the unpacked-baseline run by ≤ 1e-12
      relative.

  - id: ac-015
    summary: Neighbor-list rebuild count ≤ LAMMPS over 1000 PE-tension steps
    type: performance
    evaluator_hint: bench
    pass_when: |
      profile_1.csv contains a row whose `component` field
      is `NeighborList::build_half` and whose `n_calls`
      value is ≤ the LAMMPS rebuild count recorded in
      scaling_report.md "Rebuild cadence audit" subsection.

  - id: ac-016
    summary: NVE drift on test_pe_small_box ≤ 1e-4 after cadence change
    type: scientific
    pass_when: |
      `ctest -R test_pe_small_box --output-on-failure` passes;
      the embedded NVE-drift assertion |E(1000)-E(0)|/|E(0)|
      < 1e-4 holds.

  - id: ac-017
    summary: HaloExchangeReverse payload audit documented
    type: docs
    pass_when: |
      include/atomiverse/cpu/md/mpi/halo_exchange.h carries
      a doc comment near `HaloExchangeReverse` enumerating
      the exact ghost-row fields included in the reverse
      payload (currently fx, fy, fz) and a one-sentence
      audit conclusion ("no silent loss" or "extended to N
      fields").

  - id: ac-018
    summary: OMP_PROC_BIND warning emitted when unset
    type: runtime
    pass_when: |
      Running `OMP_NUM_THREADS=4 ./build/tests/test_md`
      with OMP_PROC_BIND unset emits a stderr line matching
      regex `OMP_PROC_BIND.*unset` exactly once during
      backend build.

  - id: ac-019
    summary: README documents OMP env + NUMA measurement
    type: docs
    pass_when: |
      benchmarks/md/reference/polyethylene/README.md
      contains a "Thread binding" subsection with at least
      one numerical before/after timing measurement (or an
      explicit "single-socket — no NUMA effect observed"
      note with the measurement command).

  - id: ac-020
    summary: Allocation-free inner loop test passes on Linux
    type: code
    pass_when: |
      tests/cpu/test_md_no_alloc.cpp exists; on Linux CI,
      `ctest -R test_md_no_alloc` passes asserting zero
      malloc/new calls during a 100-step inner loop after
      10 warm-up steps; on Darwin the test reports a
      GTEST_SKIP with a clear "malloc-hook unavailable"
      reason.

  - id: ac-021
    summary: -ffast-math policy documented and gated
    type: docs
    pass_when: |
      CMakeLists.txt contains a comment block (regex
      `-ffast-math.*policy`) stating the policy and uses
      `set_source_files_properties` (or equivalent) to
      restrict -ffast-math to FF-kernel sources; integrator
      sources (velocity_verlet.cpp, nose_hoover.cpp,
      init.cpp) are NOT in the gated set.

  - id: ac-022
    summary: bench_md_step results annotated with build flags in scaling_report.md
    type: docs
    pass_when: |
      scaling_report.md contains a row "single-core PE-
      tension, Δt=0.5fs, 1000 steps" and an adjacent
      annotation listing the exact compiler version and
      flag string used to produce the number.

  - id: ac-023
    summary: Single-core LAMMPS ratio ≤ 1.2× or written-down structural-gap analysis
    type: performance
    evaluator_hint: bench
    pass_when: |
      scaling_report.md row "single-core PE-tension"
      reports either (a) Atomiverse/LAMMPS ratio ≤ 1.2,
      or (b) a ≥ one-paragraph "structural-gap analysis"
      naming the remaining bottleneck and citing its
      profile_1.csv row.

  - id: ac-024
    summary: 8-rank LAMMPS ratio ≤ 2.0× recorded
    type: performance
    evaluator_hint: bench
    pass_when: |
      scaling_report.md row "8-rank PE-tension" reports
      Atomiverse/LAMMPS ratio ≤ 2.0; both numbers were
      produced via tools/run_msu_atomiverse.sh 8 and
      tools/run_msu_lammps.sh 8 on the same machine.

  - id: ac-025
    summary: Full check + test suite green (cpu + mpi labels)
    type: runtime
    pass_when: |
      `cmake --build build -j && ctest --test-dir build
      -L "cpu|mpi" --output-on-failure` exits 0 with no
      regressions vs the pre-Phase-5 test count
      (the 20 existing MPI tests + cpu tests + the two
      new tests test_nose_hoover_t_ramp and test_md_no_alloc
      all green; test_md_no_alloc reports SKIPPED on Darwin).
---

# Acceptance criteria

The 25 criteria above bind Phase 5 closure. They split as:

- **Track A completion (ac-001 .. ac-007):** existence + composition + LAMMPS reference data + tooling.
- **Track B profile (ac-008 .. ac-010):** profile harness, CSVs, plotter.
- **Track B kernel optimizations (ac-011 .. ac-016):** LJ + bonded speedup with bit-identicality guards; neighbor-list cadence audit with no NVE regression.
- **Track B audits + policy (ac-017 .. ac-022):** halo payload, OMP binding, allocation-freeness, -ffast-math gating, scaling report annotations.
- **Track B closure (ac-023 .. ac-024):** the two LAMMPS-ratio targets that define "shippable."
- **Cross-cutting (ac-025):** full suite green.

`runtime` and `performance` criteria require an actual run on the recorded reference machine; `code` and `docs` criteria are static checks against the repo state.
