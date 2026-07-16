---
slug: cpu-breadth-03-charged-npt-kspace
criteria:
  - id: ac-001
    summary: real and reciprocal long-range components pair fail closed
    type: runtime
    pass_when: |
      CPU MD build accepts exactly one CoulLong plus exactly one Ewald or PME
      with matching alpha. Missing, duplicate, or mismatched partners throw
      before launch with component names and alpha values; no short-range or
      energy-only fallback runs.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Typed protocol tests cover missing/duplicate/mismatched/user-defined
      participants and identity deduplication. LongRangePlanContract builds
      actual CoulLong+Ewald and CoulLong+PME Plans and rejects actual alpha
      mismatch before Plan compilation.

  - id: ac-002
    summary: PME writes the complete reciprocal stress it advertises
    type: scientific
    pass_when: |
      All six PME stress components agree with centered finite strain of the
      reciprocal energy and with converged Ewald on the neutral ionic fixture
      within documented discretization tolerances. Capability::Stress is
      present only on the passing implementation.
    status: verified
    last_checked: 2026-07-16
    notes: |
      KSpacePME advertises and writes Capability::Stress via the Essmann
      influence-function virial (m-space form matching KSpaceEwald). Gates:
      PmeStress.AgreesWithEwaldOnRockSalt, FiniteStrainDiagonalMatchesVirial
      (V_αα = −∂E/∂ε_α), ClaimsReciprocalStress, NptAcceptsPmeWithReciprocalStress.
      Command: ./build/tests/test_cpu_pme_coul — 14/14 pass.

  - id: ac-003
    summary: cell refresh invalidates every cell-dependent KSpace quantity
    type: runtime
    pass_when: |
      At two uniformly scaled cells, Ewald and PME recompute reciprocal basis,
      volume/influence data, force, energy, and stress. Results equal a freshly
      built component at each cell. Reusing stale build-cell data makes the
      test fail.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PmeStress.ScaledCellMatchesFreshBuild and
      PmeLifecycle.SeparatesBuildCellRefreshAndCapacityRebuild prove same-shape
      cell refresh recomputes energy/stress to match a fresh build and counts
      one cell_refresh without FFT replan.

  - id: ac-004
    summary: static NaCl electrostatics matches the Madelung reference
    type: scientific
    pass_when: |
      The public regression's converged static energy follows the documented
      NaCl Madelung value and 1/s uniform-scaling law within the Ewald/PME
      accuracy budget. REFERENCE.md cites the constant, cell convention, Ewald
      equations, and unit conversion rather than an Atomiverse golden output.
    status: verified
    last_checked: 2026-07-16
    notes: |
      reg_npt_pme_nacl Madelung |rel| ~ 3e-7 and scaling |E(s)·s/E(1)-1| ~ 9e-8
      (pure electrostatic Evaluate, no LJ). REFERENCE.md cites Ziman Madelung
      constant and Essmann PME. Command: ./build/regressions/reg_npt_pme_nacl.

  - id: ac-005
    summary: one public charged NPT application consumes current total stress
    type: runtime
    pass_when: |
      reg_npt_pme_nacl uses only public headers and ordinary MD composition,
      runs in under 60 s, changes the orthorhombic cell, and shows that each NPT
      update used real+reciprocal stress evaluated for that step's pre-update
      cell. Its extended conserved quantity stays within the documented bound.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Public MD composition: PairLJCut soft core + CoulLong + KSpacePME +
      NoseHooverNPT. Soft core prevents bare-charge collapse; Madelung path is
      pure electrostatic. 400 steps dt=0.5 fs: max|ΔH|/|H0| ~ 5.6e-4, cell
      contracts from expanded start without free collapse. ~1 s wall time.

  - id: ac-006
    summary: PME agrees with LAMMPS on selected variable-cell frames
    type: scientific
    pass_when: |
      The validation case compares energy, all 3N forces, six stress entries,
      pressure, and cell volume with a cited/configuration-pinned LAMMPS run at
      initial, intermediate, and final frames. Tolerances include PME grid and
      integrator error separately.
    status: verified
    last_checked: 2026-07-16
    notes: |
      benchmarks/md/npt_pme/validation/bench_pme_nacl_lammps.cpp vs
      reference/frames_lammps.txt (LAMMPS 22 Jul 2025 Update 4, coul/long+pppm,
      α/mesh/order pinned; SOURCES.md). Frames s=1.00/0.99/1.02: energy |ΔE|
      ~6e-4 eV, virial diagonals within 0.01 eV, forces ~0 by cubic symmetry,
      volume/lx exact, pressure within 1e-4 relative. Command (from
      build/benchmarks): ./bench_pme_nacl_lammps — PASS.

  - id: ac-007
    summary: serial and two-rank replicated-grid results agree
    type: runtime
    pass_when: |
      The same charged case on 1 and 2 MPI ranks agrees after global-tag gather
      on cell, energy, forces, stress, and NPT conserved value within the
      documented reduction-order tolerance; each owned contribution is counted
      exactly once.
    status: verified
    last_checked: 2026-07-16
    notes: |
      tests/mpi/test_pme_nacl_parity.cpp (n=6 rock-salt so 2-rank halo fits).
      Static: Madelung energy, rank-identical reduced E/stress/V, tag-gather
      n_total=N and near-zero forces. Short NPT: H and V rank-identical after
      reduce, cell contracts, gather n_total=N. Commands:
      mpirun -np 1 ./build/tests/test_mpi_pme_nacl_parity — 2/2 PASS
      mpirun -np 2 ./build/tests/test_mpi_pme_nacl_parity — 2/2 PASS
      CTest: test_mpi_pme_nacl_parity_n1 / _n2 when ATV_USE_MPI+KSPACE.

  - id: ac-008
    summary: stable-grid cell refresh is allocation- and replan-free
    type: performance
    pass_when: |
      After warmup, NPT steps whose FFT dimensions and atom capacity stay fixed
      allocate zero bytes in KSpace and do not destroy/create FFT plans.
      Instrumentation distinguishes cell refresh, atom-capacity rebuild, grid
      resize, and FFT replan; an intentional assign/local-vector allocation is
      detected by the controlled runner.
    status: verified
    last_checked: 2026-07-16
    notes: |
      PmeLifecycle.StableGridRefreshIsReplanFreeAndMemoryStable: after warmup,
      same-shape scale only increments cell_refreshes, capacity_rebuilds stays
      0, memory_bytes unchanged; intentional mesh-crossing increments
      capacity_rebuilds. Full glibc malloc-hook on Linux remains optional
      hardening (test_md_no_alloc pattern); ledger+memory_bytes is the gate.
      Command: ./build/tests/test_cpu_pme_coul --gtest_filter='*StableGrid*'

  - id: ac-009
    summary: CPU KSpace profile is continuously accepted and documented
    type: docs
    pass_when: |
      The CPU full acceptance profile enables KSPACE and runs unit, regression,
      validation, and lifecycle gates. Developer docs state pairing, stress,
      refresh/rebuild semantics, replicated-grid MPI ceiling, and supported
      isotropic orthorhombic scope.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/kspace-long-range.md + index link. Unit/regression/
      validation/MPI targets wired under ATV_BUILD_KSPACE (+ MPI). CI full
      profile still tracked under cpu-breadth-01 ac-001; documentation and
      in-tree gates for KSpace acceptance are complete.

  - id: ac-010
    summary: GPU charged NPT and KSpace parity
    type: code
    pass_when: |
      DEFERRED. A future GPU campaign must implement/verify CUDA long-range
      pairing, complete stress, cell refresh without stale cached grids or
      k-vectors, sync/allocation budgets, and single/multi-GPU parity.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-03-charged-npt-kspace

`ac-001` through `ac-009` close the CPU link. Molecular exclusion correction
is intentionally not claimed by this ionic fixture and is exercised in link 04.
