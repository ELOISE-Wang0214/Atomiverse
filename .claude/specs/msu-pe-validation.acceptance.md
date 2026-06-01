---
slug: msu-pe-validation
criteria:
  - id: ac-001
    summary: CPU Phase 0/A/B/C exit gates remain satisfied (no regression)
    type: runtime
    pass_when: |
      The CPU MD test suite stays green: the single-process
      tests/cpu/md/test_msu_reference_inputs (force-field + thermo schema)
      and the 8-rank tests/mpi/test_pe_msu_pipeline end-to-end stack
      (bonded + LJ special-bonds + halo + migration + NVT + barostat +
      deform) both pass.
    status: passed
  - id: ac-002
    summary: CPU ZeroMomentum zeroes COM linear momentum
    type: code
    pass_when: |
      tests/cpu/md/test_zero_momentum.cpp asserts that after a CPU
      cpu::md::ZeroLinearMomentum (alias ZeroMomentum) launch the
      center-of-mass linear momentum of a synthetic system is zero within
      FP tolerance, multi-rank-correct via MPI_Allreduce; test is green.
    status: passed
  - id: ac-003
    summary: Benchmark emits real per-axis Pxx/Pyy/Pzz and energy decomposition; no faked columns
    type: code
    pass_when: |
      benchmarks/md/polyethylene_tension.cpp writes thermo CSV columns
      Pxx/Pyy/Pzz from the virial-tensor diagonal and E_pair/E_bond/E_angle/
      E_dihedral from per-term accumulation; tools/run_msu_atomiverse.sh no
      longer copies isotropic P into the per-axis columns.
    status: passed
  - id: ac-004
    summary: Benchmark protocol is byte-faithful in physics to lammps_msu_paper.in
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      The benchmark applies a fix-momentum equivalent (ZeroMomentum) in
      stages 2/3/4, uses nve/limit+Langevin in stage 1, NPT iso in stages
      2-4 and NPT aniso y/z in stage 5, with the stage timesteps and step
      counts (stage 1: dt=1.0 fs, 10000 steps; stages 2-4: dt=0.5 fs, 50000
      each; stage 5: dt=1.0 fs, 171800 steps) matching the reference deck.
    status: passed
  - id: ac-005
    summary: CPU Thermo runs off the integration thread via AsyncComponentBase
    type: code
    pass_when: |
      cpu::recorder::Thermo subclasses AsyncComponentBase and routes its
      reduction + I/O through AnalysisExecutor; an integration test confirms
      Thermo output is produced without blocking the integration thread.
    status: passed
  - id: ac-006
    summary: Async-correctness gate -- CPU step time is independent of recorder frequency
    type: performance
    evaluator_hint: mol:perf
    pass_when: |
      The async-correctness micro-benchmark reports CPU steps/s with Thermo
      on (high frequency) within 5% of steps/s with Thermo off.
    status: passed
  - id: ac-007
    summary: AV-CPU completes the full 5-stage protocol on msu_pe.data with correct density
    type: runtime
    pass_when: |
      An operator run of polyethylene_tension --mode full --scale 1.0 on the
      real 10000-atom msu_pe.data (PDC Dardel main partition) completes all 5
      stages with no NaN in any thermo column, and the post-stage-4 density is
      within 5% of the MSU amorphous-PE reference (~0.94 g/cm3 at 100 K).
    status: pending
  - id: ac-008
    summary: Atomiverse-CPU agrees with LAMMPS-CPU on block-averaged ensemble means
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      tools/compare_msu.py compares the Atomiverse-CPU run against the
      lmp_cuda CPU-mode reference run on block-averaged ensemble means
      (Flyvbjerg-Petersen 1989, per msu-pe-align-03): <T>, <E_pot> + energy
      decomposition over stage 4, density, and the three stress-strain
      descriptors (E, sigma_y, hardening slope) agree within the stated
      statistical tolerances (~3% on means, +-5% on descriptors); exits 0.
    status: pending
  - id: ac-009
    summary: PDC Dardel reproducibility artifacts -- CPU config sbatch scripts + README
    type: runtime
    pass_when: |
      tools/run_msu_lammps.sh launches lmp_cuda in CPU mode;
      tools/run_msu_atomiverse.sh and tools/run_msu_pdc.sh submit via
      sbatch/srun with an -A account flag to the main (CPU) partition;
      README.md documents the CPU-config reproduction.
    status: passed
---

# Acceptance criteria

CPU-only scope (descoped to CPU MD on 2026-06-01; the CUDA/GPU and async-CUDA
arms were removed -- see `Out of scope` in `msu-pe-validation.md`). All criteria
restart at `ac-001`; the prior dual-backend contract (21 criteria) is retained
in git history.

- **ac-001** -- no-regression gate: the landed CPU Phase 0/A/B/C work stays green.
  Verified 2026-06-01: `test_msu_reference_inputs` 11 passed / 2 skipped (the 2
  skips only await the real run-output reference CSVs), `test_pe_msu_pipeline`
  8-rank PASSED.
- **ac-002 -- ac-006** -- landed CPU functional + async work: ZeroMomentum,
  real per-axis P + energy decomposition, byte-faithful protocol physics, CPU
  `Thermo` async, and the measurable CPU async-correctness gate.
- **ac-007** -- the headline remaining gate: a full-scale (`--scale 1.0`)
  Atomiverse-CPU 5-stage run on `msu_pe.data` completes with correct density.
  Operator-action; runnable on a CPU allocation, no GPU required. A
  `--scale 0.001` smoke (331 steps) already runs all 5 stages on the real
  10k-atom data with no NaN and the full thermo schema (verified 2026-06-01).
- **ac-008** -- CPU<->LAMMPS-CPU statistical comparison (Phase D-CPU verdict).
  Operator-action; needs an `lmp` CPU reference run. The block-average +
  invariant verdict generator (`tools/compare_msu.py`) is implemented and
  unit-tested.
- **ac-009** -- PDC Dardel reproducibility artifacts for the CPU config (tooling
  landed; the operator reproduction runs themselves remain external).

Phase 0/A/B/C completed items are folded into ac-001 as a no-regression gate.
When ac-007 and ac-008 land, the CPU MD MSU validation is complete and the spec
is eligible for `done`.
