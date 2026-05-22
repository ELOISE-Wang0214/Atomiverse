---
slug: msu-pe-validation
criteria:
  - id: ac-001
    summary: CPU Phase 0/A/B/C exit gates remain satisfied (no regression)
    type: runtime
    pass_when: |
      The full test + check suite passes; the pre-existing Phase 0/A/B/C
      tests (32-rank 1e-10 energy match, A.2-A.5, B.2/B.3, end-to-end
      8-rank pipeline) all stay green after the dual-backend changes.
    status: pending
  - id: ac-002
    summary: ZeroMomentum (CPU+CUDA) zeroes COM linear momentum
    type: code
    pass_when: |
      tests/cpu/md/test_zero_momentum.cpp and tests/cuda/md/test_cuda_md_components.cpp
      assert that after a ZeroMomentum launch the center-of-mass linear
      momentum of a synthetic system is zero within FP tolerance, on both
      backends; tests are green.
    status: pending
  - id: ac-003
    summary: CUDA NVELimit caps per-step displacement and survives the overlapped data file
    type: code
    pass_when: |
      tests/cuda/md/test_cuda_md_components.cpp asserts cuda::md::NVELimit
      caps any atom displacement at the configured bound (0.05 A) and a
      short run on the overlapped msu_pe.data configuration produces no
      NaN; test is green.
    status: pending
  - id: ac-004
    summary: CUDA FixLangevin drives temperature to setpoint with a fixed seed
    type: code
    pass_when: |
      tests/cuda/md/test_cuda_md_components.cpp asserts cuda::md::FixLangevin
      with a fixed RNG seed drives a synthetic system's temperature to the
      500 K setpoint within a statistical tolerance over a short run; test
      is green.
    status: pending
  - id: ac-005
    summary: CUDA NoseHooverNPT relaxes box to target pressure in iso and aniso y/z modes
    type: code
    pass_when: |
      tests/cuda/md/test_cuda_md_components.cpp asserts cuda::md::NoseHooverNPT
      relaxes a synthetic box toward the target pressure in both isotropic
      and anisotropic y/z-only coupling modes and honors a temperature ramp;
      test is green.
    status: pending
  - id: ac-006
    summary: CUDA LJ applies 1-2/1-3/1-4 special-bonds exclusions
    type: code
    pass_when: |
      tests/cuda/md/test_cuda_md_components.cpp asserts the CUDA LJ kernel
      with set_special_bonds excludes bonded 1-2/1-3/1-4 pairs (energy/force
      match a CPU reference with the same exclusions), and with an empty
      exclusion list matches the un-masked kernel; test is green.
    status: pending
  - id: ac-007
    summary: CUDA categorical setter API reaches CPU parity and a smoke MD run completes
    type: code
    pass_when: |
      cuda::md::MD exposes set_pair/set_bond/set_angle/set_dihedral/
      set_neighbor/set_integrator/set_thermostat/set_barostat/set_observer/
      set_record/set_thermo; a CUDA MD smoke test built via the categorical
      API compiles and completes a short run without error.
    status: pending
  - id: ac-008
    summary: Benchmark emits real per-axis Pxx/Pyy/Pzz and energy decomposition; no faked columns
    type: code
    pass_when: |
      benchmarks/md/polyethylene_tension.cpp writes thermo CSV columns
      Pxx/Pyy/Pzz from the virial-tensor diagonal and E_pair/E_bond/E_angle/
      E_dihedral from per-term accumulation; tools/run_msu_atomiverse.sh no
      longer copies isotropic P into the per-axis columns.
    status: pending
  - id: ac-009
    summary: Benchmark protocol is byte-faithful in physics to lammps_msu_paper.in
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      The benchmark applies a fix-momentum equivalent (ZeroMomentum) in
      stages 2/3/4, uses nve/limit+Langevin in stage 1, NPT iso in stages
      2-4 and NPT aniso y/z in stage 5, with the stage timesteps and step
      counts (stage 5: dt=1.0 fs, 171800 steps) matching the reference deck.
    status: pending
  - id: ac-010
    summary: CUDA path of polyethylene_tension.cpp builds and runs a short protocol without NaN
    type: runtime
    pass_when: |
      An #ifdef-guarded CUDA path in benchmarks/md/polyethylene_tension.cpp
      builds, and a short (<=1000-step) 5-stage run on a single GPU reading
      msu_pe.data completes with no NaN in any thermo column.
    status: pending
  - id: ac-011
    summary: CPU Thermo runs off the integration thread via AsyncComponentBase
    type: code
    pass_when: |
      cpu::recorder::Thermo subclasses AsyncComponentBase and routes its
      reduction + I/O through AnalysisExecutor; an integration test confirms
      Thermo output is produced without blocking the integration thread.
    status: pending
  - id: ac-012
    summary: CUDA driver does not sync defer_stream on the hot path
    type: code
    pass_when: |
      src/cuda/driver.cu Backend::step() no longer calls
      cudaStreamSynchronize(defer_stream) per phase slot; deferred recorder/
      observer work lags behind the compute stream and a bounded backlog
      drains at run teardown; recorder correctness tests stay green.
    status: pending
  - id: ac-013
    summary: CUDA VirialStress is deferrable when observed-only, hot-path when consumed by a barostat
    type: code
    pass_when: |
      cuda::VirialStress launches on the compute stream when a same-step
      barostat consumes its output, and on defer_stream when stress is only
      observed/recorded; covered by a test asserting both code paths.
    status: pending
  - id: ac-014
    summary: Async-correctness gate — GPU step time is independent of recorder frequency
    type: performance
    evaluator_hint: mol:perf
    pass_when: |
      The Phase H.4 micro-benchmark reports CUDA steps/s with recorders/
      observers at high frequency >= 95% of steps/s with recorders off; the
      CPU counterpart reports Thermo-on within 5% of Thermo-off.
    status: pending
  - id: ac-015
    summary: CUDA force field matches CPU double-precision energies within single-precision tolerance
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      tests/cuda/md/test_pe_cuda_ff.cpp confirms each of E_pair/E_bond/
      E_angle/E_dihedral on the C20H42 perturbed geometry agrees with the
      CPU double-precision reference to <= 1e-3 relative error.
    status: pending
  - id: ac-016
    summary: CUDA completes the full 5-stage protocol on a single GPU with correct density
    type: runtime
    pass_when: |
      An operator run of the CUDA benchmark path on the PDC gpugh partition
      completes all 5 stages on msu_pe.data with no NaN, and the post-stage-4
      density is within 5% of the MSU amorphous-PE reference (~0.94 g/cm3).
    status: pending
  - id: ac-017
    summary: Atomiverse-CUDA and Atomiverse-CPU agree on block-averaged ensemble means
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      tools/validate_polyethylene.py CPU<->CUDA cross-check track reports
      block-averaged <T>, <E_pot>, energy decomposition and density agreeing
      within ~2-3% relative and the three stress-strain descriptors within
      +-5%, and exits 0.
    status: pending
  - id: ac-018
    summary: Atomiverse-CUDA agrees with lmp_cuda GPU-mode on block-averaged ensemble means
    type: scientific
    evaluator_hint: mol:litrev
    pass_when: |
      tools/validate_polyethylene.py GPU-reference track compares the
      Atomiverse-CUDA run against the lmp_cuda GPU-mode reference run on
      block-averaged ensemble means (statistical methodology per
      msu-pe-align-03) and exits 0 within the stated tolerances.
    status: pending
  - id: ac-019
    summary: validate_polyethylene.py runs all three 2x2 tracks multi-engine
    type: runtime
    pass_when: |
      tools/validate_polyethylene.py accepts the four engine CSVs and runs
      the CPU-reference, GPU-reference and CPU<->CUDA cross-check tracks
      (with energy-decomposition columns), keeping the msu-pe-align-03
      block-average methodology; tools/test_validate_polyethylene.py covers
      all three tracks and is green.
    status: pending
  - id: ac-020
    summary: Per-stage wall-clock timing is compared across all four engine configs
    type: runtime
    pass_when: |
      benchmarks/md/polyethylene_tension.cpp emits a per-stage wall-clock
      timing CSV and tools/validate_polyethylene.py --check timing produces a
      four-config (AV-CPU / AV-CUDA / lmp_cuda CPU-mode / lmp_cuda GPU-mode)
      per-stage timing table without error.
    status: pending
  - id: ac-021
    summary: PDC Dardel reproducibility artifacts — lmp_cuda + gpugh sbatch scripts
    type: runtime
    pass_when: |
      tools/run_msu_lammps.sh launches lmp_cuda in CPU mode and GPU mode;
      tools/run_msu_atomiverse.sh and tools/run_msu_pdc.sh submit via
      sbatch/srun with an -A account flag to the main (CPU) and gpugh (GPU)
      partitions; README.md documents the four-config reproduction.
    status: pending
---

# Acceptance criteria

This contract is regenerated fresh for the dual-backend supersede of `msu-pe-validation`; all criteria restart at `ac-001` with `status: pending`.

- **ac-001** — no-regression gate: the landed CPU Phase 0/A/B/C work must stay green after the dual-backend changes.
- **ac-002 – ac-007** — Phase G functional gaps, unit-level (`code`): the missing CUDA components and the categorical setter API.
- **ac-008 – ac-009** — Phase G reference-fidelity defects: real per-axis P and energy decomposition (`code`); byte-faithful protocol physics (`scientific`).
- **ac-010** — Phase G.9 CUDA benchmark path smoke run (`runtime`).
- **ac-011 – ac-014** — Phase H async correctness: CPU `Thermo` async, CUDA driver hot-path desync, deferrable `VirialStress` (`code`), and the measurable async-correctness gate (`performance`).
- **ac-015 – ac-018** — Phase I CUDA single-GPU correctness: force-field correctness (`scientific`), full 5-stage completion (`runtime`), and the CUDA↔CPU / CUDA↔LAMMPS-GPU statistical comparisons (`scientific`).
- **ac-019 – ac-020** — Phase D.5 / E.8: the multi-engine 2×2 validator and per-stage timing (`runtime`).
- **ac-021** — Phase F PDC Dardel reproducibility artifacts (`runtime`).

Phase 0/A/B/C completed items are folded into ac-001 as a no-regression gate rather than re-listed as fresh criteria — they are already verified by their original (now landed) work; the supersede must not silently break them. The Phase D-CPU arm criteria are owned by the `msu-pe-phase-d-*` / `msu-pe-align-03` sub-spec acceptance contracts and are not duplicated here.
