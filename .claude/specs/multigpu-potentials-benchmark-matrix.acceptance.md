---
slug: multigpu-potentials-benchmark-matrix
criteria:
  - id: ac-001
    summary: Multi-GPU domain decomposition carries bonded topology + cross-domain tag resolution
    type: code
    pass_when: |
      The multi-GPU DD/halo ships, for every ghost atom, the global tags of its
      bond/angle/dihedral members, and exposes a bounded tag->local resolve that
      returns -1 for partners absent on a device (analog of the CPU ResolveView).
      A unit test resolves a known cross-device bonded tuple correctly.
    status: pending
  - id: ac-002
    summary: Cross-GPU bonded forces — multi-GPU PE energy matches single-GPU within tolerance
    type: scientific
    evaluator_hint: operator-action on aarch64 GH200 (ATV_BUILD_CUDA=ON), >=2 GPUs
    pass_when: |
      test_cuda_multi_gpu_bonded: for the gen_pe melt, the 2- and 4-GPU energy
      decomposition (E_pair/E_bond/E_angle/E_dihedral), kinetic and virial match
      the single-GPU run from the same initial frame to <=1e-6 relative (or the
      documented mixed-precision tolerance), AND match the --backend cpu single
      reference.
    status: pending
  - id: ac-003
    summary: reverse-once / momentum conservation on multi-GPU PE
    type: scientific
    pass_when: |
      Sum of forces over all atoms stays ~0 each step on 2/4-GPU PE (bonded ghost
      forces are reversed home exactly once, sharing the single halo with LJ); no
      energy drift attributable to double force application.
    status: pending
  - id: ac-004
    summary: existing LJ multi-GPU tests stay green
    type: code
    pass_when: |
      test_cuda_multi_gpu (LJ, 13/13) and the no-d2h / special-bonds CUDA tests
      remain green after the bonded multi-GPU work lands.
    status: pending
  - id: ac-005
    summary: Multi-GPU PE benchmark entry emits throughput + per-card VRAM
    type: code
    pass_when: |
      A PE multi-GPU benchmark path exists (bench_multi_gpu_scaling --system pe
      or bench_multi_gpu_pe) that builds under the ATV_BUILD_CUDA gate and emits
      a CSV row per (atoms,gpus) with matom_step_s and peak VRAM per card.
    status: pending
  - id: ac-006
    summary: Raw-data schema + aggregation with mean/std and outlier rule
    type: code
    pass_when: |
      raw.csv carries one row per (code,system,platform,scaling,point,rep) with
      rep in {1,2,3}; aggregate.py groups and writes summary.csv with
      matom_step_s_mean/_std, mem_mb_mean/_std, n_rep, n_dropped; the outlier
      rule (wall_s > 1.5x median dropped, counted) is implemented and documented.
    status: pending
  - id: ac-007
    summary: Complete CPU matrix (C1-C4) collected, 3 reps each
    type: performance
    evaluator_hint: operator-action — EPYC cpu partition SLURM, srun --mpi=pmi2, up to 4 nodes
    pass_when: |
      raw.csv contains 3 reps for every point of C1 (LJ weak), C2 (LJ strong),
      C3 (PE weak), C4 (PE strong), 1->1024 cores, both AV and LAMMPS; summary.csv
      reports mean+/-std with steady-state std/mean < ~5% (or contention noted).
    status: pending
  - id: ac-008
    summary: Complete single-GPU + multi-GPU-LJ matrix (G1-G4) collected, 3 reps
    type: performance
    evaluator_hint: operator-action — GH200 gpu partition (naiss2026-4-715-gpu)
    pass_when: |
      raw.csv contains 3 reps for G1 (LJ single size-scaling), G2 (LJ multi weak),
      G3 (LJ multi strong), G4 (PE single size-scaling), both AV and LAMMPS-KOKKOS;
      summary.csv with mean+/-std.
    status: pending
  - id: ac-009
    summary: Multi-GPU PE matrix (G5-G6) collected vs LAMMPS-KOKKOS, 3 reps
    type: performance
    evaluator_hint: operator-action — GH200, depends on Phase 1 (ac-002..004)
    pass_when: |
      raw.csv contains 3 reps for G5 (PE multi weak) and G6 (PE multi strong),
      1/2/4 GPU, AV vs LAMMPS-KOKKOS (-sf kk; dihedral multi/harmonic/kk confirmed
      on device, so the comparison is fair); summary.csv with mean+/-std.
    status: pending
  - id: ac-010
    summary: Multi-GPU bonded optimization shows a measured speedup
    type: performance
    evaluator_hint: operator-action — GH200, before/after profile
    pass_when: |
      The Phase-1 optimizations (merged P2P, launch bounds, owner-grouped atomics,
      shared halo) yield a measured multi-GPU PE step-time improvement over the
      first correct version, reported with before/after numbers.
    status: pending
  - id: ac-011
    summary: Unified-style figure with std error bars + README
    type: code
    pass_when: |
      The campaign figure reuses the summary unified style (throughput left/solid,
      memory right/dotted, blue=AV/red=LAMMPS) with error bars = std; a README
      documents the matrix, the 3-rep/outlier method, and the headline mean+/-std.
    status: pending
---

# Acceptance criteria

- **ac-001** — multi-GPU DD ghosts carry bonded tags + bounded cross-device resolve (code).
- **ac-002** — cross-GPU PE energy == single-GPU == CPU reference within tolerance (scientific; operator GH200).
- **ac-003** — reverse-once: ΣF≈0 on multi-GPU PE, no double-application drift (scientific).
- **ac-004** — existing LJ multi-GPU tests stay green (code).
- **ac-005** — PE multi-GPU bench entry emits throughput + per-card VRAM (code).
- **ac-006** — raw 3-rep schema + aggregate.py mean/std + documented outlier rule (code).
- **ac-007** — CPU matrix C1–C4, 3 reps, mean±std (performance; operator EPYC SLURM).
- **ac-008** — GPU matrix G1–G4, 3 reps, mean±std (performance; operator GH200).
- **ac-009** — multi-GPU PE matrix G5–G6 vs LAMMPS-KOKKOS, 3 reps (performance; depends on Phase 1).
- **ac-010** — multi-GPU bonded optimization measured speedup, before/after (performance).
- **ac-011** — unified-style figure with std error bars + README documenting the method (code).
