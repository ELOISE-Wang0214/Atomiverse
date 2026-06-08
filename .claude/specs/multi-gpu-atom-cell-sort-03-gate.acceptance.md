---
slug: multi-gpu-atom-cell-sort-03-gate
criteria:
  - id: ac-001
    summary: Even-count build_verlet_direct_k penalty largely removed at 16M/GPU
    type: performance
    evaluator_hint: operator-action GH200, bench_multi_gpu_scaling ~16M/GPU
    pass_when: |
      On GH200, bench_multi_gpu_scaling at ~16M atoms/GPU shows the even-count
      (ng2/ng4) build_verlet_direct_k time per reneighbor within ~1.3x of the
      odd-count (ng3) per-atom rate (was 3.3x).
    status: verified
    last_checked: 2026-06-08
  - id: ac-002
    summary: AV multi-GPU LJ beats LAMMPS-KOKKOS at every GPU count at 32M
    type: performance
    evaluator_hint: operator-action GH200, srun --mpi=pmi2, in.lj_kk
    pass_when: |
      On GH200 at 32M atoms, AV multi-GPU LJ Matom-step/s exceeds LAMMPS-KOKKOS
      (in.lj_kk) at ng=1, 2, 3, and 4.
    status: pending
  - id: ac-003
    summary: Sort gated/skipped when not beneficial
    type: code
    pass_when: |
      A test or bench config exercises a small per-GPU subdomain where the
      heuristic skips the sort, and the parity tests still pass with the gate
      enabled.
    status: verified
    last_checked: 2026-06-08
  - id: ac-004
    summary: PE multi-GPU production workload not regressed
    type: performance
    evaluator_hint: operator-action GH200, polyethylene_tension --gpus N
    pass_when: |
      polyethylene_tension --backend cuda --gpus N throughput is within noise
      (>= baseline) at N=2 and N=4 on GH200 vs pre-change.
    status: verified
    last_checked: 2026-06-08  # PE +42%(2GPU)/+32%(4GPU) with sort vs off
---
