---
slug: benchmark-matrix
criteria:
  - id: ac-001
    summary: raw-data schema + aggregate.py with mean/std and outlier rule
    type: code
    pass_when: |
      raw.csv has one row per (code,system,platform,scaling,point,rep) with rep
      in {1,2,3}; aggregate.py emits summary.csv with matom_step_s_mean/_std,
      mem_mb_mean/_std, n_rep, n_dropped; the wall_s>1.5x median outlier rule is
      implemented and documented; unit-checked on synthetic rows.
    status: pending
  - id: ac-002
    summary: CPU matrix C1-C4 collected, 3 reps each
    type: performance
    evaluator_hint: operator-action EPYC SLURM, srun --mpi=pmi2, up to 4 nodes
    pass_when: raw.csv has 3 reps for every C1-C4 point (AV+LAMMPS); summary mean±std.
    status: pending
  - id: ac-003
    summary: GPU matrix G1-G4 collected, 3 reps each
    type: performance
    evaluator_hint: operator-action GH200 (naiss2026-4-715-gpu)
    pass_when: raw.csv has 3 reps for G1-G4 (AV+LAMMPS-KOKKOS); summary mean±std.
    status: pending
  - id: ac-004
    summary: multi-GPU PE matrix G5-G6 vs LAMMPS-KOKKOS, 3 reps
    type: performance
    evaluator_hint: operator-action GH200, depends on multigpu-bonded chain
    pass_when: |
      raw.csv has 3 reps for G5 (PE multi weak) and G6 (PE multi strong), 1/2/4
      GPU, AV vs LAMMPS-KOKKOS (-sf kk; multi/harmonic/kk on device); summary mean±std.
    status: pending
  - id: ac-005
    summary: unified-style figure with std error bars + README
    type: code
    pass_when: |
      The figure reuses the summary unified style (throughput left/solid, memory
      right/dotted, blue=AV/red=LAMMPS) with error bars = std; README documents
      the matrix, the 3-rep/outlier method, and headline mean±std.
    status: pending
---

# Acceptance criteria

- **ac-001** — raw schema + aggregate.py mean/std + outlier rule (code).
- **ac-002** — CPU matrix C1–C4, 3 reps (performance; operator EPYC).
- **ac-003** — GPU matrix G1–G4, 3 reps (performance; operator GH200).
- **ac-004** — multi-GPU PE G5–G6, 3 reps (performance; operator GH200, depends on chain).
- **ac-005** — unified figure with std error bars + README (code).
