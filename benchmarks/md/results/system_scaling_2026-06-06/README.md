# System scaling benchmark — Atomiverse vs LAMMPS (2026-06-06)

A 2×2 scaling study of the Atomiverse MD engine against LAMMPS as a benchmark,
all on **one shared LJ system** so the CPU and GPU panels are directly
comparable.

![2x2 scaling](sysbench_2x2.png)

## The system

Simple-cubic LJ fluid, identical parameters across every panel and both codes:

| quantity | value |
|----------|-------|
| σ        | 2.5 Å |
| ε        | 1.0 |
| lattice spacing | 1.12 σ → reduced density ρ\* = (1/1.12)³ ≈ 0.712 |
| cutoff   | 2 σ = 5.0 Å |
| steps    | 200 timed (+10 warmup) |

Built by `make_sc_lj_total` in `benchmarks/md/scaling.cpp` (CPU) and
`bench_multi_gpu_scaling.cu` (GPU) — same lattice, jitter and thermal
velocities, so the throughput numbers measure the same physics.

## Panels

- **(a) CPU weak scaling** — 16k atoms/core, 1→128 MPI ranks. Throughput (left,
  log) and peak RSS summed over ranks (right, log dashed).
- **(b) CPU strong scaling** — 1M atoms fixed, 1→128 ranks. Same axes.
- **(c) Single GH200** — throughput and peak VRAM vs system size (32k→64M
  atoms), mixed precision.
- **(d) Multi-GH200 strong scaling** — at the Pareto-optimal single-GPU size
  (16M atoms), 1/2/4 GPUs; AV uses single-process CUDA P2P over NVLink, LAMMPS
  uses KOKKOS/CUDA + MPI (`-pk kokkos gpu/aware on`).

## Hardware / precision

- **CPU:** AMD EPYC 9005 (`cpu` partition), MPICH 4.3.2, `FI_PROVIDER=shm`,
  `OMP_NUM_THREADS=1`, **double precision**.
- **GPU:** single NVIDIA GH200 node (up to 4 GPUs, NVLink). AV CUDA and LAMMPS
  KOKKOS/CUDA both **mixed precision**.

## Engine optimizations applied (this data)

This run uses the optimized LJ-MD path:
- **LAMMPS-style vflag virial gating** (`set_compute_stress(false)`, CPU + CUDA):
  a plain NVE run never consumes the pair virial, so the 6-component stress is
  skipped. On CPU it trims the inner loop; on CUDA it additionally frees the
  6 double accumulator registers under the kernel's launch-bounds budget.
  Verified forces/energy bit-unchanged (n=1/n=8 MPI bruteforce; GH200
  multi-GPU 13/13). Single-GPU **+~17%** across 0.26M–33M atoms.
- **Verlet skin matched to LAMMPS (0.6 Å):** the engine auto-skin was
  `0.3·cutoff = 1.5 Å`, 2.5× fatter than the `neighbor 0.6 bin` LAMMPS uses
  here. The fat skin inflated both the per-step kernel and the list. Matching
  it is the fair, correctness-preserving choice. CPU per-core throughput rose
  1.3–3.3× (single-core weak now **6.46 vs the prior 4.84**, *ahead of* LAMMPS
  4.67).

## Honest reading of the result

- **Per-core / single-device, AV is now competitive.** CPU weak n=1: AV 6.46 vs
  LAMMPS 4.67 Matom-step/s (AV ahead). Single GH200 at 8M: AV 967 vs LAMMPS
  1957 (~2× behind, was 2.36×). The remaining single-GPU gap is the
  uncoalesced gather — atoms are never physically resorted into spatial order
  (KOKKOS resorts; AV only sorts transient neighbor buffers). That is the next
  CUDA lever.
- **The residual CPU gap is parallel scaling, not the kernel.** At 128 ranks
  AV is still ~9× behind (weak 40 vs 366; strong 45 vs 445). Two causes:
  (1) the per-rank frame is **replicated**, so weak-scaling peak RSS explodes to
  ~188 GB at 128 ranks vs LAMMPS <1 GB (log dashed axis); (2) the
  **`HaloSend{cutoff, every=1}` workaround** (full halo rebuild every step,
  the position-only-halo desync fix) costs a per-step count exchange — the
  weak-scaling dip at n=2 (single rank pays no halo). See the project note
  `atomiverse-cpu-mpi-desync-hang`.
- The next CPU compute lever is the neighbor-build over-scan: cells sized at
  `cutoff+skin` with a ±1 stencil scan ~6.7× too many candidates (profile:
  build ~27% of single-rank time). Binning at half that with a cutoff-clipped
  stencil is the open structural fix.
- CPU strong-scaling AV points start at n=4: the n=1/2 single/dual-core 1M-atom
  runs exceed the benchmark time budget and add nothing to the (clear, log-log
  linear) trend.

## Reproduce the figure

```bash
# from this directory, with a Python that has matplotlib + scienceplots
python plot_2x2.py     # reads the *.csv here, writes sysbench_2x2.{png,pdf}
```

## Data files

| file | contents |
|------|----------|
| `cpu_scaling.csv` | CPU weak+strong, AV+LAMMPS: cores, atoms, wall_s, Matom-step/s, peak RSS |
| `av_single.csv` / `lmp_single.csv` | single-GPU sweep vs atom count |
| `av_multi.csv` / `lmp_multi.csv` | multi-GPU (1/2/4) at several sizes |
| `plot_2x2.py` | the plotting script (co-located paths) |
