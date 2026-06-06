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

## Engine fixes + optimizations applied (this data)

- **Multi-rank force correctness fix (commit d8ea6a5) — the big one.**
  `MD::install_defaults` auto-injected a HaloSend+HaloRecv even when the caller
  already registered them, so the reverse halo ran TWICE and re-applied the
  already-reduced ghost forces — a per-step momentum leak (Σ F ≠ 0) that pumped
  energy and blew the NVE trajectory up. Fixed by skipping the auto-inject when
  a halo is present. Multi-rank LJ NVE is now **bit-identical to single-rank**
  (n=8 == n=1 = -468123.033977 over 500 steps) with momentum conserved.
- **Stable timestep (commit 99643aa).** The bench used dt=1.0 (dt/τ≈0.4, ~80×
  too large) → the integrator blew up (single-rank too). Set dt=0.0125 to match
  LAMMPS `units lj` / `timestep 0.005` (0.005·τ, τ=σ=2.5). Both codes now run a
  stable, dense system, so the throughput numbers are physical and fair.
- **Reneighbor every 20 + comm cutoff = pair cutoff + skin** (commit 8b6ba8f),
  matching LAMMPS `neigh_modify every 20` / `neighbor 0.6 bin`. The halo + NL
  rebuild only every 20 steps; the per-step cost is the LJ kernel.
- **vflag virial gating** (`set_compute_stress(false)`, CPU + CUDA, ac98851 /
  d854b6f): NVE skips the unused pair virial. Forces/energy bit-unchanged.
  Single-GPU **+~17%**.

## Honest reading of the result

- **CPU multi-core: AV now BEATS LAMMPS at every core count, 1.5–2.9×** (same
  EPYC node, double precision, LAMMPS-matched skin/cadence/timestep):
  | cores | weak AV vs LMP | strong AV vs LMP |
  |---|---|---|
  | 1 | 13.5 vs 4.67 (2.9×) | 13.4 vs 4.61 (2.9×) |
  | 8 | 82.5 vs 34.5 (2.4×) | 90.8 vs 36.3 (2.5×) |
  | 32 | 288 vs 126 (2.3×) | 297 vs 132 (2.3×) |
  | 128 | 617 vs 366 (1.7×) | 794 vs 445 (1.8×) Matom-step/s |
- **AV's remaining CPU weakness is memory, not speed:** the per-rank frame is
  replicated, so weak-scaling peak RSS still reaches ~188 GB at 128 ranks vs
  LAMMPS <1 GB (log dashed axis). Right-sizing per-rank storage is the open
  CPU work — it does not affect the throughput win.
- **Single GH200:** AV 967 vs LAMMPS 1957 Matom-step/s at 8M (~2× behind). The
  gap is the uncoalesced gather — atoms are never physically resorted into
  spatial order (KOKKOS resorts; AV only sorts transient neighbor buffers).
  That is the next CUDA lever.

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
