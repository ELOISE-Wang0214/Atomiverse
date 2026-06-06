# CPU multicore scaling — Atomiverse beats LAMMPS (2026-06-06)

Latest CPU multicore weak-scaling results (LJ + PE-melt) after two
correctness-preserving optimizations (commits `17032e9`, `5c47ef6` on `dev`).
This dir supersedes the old `multicore_2026-06-06` (which compared the PE
system with a wrong 18-atom-chain build).

## Weak scaling — atom count grows with cores

These are **weak-scaling** runs: each rank holds a fixed slice, so the *total*
atom count grows with the core count (annotated on the top x-axis of the
figure):

- **LJ**: ~16 000 atoms/rank → 15 625 (1 core) … 16.4 M (1024 cores).
- **PE**: 8 000 atoms/rank → 8 000 (1 core) … 8.19 M (1024 cores).

## The two optimizations

1. **Single-thread direct-write force path** (`17032e9`). For pure-MPI runs
   (one rank/core, OMP=1) forces accumulate straight into the frame arrays
   instead of per-thread scratch — eliding the per-step zero-fill + reduce
   sweep (two O(n_force) passes/kernel LAMMPS never pays). Bonded −26 %,
   Pair −2.7 %. Helps both LJ (pair) and PE (pair + 3 bonded kernels).
2. **Opt-in global-Allreduce cadence** (`5c47ef6`). The per-step 18-double
   energy/kinetic/stress `MPI_Allreduce` is pure diagnostic overhead in
   NVE+Langevin; `MD::set_scalar_reduce_every(N)` fires it only on sampled
   steps (gated on the all-ranks-identical `ctx.step` — no deadlock). NoseHoover
   NVT/NPT stages keep every-step reduction. PE Reduce −2.9 %.

Both verified physics-preserving: full `cpu_md` (40) + MPI suite (10, incl.
NVT/NPT multirank + NVE→NPT `pe_msu_pipeline`) green at OMP=1 and OMP=4; PE
thermo byte-identical at every sampled step.

## Result (weak-scaling throughput, AV / LAMMPS)

| cores | 1 | 8 | 64 | 256 | 512 | 1024 |
|-------|---|---|----|-----|-----|------|
| **LJ**  | 3.0× | 2.5× | 1.8× | 2.2× | 1.2× | 1.2× |
| **PE**  | 1.04× | 1.05× | 1.03× | 1.03× | 1.03× | 1.02× |

PE flipped from **0.92× (8 % behind)** at 256 cores to **1.02–1.06× ahead**
across the whole range. LJ leads everywhere (1.2–3.0×). Peak per-rank memory
(fair `/usr/bin/time -v` max RSS, both codes) is within ~1.5× of LAMMPS.

## Files

| file | contents |
|------|----------|
| `lj_pe_weak_matched.csv` | LJ + PE weak scaling 1→1024 cores: throughput + max RSS/rank, both codes |
| `pe_256_ab.csv` | 256-core PE per-component A/B: baseline → +direct-write → +reduce-gate → LAMMPS |
| `plot_multicore_opt.py` / `multicore_weak_opt.png` | the combined throughput + memory figure (atoms on top x-axis) |

## Method

EPYC `cpu` partition (256 cores/node), 1→4 nodes, OMP=1, `srun --mpi=pmi2`.
LJ: `scaling --atoms-per-rank 16000` vs LAMMPS `lj.in` (ρ*=0.712, rc=2σ,
skin 0.6, reneigh 20). PE: `polyethylene_tension --pe-melt-grid …
--n-atoms-per-chain 20 --mode stage:1` vs LAMMPS `lammps_bench.in` on the
matched `pe_NxNxN.data`. PE throughput from the clean sweep; PE memory and all
LJ data from the `/usr/bin/time -v` sweep (two contention-stragglers re-run
best-of-3). **Gotcha:** AV `--pe-melt-grid` defaults to 18-atom chains; always
pass `--n-atoms-per-chain 20` to match the LAMMPS gen_pe data.

## Remaining lever

PE Comm is still ~1.8× LAMMPS at 256 cores (3-axis cascade halo vs pipelined
6-swap `comm_brick`) — would widen the lead; not required to beat LAMMPS.
