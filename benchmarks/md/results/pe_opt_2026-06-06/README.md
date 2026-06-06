# PE-melt CPU optimization — Atomiverse now beats LAMMPS (2026-06-06)

Two correctness-preserving optimizations (commits `17032e9`, `5c47ef6` on
`dev`) flipped the polyethylene-melt benchmark from **0.92× LAMMPS** (8 %
behind at 256 cores) to **1.02–1.06× at every scale 1→1024 cores**.

## What changed

1. **Single-thread direct-write force path** (`17032e9`). For pure-MPI runs
   (one rank per core, `OMP_NUM_THREADS=1` — the production multi-node layout)
   the per-thread force scratch is pure overhead. Forces now accumulate
   straight into the frame arrays (already cleared by `ZeroForces` +
   `HaloSend::zero_ghost_forces`), eliding the per-step zero-fill **and** the
   reduce sweep — two `O(n_force)` passes per kernel per step that LAMMPS never
   pays. Applied to the three bonded kernels (`BondedForceScratch`) and
   `PairLJCut`. Multi-thread (hybrid MPI+OMP) keeps the private-buffer + reduce
   path. Effect: **Bonded −26 %**, Pair −2.7 %, Comm −10 %.

2. **Opt-in global-Allreduce cadence** (`5c47ef6`). The end-of-step `Reduce`
   (18-double energy/kinetic/stress `MPI_Allreduce`) fired every step; a plain
   NVE / per-atom-Langevin run consumes none of it between thermo samples
   (LAMMPS likewise only reduces at thermo output). `Reduce(every)` +
   `MD::set_scalar_reduce_every(N)` gate the collective on the all-ranks-
   identical `ctx.step` (lockstep — no deadlock). `polyethylene_tension` stage 1
   sets the cadence to the recorder stride; **stages 2–3 (NoseHoover NVT/NPT)
   keep every-step reduction** — they integrate the global T/P each step.
   Effect: Reduce_s 0.34 → 0.027 s, **−2.9 %**.

Both verified physics-preserving: full `cpu_md` suite (40) + MPI suite (10,
incl. NVT/NPT multirank and the NVE→NPT `pe_msu_pipeline`) green at OMP=1
(direct path) and OMP=4 (scratch path). Thermo at every sampled step is
byte-identical to the every-step-reduce baseline (E_pot=139981, E_pair=−58486.5,
E_bond=46409.3 at step 500 of the 2.048M-atom / 256-rank run).

## Benchmark

Weak scaling, 8000 atoms/rank, gen_pe Dreiding-UA melt (`--pe-melt-grid`,
**20-atom chains** — `--n-atoms-per-chain 20`, matched to LAMMPS
`pe_NxNxN.data`), stage 1 (NVE/limit + Langevin @ 500 K), 1000 steps, EPYC
`cpu` partition (256 cores/node), `srun --mpi=pmi2`, OMP=1. AV
`polyethylene_tension --mode stage:1`; LAMMPS `lammps_bench.in`.

**Gotcha:** AV `--pe-melt-grid` defaults to **18-atom chains**; the LAMMPS
gen_pe data files are **20-atom**. Always pass `--n-atoms-per-chain 20` for a
matched comparison (the earlier figure data compared 18 vs 20 chains).

## Data files

| file | contents |
|------|----------|
| `pe_weak_1-1024_matched.csv` | weak scaling 1→1024 cores, AV vs LAMMPS, Matom-step/s (the headline sweep) |
| `pe_256_ab.csv` | 256-core A/B: baseline → +direct-write → +reduce-gate → LAMMPS, per-component breakdown |
| `plot_pe_opt.py` / `pe_weak_scaling_opt.png` | throughput + speedup-ratio figure |

## Headline (Matom-step/s, AV / LAMMPS)

| cores | 1 | 32 | 64 | 128 | 256 | 512 | 1024 |
|-------|---|----|----|-----|-----|-----|------|
| AV    | 1.4 | 36.6 | 75.1 | 142.6 | 236.4 | 455.1 | 880.0 |
| LMP   | 1.3 | 36.3 | 72.9 | 138.8 | 229.2 | 440.4 | 864.7 |
| ratio | 1.04× | 1.01× | 1.03× | 1.03× | 1.03× | 1.03× | 1.02× |

## Remaining lever

Comm is still ~1.8× LAMMPS at 256 cores (3-axis cascade halo with per-axis
`Waitall` vs LAMMPS's pipelined 6-swap `comm_brick`). Closing it would widen
the lead further; not required to beat LAMMPS.
