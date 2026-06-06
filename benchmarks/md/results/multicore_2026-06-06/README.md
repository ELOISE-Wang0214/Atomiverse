# Multi-core CPU scaling: Atomiverse vs LAMMPS (2026-06-06)

Raw weak/strong-scaling data for the CPU MD engine vs LAMMPS, 1 → 1024 cores
(1 → 4 nodes), on the EPYC `cpu` partition (256 cores / 755 GB per node, Cray
Slingshot `cxi` fabric). Cross-node launched with `srun --mpi=pmi2`.

All speeds in **Matom-step/s** (LJ) or **steps/s** (PE). Memory is **max
per-rank resident set** via `/usr/bin/time -v` (full process RSS) for BOTH
codes — a fair comparison (earlier figures had compared AV full RSS vs LAMMPS
self-reported allocation, which overstated the gap).

## Systems

- **LJ** — simple-cubic Lennard-Jones, ρ\*=0.712, rc=2σ (σ=2.5), skin 0.6,
  reneighbor every 20, dt=0.0125 (= LAMMPS `units lj` dt\*=0.005). Matched to
  `benchmarks/md/scaling.cpp` + `lj.in`.
- **PE-melt** — polyethylene melt, Dreiding-UA united-atom (bond/angle/dihedral
  + lj/cut), 20-atom zig-zag chains (`gen_pe.py`), nve/limit + Langevin @ 500 K.
  AV `polyethylene_tension --mode stage:1`; LAMMPS `lammps_bench.in`; both read
  the same gen_pe data file. 8000 atoms/rank weak scaling.

## Data files

| file | contents |
|------|----------|
| `lj_1-128_fairmem.csv`  | LJ weak+strong, 1–128 cores (1 node), fair RSS |
| `lj_256-1024_fairmem.csv` | LJ weak+strong, 256/512/1024 cores (1/2/4 nodes) |
| `pe_1-1024_weak.csv`    | PE-melt weak scaling, 1–1024 cores, fair RSS |
| `lj_1-128.csv`          | LJ 1–128 (earlier run; speed only, old RSS metric) |
| `pe_1-256_prior.csv`    | PE weak 1–256 from the pre-fix run (kept for history) |

## Headline results

**LJ — AV beats LAMMPS at every scale, on-par memory:**
| cores | weak AV/LMP | strong AV/LMP | mem (weak) |
|-------|-------------|---------------|------------|
| 256   | 2.13×       | 1.54×         | 1.7× LMP   |
| 512   | 1.15×       | 1.19×         | 1.4×       |
| 1024  | 1.20×       | 1.22×         | 1.5×       |

**PE-melt — AV competitive ≤64 cores, behind 4–8 % at 128+:**
| cores | AV/LMP speed | AV mem / LMP mem |
|-------|--------------|------------------|
| 1–64  | 1.00–1.05×   | 34→773 / 30→51 MB |
| 128   | 0.94×        | 1541 / 56 MB |
| 256   | 0.92×        | 3082 / 70 MB |
| 1024  | 0.95×        | 12295 / 150 MB |

PE is comm-bound sooner (cheap bonded forces), so the high-core halo gap that
only *thins* AV's LJ lead actually flips it behind for PE. PE memory grows with
global N because `polyethylene_tension` reads the data file on rank 0 +
`scatter` (rank 0 holds the whole frame); the per-rank-local construction that
fixed this for LJ has not yet been applied to PE.

## Engine fixes this round (see git log)

- `d8ea6a5` double-halo (auto-inject 2nd HaloRecv) → momentum leak / blow-up.
- `3665141` scatter int-overflow at >2 GB (segfault at 256+ ranks).
- `4e115e2` per-rank-local LJ construction (no rank holds the global frame).

## Open levers to win PE-melt at high cores

1. **High-core halo/comm efficiency** — close the 4–8 % gap (overlap
   comm/compute, persistent requests). Same lever that thins AV's LJ lead
   cross-node.
2. **PE memory** — per-rank-local / streaming construction so rank 0 doesn't
   hold the global frame (the LJ fix, ported to the bonded chain builder).
