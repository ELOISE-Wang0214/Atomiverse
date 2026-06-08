# LJ argon (reduced-unit Lennard-Jones) scaling config. Sourced by
# bench_scaling.sbatch. Single source of truth for the LJ benchmark sizes.
#
# AV driver : bench_multi_gpu_scaling  (sc lattice rho*=0.712, rc=2 sigma,
#             skin auto, reneigh 20, NVE)  — prints throughput + VRAM
# LAMMPS    : lj/in.lj_kk  (units lj, lattice sc 0.712, lj/cut 2.0, neigh 20)
# Per-GPU saturation N* (Phase-0 single-card sweep): ~8M atoms/GPU.

NGS="1 2 3 4"
STEPS=200
RENEIGH=20
WEAK_PER_GPU=8000000      # weak: this many atoms per GPU (totals = N* * ng)
STRONG_TOTAL=32000000     # strong: fixed total, split across ng GPUs

# --- CPU scaling (cores). 512/1024 span 2/4 nodes (256 cores/node). ---
CPU_CORES="1 8 32 128 256 512 1024"
CPU_WEAK_PER_CORE=200000      # weak: atoms per rank
CPU_STRONG_TOTAL=8000000      # strong: fixed total
