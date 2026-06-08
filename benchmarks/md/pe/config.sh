# PE melt (united-atom polyethylene) scaling config. Sourced by
# bench_scaling.sbatch. Single source of truth for the PE benchmark sizes.
#
# AV driver : polyethylene_tension --backend cuda --gpus N --pe-melt-grid
#             --mode stage:1  (lj/cut 10.5 + harmonic bond/angle + multi/harmonic
#             dihedral, special_bonds 0 0 0, reneigh 10)
# LAMMPS    : pe/in.pe_kk  (read_data pe_base.data + replicate nx ny nz — SAME
#             nx,ny,nz grid as AV, so identical structures; MSU Dreiding FF)
# Per-GPU saturation N* (Phase-0): ~4M atoms/GPU. Grids: nx*ny*nz*18 atoms.

NGS="1 2 3 4"
STEPS=200          # LAMMPS timed steps
PESCALE=0.05       # AV: stage1 10000 * scale -> 500 steps (throughput is per-step)
RENEIGH=10
# weak: fixed ~4M/GPU. grid + atom count per ng (atoms = nx*ny*nz*18):
PE_WEAK_GRID_1="60 60 60";  PE_WEAK_ATOMS_1=3888000
PE_WEAK_GRID_2="76 76 76";  PE_WEAK_ATOMS_2=7901568
PE_WEAK_GRID_3="87 87 88";  PE_WEAK_ATOMS_3=11992176
PE_WEAK_GRID_4="96 96 96";  PE_WEAK_ATOMS_4=15925248
# strong: fixed total (96^3 ~ 16M) across ng GPUs
PE_STRONG_GRID="96 96 96";  PE_STRONG_ATOMS=15925248

# --- CPU scaling (cores). weak grids ~73728 atoms/core (4096 chains/core) ---
CPU_CORES="1 2 8 32 128 256"
declare -A PE_CPU_WEAK_GRID=( [1]="16 16 16" [2]="16 16 32" [8]="32 32 32" [32]="32 64 64" [128]="64 64 128" [256]="64 128 128" )
declare -A PE_CPU_WEAK_ATOMS=( [1]=73728 [2]=147456 [8]=589824 [32]=2359296 [128]=9437184 [256]=18874368 )
PE_CPU_STRONG_GRID="38 38 38"; PE_CPU_STRONG_ATOMS=987696
