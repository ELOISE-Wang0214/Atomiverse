#!/bin/bash
#SBATCH --account=naiss2025-1-25-cpu
#SBATCH --partition=cpu
#SBATCH --time=04:00:00

# Use MPICH's slurm-aware Hydra launcher (mpiexec): this cluster's MpiDefault is
# "none", so bare `srun` leaves MPICH ranks unwired (every rank a 1-process
# world). Same launcher as job_lammps.sbatch; the binary is MPICH-linked.
MPIEXEC=/software/sse2/el9_epyc9005/manual/MPICH/4.3.2/g14/hpc2/bin/mpiexec
OMP_NUM_THREADS=1 $MPIEXEC -n $SLURM_NTASKS /home/jicli594/work/molcrafts/Atomiverse/build-cpu/benchmarks/polyethylene_tension --out-dir . --data /home/jicli594/work/molcrafts/Atomiverse/benchmarks/md/reference/polyethylene/msu_pe.data \
  --mode full --traj-every 0 --print-every 100
