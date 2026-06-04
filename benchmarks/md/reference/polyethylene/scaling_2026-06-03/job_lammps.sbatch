#!/bin/bash
#SBATCH --account=naiss2025-1-25-cpu
#SBATCH --partition=cpu
#SBATCH --time=04:00:00

# NOTE: this cluster's MpiDefault is "none", so bare `srun` leaves MPICH ranks
# unwired (every rank a 1-process world, all clobbering the same output). Use
# MPICH's own slurm-aware Hydra launcher (mpiexec) for real domain decomposition.
MPIEXEC=/software/sse2/el9_epyc9005/manual/MPICH/4.3.2/g14/hpc2/bin/mpiexec
OMP_NUM_THREADS=1 $MPIEXEC -n $SLURM_NTASKS /home/jicli594/.local/x86_64/bin/lmp -in lammps.in -log log.lammps
