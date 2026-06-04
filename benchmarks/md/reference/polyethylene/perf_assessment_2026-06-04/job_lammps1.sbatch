#!/bin/bash
#SBATCH --account=naiss2025-1-25-cpu
#SBATCH --partition=cpu
#SBATCH --time=00:30:00
OMP_NUM_THREADS=1 srun --mpi=pmi2 /home/jicli594/.local/x86_64/bin/lmp -in lammps_perf.in -log log.lammps_perf
