#!/bin/bash
#SBATCH --account=naiss2025-1-25-cpu
#SBATCH --partition=cpu
#SBATCH --time=00:15:00
export LD_LIBRARY_PATH="/software/sse2/el9_epyc9005/manual/GCC/14.3.0/hpc1/lib64:/software/sse2/el9_epyc9005/manual/OpenBLAS/0.3.32/g14/hpc1/openmp/lp64/lib:/software/sse2/el9_epyc9005/manual/MPICH/4.3.2/g14/hpc2/lib:/software/sse2/el9_epyc9005/manual/FFTW/3.3.11/g14/mp43/hpc1/lib:$LD_LIBRARY_PATH"
MPIEXEC=/software/sse2/el9_epyc9005/manual/MPICH/4.3.2/g14/hpc2/bin/mpiexec
B=/home/jicli594/work/molcrafts/Atomiverse/build/benchmarks/polyethylene_tension
D=/home/jicli594/work/molcrafts/Atomiverse/benchmarks/md/reference/polyethylene/msu_pe.data
OMP_NUM_THREADS=1 $MPIEXEC -n $SLURM_NTASKS $B --out-dir . --data $D --mode full --scale 0.1 --traj-every 0 --print-every 2000 --quiet
