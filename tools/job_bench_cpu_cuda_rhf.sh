#!/bin/bash -l
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gpus 1
#SBATCH -c 8
#SBATCH -t 00:50:00
#SBATCH -J atv-cpu-cuda-rhf
#SBATCH -o tools/bench_cpu_cuda_rhf_%j.out
#SBATCH -e tools/bench_cpu_cuda_rhf_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

module load GPU/buildenv-gcccuda/2026.03-cu13.0
export OPENBLAS_ROOT="${OPENBLASDIR:-}"
unset RUSTUP_TOOLCHAIN || true
export RUSTUP_TOOLCHAIN=stable-aarch64-unknown-linux-gnu
export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"

BD=build-cuda-bench
cmake -B "$BD" -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_CUDA=ON -DATV_BUILD_HF=ON -DATV_BUILD_TESTS=ON \
  -DATV_BUILD_PW=OFF -DATV_USE_MPI=OFF -DBLA_VENDOR=OpenBLAS \
  -DATV_BUILD_REGRESSIONS=OFF \
  -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build "$BD" -j"$(nproc)" --target bench_cpu_cuda_rhf

export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo '======== OMP=1 (serial CPU leg) ========'
OMP_NUM_THREADS=1 "$BD/bin/bench_cpu_cuda_rhf" 10 6

echo ''
echo '======== OMP=8 (job -c 8) ========'
OMP_NUM_THREADS=8 "$BD/bin/bench_cpu_cuda_rhf" 10 6

echo DONE
