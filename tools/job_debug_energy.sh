#!/bin/bash -l
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gpus 1
#SBATCH -c 8
#SBATCH -t 00:25:00
#SBATCH -J atv-e-debug
#SBATCH -o tools/debug_energy_%j.out
#SBATCH -e tools/debug_energy_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export OPENBLAS_ROOT="${OPENBLASDIR:-}"
unset RUSTUP_TOOLCHAIN || true
export RUSTUP_TOOLCHAIN=stable-aarch64-unknown-linux-gnu
export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
BD=build-cuda-bench
cmake --build "$BD" -j"$(nproc)" --target atomiverse test_cuda_scf_rhf
echo '=== RUN ==='
"$BD/tests/test_cuda_scf_rhf" --gtest_filter='RHF.CUDASinglePointEnergy' --gtest_color=no
echo DONE
