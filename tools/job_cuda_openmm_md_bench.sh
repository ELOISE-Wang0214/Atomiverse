#!/bin/bash -l
#SBATCH -J atv-openmm-md-bench
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gpus 1
#SBATCH -c 16
#SBATCH -t 02:00:00
#SBATCH -o logs/atv-openmm-md-bench-%j.out
#SBATCH -e logs/atv-openmm-md-bench-%j.err

# OpenMM CUDA MD wall timing: N ~ 10k → 1M (cubic LJ liquid, density 0.8).
set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$ROOT"
mkdir -p logs

module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
export OPENBLAS_ROOT="${OPENBLASDIR:-${OPENBLAS_ROOT:-}}"
# GH200 is aarch64 — never force x86_64 rustup (login-node override breaks cmake).
export RUSTUP_TOOLCHAIN=stable-aarch64-unknown-linux-gnu
unset CARGO_BUILD_TARGET || true

echo "=== host $(hostname) $(date -Is) ==="
echo "=== cuda ==="
nvidia-smi -L || true
nvcc --version | head -4 || true
echo "RUSTUP_TOOLCHAIN=$RUSTUP_TOOLCHAIN"
rustup show | head -8 || true

BUILD=build-cuda-bench

# Stale / incomplete tree from a previous failed configure: wipe and redo.
need_cfg=0
if [[ ! -f "$BUILD/CMakeCache.txt" || ! -f "$BUILD/Makefile" && ! -f "$BUILD/build.ninja" ]]; then
  need_cfg=1
fi
if [[ -f "$BUILD/CMakeCache.txt" ]] && grep -q 'x86_64-unknown-linux-gnu' "$BUILD/CMakeCache.txt" 2>/dev/null; then
  echo "=== stale x86 rust/configure in $BUILD — wiping ==="
  rm -rf "$BUILD"
  need_cfg=1
fi
if [[ "$need_cfg" -eq 1 ]]; then
  rm -rf "$BUILD"
  git submodule update --init --recursive molrs || true
  echo "=== cmake configure ==="
  cmake -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DATV_BUILD_CUDA=ON \
    -DATV_BUILD_TESTS=ON \
    -DATV_BUILD_HF=OFF \
    -DATV_BUILD_PW=OFF \
    -DATV_BUILD_GPW=OFF \
    -DBLA_VENDOR=OpenBLAS
fi

echo "=== build cuda_nonbonded_style_bench ==="
cmake --build "$BUILD" -j"$(nproc)" --target cuda_nonbonded_style_bench

BENCH="$BUILD/bin/cuda_nonbonded_style_bench"
if [[ ! -x "$BENCH" ]]; then
  BENCH=$(find "$BUILD" -name cuda_nonbonded_style_bench -type f -executable | head -1)
fi
echo "=== binary: $BENCH ==="
ls -la "$BENCH"

echo
echo "=== OpenMM cluster MD wall timing (LJ NVE) ==="
echo "density=0.8 cut=2.5 skin=0.3; first arg = target N (n_side=ceil(N^{1/3}))"
echo

run_one() {
  local N=$1 steps=$2 warm=$3
  echo "--- target N=$N steps=$steps warm=$warm ---"
  "$BENCH" "$N" "$steps" "$warm"
}

run_one 10000   100 20
run_one 50000    80 15
run_one 100000   50 10
run_one 250000   40 10
run_one 500000   30  8
run_one 1000000  20  5

echo
echo "=== done $(date -Is) ==="
