#!/bin/bash -l
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gpus 1
#SBATCH -c 8
#SBATCH -t 00:45:00
#SBATCH -J atv-bench-rhf
#SBATCH -o tools/bench_cuda_rhf_%j.out
#SBATCH -e tools/bench_cuda_rhf_%j.err

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
cmake --build "$BD" -j"$(nproc)" --target atomiverse test_cuda_scf_rhf

export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
TEST="$BD/tests/test_cuda_scf_rhf"

echo '=== functional (energy + forces + FD) ==='
"$TEST" \
  --gtest_filter='RHF.CUDASinglePointEnergy:RHF.CUDAForcesMatchCPU_H2_STO3G:RHF.CUDAForcesMatchFD_H2_STO3G' \
  --gtest_color=no

echo '=== wall clock (gtest_repeat, warm process) ==='
# First process pays CUDA context + fatbin load; repeat measures steady state.
/usr/bin/time -f 'energy x20 wall %e s' \
  "$TEST" --gtest_filter='RHF.CUDASinglePointEnergy' --gtest_repeat=20 --gtest_color=no
/usr/bin/time -f 'forcesCPU_parity x10 wall %e s' \
  "$TEST" --gtest_filter='RHF.CUDAForcesMatchCPU_H2_STO3G' --gtest_repeat=10 --gtest_color=no

echo '=== microbench via nvcc device-link ==='
INC=(
  -Iinclude -I"$BD/include" -I"$BD/_deps/xtensor-src/include"
  -I"$BD/_deps/xtl-src/include" -I"$BD/_deps/xtensor-blas-src/include"
  -I"$BD/_deps/xsimd-src/include" -I"$BD/_deps/spdlog-src/include"
  -I"$BD/corrosion_generated/cxxbridge/molrs_cxx/include"
)
OB=$(grep -oE '/[^ ]+libopenblas\.so' \
  "$BD/tests/CMakeFiles/test_cuda_scf_rhf.dir/link.txt" | head -1)
CUDA_LIBDIR=$(dirname "$(grep -oE '/[^ ]+/lib64/libcusolver\.so' \
  "$BD/tests/CMakeFiles/test_cuda_scf_rhf.dir/link.txt" | head -1)")

# Separable compilation + explicit device link against the static CUDA archive.
nvcc -O3 -std=c++20 -arch=sm_90 -rdc=true -x cu \
  -DATV_HAS_HF -DATV_HAS_CUDA \
  "${INC[@]}" \
  tools/bench_cuda_rhf.cu -c -o /tmp/bench_cuda_rhf.o

nvcc -O3 -arch=sm_90 -dlink /tmp/bench_cuda_rhf.o "$BD/libatomiverse.a" \
  -o /tmp/bench_cuda_rhf_dlink.o \
  2>&1 | tee /tmp/bench_dlink.log | tail -20 || true

if [[ -f /tmp/bench_cuda_rhf_dlink.o ]]; then
  c++ -O3 -DNDEBUG /tmp/bench_cuda_rhf.o /tmp/bench_cuda_rhf_dlink.o \
    -o /tmp/bench_cuda_rhf \
    -Wl,--start-group \
    "$BD/libatomiverse.a" \
    "$BD/_deps/spdlog-build/libspdlog.a" \
    "$OB" \
    -L"$CUDA_LIBDIR" \
    -lcusolver -lcublas -lcusparse -lcublasLt -lnvJitLink \
    -lcudadevrt -lcudart_static \
    "$BD/libmolrs_cxx.a" "$BD/libmolrs_cxxapi.a" \
    -Wl,--end-group \
    -lgomp -lrt -lpthread -ldl -lm \
    2>&1 | tee /tmp/bench_link.log | tail -30 || true
fi

if [[ -x /tmp/bench_cuda_rhf ]]; then
  /tmp/bench_cuda_rhf 40
  if command -v nsys >/dev/null 2>&1; then
    echo '=== nsys microbench ==='
    nsys profile -o tools/nsys_rhf_h2 --force-overwrite true --stats=true \
      /tmp/bench_cuda_rhf 10 2>&1 | tail -120
  fi
else
  echo 'standalone microbench link failed; nsys on gtest instead'
  if command -v nsys >/dev/null 2>&1; then
    nsys profile -o tools/nsys_rhf_h2 --force-overwrite true --stats=true \
      "$TEST" --gtest_filter='RHF.CUDASinglePointEnergy' --gtest_repeat=15 \
      --gtest_color=no 2>&1 | tail -120
  fi
fi

echo DONE
