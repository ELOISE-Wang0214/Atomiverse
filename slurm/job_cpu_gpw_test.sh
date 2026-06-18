#!/bin/bash -l
#SBATCH -A naiss2026-3-381-cpu
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -J atv_gpw_test
#SBATCH -o slurm/%x-%j.out
#
# Build + run the GPW module's CPU tests (currently the CP2K data-file parsers).
# Submit from the worktree root:  sbatch slurm/job_cpu_gpw_test.sh
set -euo pipefail

module load buildenv-gcc/2026.03-mpich
cd "$SLURM_SUBMIT_DIR"

# Drop the molcrafts venv cmake stub from PATH; keep cargo (molrs crate).
export PATH="$(printf '%s' "$PATH" | tr ':' '\n' \
    | grep -vE '/molcrafts/\.(x86|x86_64|aarch64)/bin/?$' | paste -sd:)"
export PATH="$HOME/.local/x86_64/cargo/bin:$PATH"

CMAKE="$(command -v cmake)"
CTEST="$(dirname "$CMAKE")/ctest"; [ -x "$CTEST" ] || CTEST="$(command -v ctest)"
echo "cmake: $CMAKE ($("$CMAKE" --version | head -1))"

MOLRS=/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/molrs
BUILD=build-cpu-gpw

"$CMAKE" -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DATV_BUILD_CUDA=OFF \
    -DATV_BUILD_TESTS=ON \
    -DATV_BUILD_RHF=ON \
    -DATV_BUILD_GPW=ON \
    -DBLA_VENDOR=OpenBLAS \
    -DCMAKE_PREFIX_PATH="${FFTW_ROOT:-};${OPENBLAS_ROOT:-}" \
    -DATV_MOLRS_PATH="$MOLRS"

"$CMAKE" --build "$BUILD" -j"$(nproc)" \
    --target test_cpu_gpw_cp2k_io test_cpu_gpw_basis \
            test_cpu_gpw_gth_integrals test_cpu_gpw_collocate test_cpu_xc_pade

export LD_LIBRARY_PATH="${LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"
"$CTEST" --test-dir "$BUILD" \
    -R "test_cpu_(gpw_cp2k_io|gpw_basis|gpw_gth_integrals|gpw_collocate|xc_pade)" \
    --output-on-failure
echo "=== JOB DONE rc=$? ==="
