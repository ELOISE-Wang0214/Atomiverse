#!/usr/bin/env bash
# Single source of truth for the CPU unit gate.
# Used by:
#   - .pre-commit-config.yaml  (stage: pre-push)  — THE gate; blocks a red push
#   - .github/workflows/ci.yml (workflow_dispatch only) — optional remote re-run
#
# Not run on every GitHub push. Multi-backend (CUDA) is machine-local, not GHA.
#
# Unit rule: every ctest case ≤ 10s wall (ATV_TEST_TIMEOUT_SECONDS). Slow
# method campaigns are not unit tests — keep them out of tests/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

jobs="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
BUILD_DIR="${ATV_CI_BUILD_DIR:-build}"
# Match CMake ATV_TEST_TIMEOUT_SECONDS / OMP defaults — keep the gate tight.
TEST_TIMEOUT="${ATV_TEST_TIMEOUT_SECONDS:-10}"
export OMP_NUM_THREADS="${ATV_TEST_OMP_THREADS:-2}"

if [[ -n "${ATV_MOLRS_PATH:-}" && -d "${ATV_MOLRS_PATH}" ]]; then
  MOLRS="$(cd "${ATV_MOLRS_PATH}" && pwd)"
elif [[ -d "${ROOT}/molrs" && -f "${ROOT}/molrs/Cargo.toml" ]]; then
  MOLRS="${ROOT}/molrs"
else
  echo "ci_cpu_unit: molrs missing. Init the submodule:" >&2
  echo "  git submodule update --init --recursive molrs" >&2
  echo "Or set ATV_MOLRS_PATH to a compatible molrs tree." >&2
  exit 1
fi

echo "ci_cpu_unit: molrs=${MOLRS}"
echo "ci_cpu_unit: build=${BUILD_DIR}"
echo "ci_cpu_unit: OMP_NUM_THREADS=${OMP_NUM_THREADS} test_timeout=${TEST_TIMEOUT}s"

cmake -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_CUDA=OFF \
  -DATV_BUILD_TORCH=OFF \
  -DATV_BUILD_TESTS=ON \
  -DATV_BUILD_EXAMPLES=OFF \
  -DATV_BUILD_REGRESSIONS=OFF \
  -DATV_BUILD_BENCHMARKS=OFF \
  -DATV_BUILD_HF=ON \
  -DATV_REAL_STORAGE=double \
  -DATV_TEST_TIMEOUT_SECONDS="${TEST_TIMEOUT}" \
  -DATV_TEST_OMP_THREADS="${OMP_NUM_THREADS}" \
  -DATV_TEST_HEAVY=OFF \
  -DATV_MOLRS_PATH="${MOLRS}"

cmake --build "${BUILD_DIR}" -j"${jobs}"

# Hard kill per-test at TEST_TIMEOUT.
# - Skip MPI multi-rank (needs mpirun).
# - Heavy method tests were not even built (ATV_TEST_HEAVY=OFF).
# - Tree-wide grep "scrub" scripts scan the whole monorepo + molrs and
#   routinely exceed 10s — not component unit tests; run them in a
#   separate architecture job if needed.
ctest --test-dir "${BUILD_DIR}" \
  -j"${jobs}" \
  --timeout "${TEST_TIMEOUT}" \
  --output-on-failure \
  -E 'no_alloc|^test_mpi_|scrub|module_dag|module_off|compatibility_surface|molrs_source_contract'
)
