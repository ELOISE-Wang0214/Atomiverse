#!/usr/bin/env bash
# Single source of truth for the CPU unit half of CI.
# Used by:
#   - .pre-commit-config.yaml  (stage: pre-push)  — blocks a red push
#   - .github/workflows/ci.yml (cpu job)          — same commands on Actions
#
# Flags match CI exactly. Do not fork a second matrix here.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

jobs="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
BUILD_DIR="${ATV_CI_BUILD_DIR:-build}"

if [[ -n "${ATV_MOLRS_PATH:-}" && -d "${ATV_MOLRS_PATH}" ]]; then
  MOLRS="$(cd "${ATV_MOLRS_PATH}" && pwd)"
elif [[ -d "${ROOT}/molrs" ]]; then
  MOLRS="${ROOT}/molrs"
elif [[ -d "${ROOT}/../molrs" ]]; then
  MOLRS="$(cd "${ROOT}/../molrs" && pwd)"
else
  echo "ci_cpu_unit: set ATV_MOLRS_PATH or place molrs at ./molrs or ../molrs" >&2
  exit 1
fi

echo "ci_cpu_unit: molrs=${MOLRS}"
echo "ci_cpu_unit: build=${BUILD_DIR}"

cmake -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_CUDA=OFF \
  -DATV_BUILD_TORCH=OFF \
  -DATV_BUILD_TESTS=ON \
  -DATV_BUILD_EXAMPLES=OFF \
  -DATV_BUILD_RHF=ON \
  -DATV_REAL_STORAGE=double \
  -DATV_MOLRS_PATH="${MOLRS}"

cmake --build "${BUILD_DIR}" -j"${jobs}"

# tests/ only. no_alloc is runtime-dependent (glibc/libgomp) — not a portable gate.
ctest --test-dir "${BUILD_DIR}" -j"${jobs}" --output-on-failure -E 'no_alloc'
