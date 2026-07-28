#!/usr/bin/env bash
# Regression-harness gate. Two defects made the `regressions` target lie about
# its own results; both are cheap to reintroduce and expensive to notice, so
# they are grep-enforced.
#
#  1. CWD-relative data paths. Seven cases probed "tests/data/...",
#     "../tests/data/...", ... and so PASSED from the repo root while FAILING
#     from the build directory — which is exactly how the aggregate target runs
#     them. Data is resolved through the compiled-in absolute ATV_REG_DATA_DIR
#     (or a per-case ATV_* absolute define), never by searching the CWD.
#
#  2. Silent degradation to a different Hamiltonian. reg_pw-smearing-al cleared
#     upf_path and solved a toy V_loc when the pseudopotential was missing, then
#     printed PASS against a pin recorded for the real UPF. A missing input is a
#     failure; house rule "no fallback values".
#
# Also pins the declared-thread-count harness (see regressions/CMakeLists.txt):
# every case's OpenMP width is declared at registration, so run_reg_* must go
# through `cmake -E env OMP_NUM_THREADS=...`.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
REG="${ROOT}/regressions"
SELF="regression_harness_scrub.sh"
fail=0

# 1. No CWD-relative data path in any case source.
hits="$(grep -rn --include='*.cpp' --include='*.h' \
        -e '"tests/data' -e '"\.\./tests/data' -e '"\.\./\.\./tests/data' \
        -e '"[a-z_]*/\.\./tests/data' \
        "${REG}" 2>/dev/null | grep -v "${SELF}")"
if [ -n "${hits}" ]; then
    echo "FAIL: CWD-relative data path in a regression case — use"
    echo "      ATV_REG_DATA_DIR \"/<subdir>/<file>\" (compiled-in, absolute):"
    echo "${hits}" | head -20
    fail=1
fi

# 2. No environment-variable escape hatch for locating case data. An env var is
#    a fallback by another name: it makes the case pass or fail depending on the
#    shell that launched it.
hits="$(grep -rn --include='*.cpp' \
        -e 'getenv("ATV_SOURCE_DIR")' -e 'getenv("ATV_SI_UPF")' \
        -e 'getenv("ATV_AL_UPF")' \
        "${REG}" 2>/dev/null | grep -v "${SELF}")"
if [ -n "${hits}" ]; then
    echo "FAIL: regression case locates data through the environment:"
    echo "${hits}" | head -20
    fail=1
fi

# 3. No silent substitution of a different physical model for missing input.
hits="$(grep -rn --include='*.cpp' -e 'upf_path.clear()' \
        "${REG}" 2>/dev/null | grep -v "${SELF}")"
if [ -n "${hits}" ]; then
    echo "FAIL: case clears upf_path (silently solves a toy V_loc and still"
    echo "      reports PASS). A missing pseudopotential is a failure:"
    echo "${hits}" | head -20
    fail=1
fi

# 4. Every case is registered through atv_add_regression, which is what supplies
#    ATV_REG_DATA_DIR and the declared OMP_NUM_THREADS. A hand-rolled
#    add_custom_target(run_reg_* COMMAND <binary>) bypasses both.
CM="${REG}/CMakeLists.txt"
if ! grep -q 'ATV_REG_DATA_DIR="\${CMAKE_SOURCE_DIR}/tests/data"' "${CM}"; then
    echo "FAIL: atv_add_regression no longer compiles in ATV_REG_DATA_DIR"
    fail=1
fi
if ! grep -q 'OMP_NUM_THREADS=\${REG_THREADS}' "${CM}"; then
    echo "FAIL: atv_add_regression no longer declares OMP_NUM_THREADS"
    fail=1
fi

if [ ${fail} -eq 0 ]; then
    echo "PASS: regression cases resolve data absolutely, never degrade"
    echo "      silently, and declare their thread width at registration"
else
    exit 1
fi
