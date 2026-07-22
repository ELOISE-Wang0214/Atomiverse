#!/usr/bin/env bash
# backend-mirror-04-driver gate: the loop node is Driver (not CPUDriver/
# CUDADriver), the erased nestable base is DriverNode, and CUDA aligns its
# termination API to Criterion / add_criterion (no LoopHandler / add_handler).
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
hits="$(grep -rnwE 'CPUDriver|CUDADriver|LoopHandler|add_handler' \
        "${ROOT}/include" "${ROOT}/src" "${ROOT}/tests" "${ROOT}/examples" \
        "${ROOT}/regressions" "${ROOT}/python" "${ROOT}/app" "${ROOT}/tools" \
        2>/dev/null | grep -v 'driver_names_scrub.sh')"
if [ -n "${hits}" ]; then echo "FAIL: legacy driver/handler names:"; echo "${hits}" | head; fail=1; fi
if ! grep -rq 'class DriverNode' "${ROOT}/include/atomiverse/cpu/driver/driver_node.h" 2>/dev/null; then
    echo "FAIL: cpu::DriverNode erased base not found"; fail=1
fi
if [ ${fail} -eq 0 ]; then echo "PASS: Driver/DriverNode names; CUDA Criterion API"; else exit 1; fi
