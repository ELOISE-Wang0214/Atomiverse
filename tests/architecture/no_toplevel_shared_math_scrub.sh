#!/usr/bin/env bash
# backend-mirror-02-absorb gate: no top-level shared-math module. Reciprocal /
# pair kernels are owned per backend; nothing sits at include/atomiverse/
# {longrange,potentials}, and cuda uses the full word longrange, never kspace.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
for d in include/atomiverse/longrange include/atomiverse/potentials \
         src/longrange tests/longrange \
         include/atomiverse/cuda/kspace src/cuda/kspace tests/cuda/kspace; do
    if [ -e "${ROOT}/${d}" ]; then echo "FAIL: ${d} must not exist"; fail=1; fi
done
hits="$(grep -rnF \
        -e 'atomiverse/longrange/' \
        -e 'atomiverse/potentials/' \
        -e 'cuda/kspace/' \
        "${ROOT}/include" "${ROOT}/src" "${ROOT}/tests" "${ROOT}/examples" \
        "${ROOT}/regressions" "${ROOT}/python" "${ROOT}/app" "${ROOT}/tools" \
        2>/dev/null | grep -v 'no_toplevel_shared_math_scrub.sh')"
if [ -n "${hits}" ]; then
    echo "FAIL: forbidden include/path:"; echo "${hits}" | head -20; fail=1
fi
if [ ${fail} -eq 0 ]; then echo "PASS: no top-level shared math, no cuda/kspace"; else exit 1; fi
