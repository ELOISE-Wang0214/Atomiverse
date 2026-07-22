#!/usr/bin/env bash
# backend-mirror-01-core gate: the facade/ layer is gone and nothing in the
# code trees includes atomiverse/facade/. Contracts live under core/.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
if [ -d "${ROOT}/include/atomiverse/facade" ]; then
    echo "FAIL: include/atomiverse/facade/ still exists"; fail=1
fi
hits="$(grep -rn 'atomiverse/facade/' \
    "${ROOT}/include" "${ROOT}/src" "${ROOT}/tests" "${ROOT}/examples" \
    "${ROOT}/regressions" "${ROOT}/python" "${ROOT}/app" "${ROOT}/tools" \
    2>/dev/null | grep -v 'no_facade_scrub.sh')"
if [ -n "${hits}" ]; then
    echo "FAIL: lingering atomiverse/facade/ includes:"; echo "${hits}" | head -20; fail=1
fi
if [ ${fail} -eq 0 ]; then echo "PASS: no facade references"; else exit 1; fi
