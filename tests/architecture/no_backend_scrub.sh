#!/usr/bin/env bash
# backend-mirror-05-dissolve gate: there is no CPU "Backend" god class. The CPU
# loop runtime is cpu::Runtime, composed from communication::DomainRuntime (DD),
# device::DeviceRuntime (threads/timer), asynchronous::AnalysisRuntime.
# (The CUDA backend is a separate, dead surface — out of this gate's scope.)
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
# No cpu::Backend anywhere in the code trees.
hits="$(grep -rnw 'cpu::Backend' \
        "${ROOT}/include" "${ROOT}/src" "${ROOT}/tests" "${ROOT}/examples" \
        "${ROOT}/regressions" "${ROOT}/python" "${ROOT}/app" "${ROOT}/tools" \
        2>/dev/null | grep -v 'no_backend_scrub.sh')"
if [ -n "${hits}" ]; then echo "FAIL: lingering cpu::Backend:"; echo "${hits}" | head; fail=1; fi
# No god-class definition under the CPU tree.
if grep -rnq 'class Backend' "${ROOT}/include/atomiverse/cpu" 2>/dev/null; then
    echo "FAIL: a class Backend still exists under include/atomiverse/cpu"; fail=1
fi
# The dissolved sub-runtimes live in their directories.
for f in include/atomiverse/cpu/communication/domain_runtime.h \
         include/atomiverse/cpu/device/device_runtime.h \
         include/atomiverse/cpu/asynchronous/analysis_runtime.h; do
    if [ ! -f "${ROOT}/${f}" ]; then echo "FAIL: missing ${f}"; fail=1; fi
done
if [ ${fail} -eq 0 ]; then echo "PASS: Backend dissolved into communication/device/asynchronous runtimes"; else exit 1; fi
