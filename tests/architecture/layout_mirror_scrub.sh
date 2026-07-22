#!/usr/bin/env bash
# backend-mirror-06-gates: the cross-backend layout invariants that the earlier
# per-phase scrubs don't cover.
#   1. no runtime/ , loop/ , shared/ directories anywhere under include/atomiverse
#      (the loop node is Driver; the coordinator is cpu::Runtime in driver/, a
#      file — never a directory; shared math is backend-owned, not top-level).
#   2. cuda/ mirrors cpu/: every cuda backend subdirectory has a same-named cpu
#      twin (cuda is a subset of the cpu directory vocabulary — a CPU-only domain
#      may be absent from cuda, but cuda must never invent a directory cpu lacks).
#   3. tests mirror source: every tests/cpu/<dir> has include/atomiverse/cpu/<dir>
#      (no stale old-layout test directories).
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0

for bad in $(find "${ROOT}/include/atomiverse" -type d \
             \( -name runtime -o -name loop -o -name shared \) 2>/dev/null); do
    echo "FAIL: forbidden directory ${bad#${ROOT}/}"; fail=1
done

for d in "${ROOT}"/include/atomiverse/cuda/*/; do
    [ -d "$d" ] || continue
    b="$(basename "$d")"
    if [ ! -d "${ROOT}/include/atomiverse/cpu/${b}" ]; then
        echo "FAIL: cuda/${b} has no cpu/ twin (cuda must mirror cpu)"; fail=1
    fi
done

for d in "${ROOT}"/tests/cpu/*/; do
    [ -d "$d" ] || continue
    b="$(basename "$d")"
    if [ ! -d "${ROOT}/include/atomiverse/cpu/${b}" ]; then
        echo "FAIL: tests/cpu/${b} has no include/atomiverse/cpu/${b} source dir"; fail=1
    fi
done

if [ ${fail} -eq 0 ]; then echo "PASS: no runtime/loop/shared dirs; cuda mirrors cpu; tests mirror source"; else exit 1; fi
