#!/usr/bin/env bash
# backend-mirror-03-paths gate: backend subdirectories are full English words,
# mirrored across cpu/ and cuda/. No abbreviated legacy dir names, no legacy
# include paths, no alias directories.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
OLD=(comm async md mc minimizer vib pw observer recorder)
for be in cpu cuda; do
    for old in "${OLD[@]}"; do
        for base in "include/atomiverse/${be}" "src/${be}" "tests/${be}"; do
            if [ -d "${ROOT}/${base}/${old}" ]; then
                echo "FAIL: abbreviated dir ${base}/${old} must be a full word"; fail=1
            fi
        done
    done
done
# cpu must carry the full-word directories
for w in dynamics communication asynchronous montecarlo minimize vibration planewave observe record; do
    if [ ! -d "${ROOT}/include/atomiverse/cpu/${w}" ]; then
        echo "FAIL: missing full-word dir include/atomiverse/cpu/${w}"; fail=1
    fi
done
# no legacy include paths anywhere in the code trees
hits="$(grep -rnEl 'atomiverse/(cpu|cuda)/(comm|async|md|mc|minimizer|vib|pw|observer|recorder)/' \
        "${ROOT}/include" "${ROOT}/src" "${ROOT}/tests" "${ROOT}/examples" \
        "${ROOT}/regressions" "${ROOT}/python" "${ROOT}/app" "${ROOT}/tools" \
        2>/dev/null | grep -v 'path_words_scrub.sh')"
if [ -n "${hits}" ]; then echo "FAIL: legacy include paths in:"; echo "${hits}" | head; fail=1; fi
if [ ${fail} -eq 0 ]; then echo "PASS: full-word mirrored backend directories"; else exit 1; fi
