#!/usr/bin/env bash
# Grep gate: the migrated CUDA backend must carry no `scf::Phase` token.
#
# After the CUDA driver-taxonomy migration (cuda-driver-taxonomy chain,
# links 01-05) CUDA SCF is the `cuda::scf::SCFSolver` inner solver and no
# CUDA driver is parameterized on `ATV::scf::Phase`. This test asserts the
# enum's name does not survive anywhere under the CUDA source tree —
# neither as code nor in doc-comments.
#
# It deliberately scopes ONLY include/atomiverse/cuda/ and src/cuda/.
# The CPU plane-wave SCF (`cpu::pw`) legitimately keeps `scf::Phase` as
# its outer-driver phase enum, so cpu/ is out of scope.
#
# Usage: scf_phase_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: scf_phase_scrub.sh <repo-root>}"

status=0
for dir in "include/atomiverse/cuda" "src/cuda"; do
    target="${ROOT}/${dir}"
    if [ ! -d "${target}" ]; then
        echo "FAIL: expected directory not found: ${target}"
        status=1
        continue
    fi
    matches="$(grep -rn 'scf::Phase' "${target}" 2>/dev/null)"
    if [ -n "${matches}" ]; then
        echo "FAIL: 'scf::Phase' token(s) survive under ${dir}:"
        echo "${matches}"
        status=1
    fi
done

if [ "${status}" -eq 0 ]; then
    echo "PASS: no 'scf::Phase' token under include/atomiverse/cuda/ or src/cuda/"
fi
exit "${status}"
