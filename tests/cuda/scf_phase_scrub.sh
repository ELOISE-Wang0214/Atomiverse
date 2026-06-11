#!/usr/bin/env bash
# Grep gate: no `scf::Phase` token anywhere in the source tree.
#
# After the CUDA driver-taxonomy migration (cuda-driver-taxonomy chain,
# links 01-05) CUDA SCF became the `cuda::scf::SCFSolver` inner solver,
# and after the cpu::pw migration onto PWSolver/PWEvaluator the CPU
# plane-wave SCF no longer carries an outer phase enum either. The
# `ATV::scf::Phase` enum itself was deleted from facade/stage.h —
# SCF (gaussian or plane-wave, CPU or CUDA) is an inner solver behind
# an Evaluator component and never a Driver. This gate asserts the
# token does not survive as code in include/, src/, examples/, or
# tests/ (the facade/stage.h comment documenting the deletion is the
# single allowed mention).
#
# Usage: scf_phase_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: scf_phase_scrub.sh <repo-root>}"

status=0
for dir in "include/atomiverse" "src" "examples" "tests"; do
    target="${ROOT}/${dir}"
    if [ ! -d "${target}" ]; then
        echo "FAIL: expected directory not found: ${target}"
        status=1
        continue
    fi
    matches="$(grep -rn 'scf::Phase' "${target}" 2>/dev/null \
        | grep -v 'facade/stage.h' \
        | grep -v 'scf_phase_scrub.sh')"
    if [ -n "${matches}" ]; then
        echo "FAIL: 'scf::Phase' token(s) survive under ${dir}:"
        echo "${matches}"
        status=1
    fi
done

if [ "${status}" -eq 0 ]; then
    echo "PASS: no 'scf::Phase' token under include/, src/, examples/, tests/"
fi
exit "${status}"
