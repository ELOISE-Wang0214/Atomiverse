#!/usr/bin/env bash
# Isolation gate: the CUDA SCF inner solver must not see the component
# layer (arch-refactor-04-engines). Everything under cuda/scf/ except the
# SCFEvaluator component wrapper is forbidden from including
# component/component.cuh or component/context.cuh, and from naming
# cuda::Context / BuildContext / LaunchContext / DestroyContext — the
# solver consumes only (void* stream, cuda::Frame&) via
# internal::SolverContext. Pure grep, runs on any node (no CUDA needed).
# Mirrors tests/cuda/scf_phase_scrub.sh.
set -u

root="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0

scan() {
    local pattern="$1"
    local label="$2"
    local hits
    hits=$(grep -rn "$pattern" \
        "$root/include/atomiverse/cuda/scf" \
        "$root/src/cuda/scf" 2>/dev/null \
        | grep -v "scf_evaluator\.cuh" \
        | grep -v "scf_evaluator\.cu")
    if [ -n "$hits" ]; then
        echo "FAIL: $label leaked into the CUDA SCF solver layer:"
        echo "$hits"
        fail=1
    fi
}

scan "component/component\.cuh" "component base include"
scan "component/context\.cuh" "component context include"
scan "cuda::Context" "cuda::Context token"
scan "BuildContext" "BuildContext token"
scan "LaunchContext" "LaunchContext token"
scan "DestroyContext" "DestroyContext token"

if [ "$fail" -ne 0 ]; then
    exit 1
fi
echo "OK: CUDA SCF solver layer is component-free (evaluator excepted)."
