#!/usr/bin/env bash
# Gate: CUDA classical MD force path is OpenMM cluster/tile only — no dual
# NonbondedStyle architecture and no default Verlet (vl_*) force source.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0

report() {
    local label="$1" re="$2"
    local hits
    hits="$(grep -rnE "${re}" \
        "${ROOT}/include/atomiverse/cuda" \
        "${ROOT}/src/cuda" \
        2>/dev/null || true)"
    if [ -n "${hits}" ]; then
        echo "FAIL: ${label}"
        echo "${hits}" | head -30
        fail=1
    fi
}

# Dual-path architecture switches
report "NonbondedStyle dual path" \
    'NonbondedStyle|ATV_CUDA_NONBONDED|nonbonded_style'

# Force styles must not call Verlet launcher (definition may remain quarantined)
hits="$(grep -rnE 'launch_pair_compute\s*\(' \
    "${ROOT}/src/cuda/potentials" \
    "${ROOT}/include/atomiverse/cuda/potentials" \
    2>/dev/null | grep -v 'pair_compute.cuh' || true)"
if [ -n "${hits}" ]; then
    echo "FAIL: MD pair force uses launch_pair_compute (Verlet)"
    echo "${hits}" | head -20
    fail=1
fi

# PairList must not reintroduce public vl_* force fields
hits="$(grep -rnE '\bvl_nbr\b|\bvl_cnt\b' \
    "${ROOT}/include/atomiverse/cuda/neighbor" 2>/dev/null || true)"
if [ -n "${hits}" ]; then
    echo "FAIL: PairList.vl_* field"
    echo "${hits}" | head -20
    fail=1
fi

# Deleted direct Verlet rebuild entry
report "launch_build_verlet_direct still present" \
    'launch_build_verlet_direct'

# Production mask all-ones ban (neighbor rebuild)
if grep -rnE 'cudaMemset(Async)?\([^)]*interaction_mask[^)]*0xFF' \
    "${ROOT}/src/cuda/neighbor" 2>/dev/null | grep -v '^\s*//'; then
    echo "FAIL: interaction_mask memset 0xFF in production neighbor"
    fail=1
fi

if [ "$fail" -ne 0 ]; then
    echo "cuda_openmm_md_scrub: FAILED"
    exit 1
fi
echo "cuda_openmm_md_scrub: OK"
exit 0
