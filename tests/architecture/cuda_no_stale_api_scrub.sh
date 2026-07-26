#!/usr/bin/env bash
# cuda-resurrect gate: the CUDA backend must not reference the deleted host CPU
# API it rotted against. Pure grep over the cuda trees only — runs on the login
# node (no CUDA toolkit needed), every build configuration.
#
# Three ban families, each an unambiguous marker of the pre-restructure API:
#   1. the deleted catch-all `Real` type and `core/real.h` header (host reals are
#      StorageReal/AccReal in core/precision.h; device reals are real_t/work_t/
#      prec_t in cuda/frame/types.cuh — both keep their qualifiers, so a BARE
#      `Real` word is always the deleted type).
#   2. the renamed LAMMPS unit constants (now core/units.h descriptive names:
#      units::force_to_accel / mass_v2_to_energy / kB_ev_per_k / pi / coulomb).
#   3. stale access on a HOST cpu::Frame. cpu::Frame is now blocks+meta+box,
#      read via .span(k...)/.meta.get(k...)/.box. The device cuda::Frame
#      LEGITIMATELY owns .atoms/.box/.scalars (raw device pointers) — so we ban
#      only the accessor SHAPES that can appear on a host frame and NEVER on the
#      device pointer frame:
#        .atoms.<col>.data()   host span buffer   (device member is a raw ptr)
#        .atoms.size()         host count         (device uses n_atoms())
#        .bond_atom( / .angle_atom( / .dihedral_atom(   host METHOD form
#                              (device topology is an int* indexed with [])
#        .scalars.<field>      host VALUE access   (device is a pointer: ->)
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
DIRS=("${ROOT}/include/atomiverse/cuda" "${ROOT}/src/cuda" "${ROOT}/tests/cuda")
fail=0

report() { # <label> <extended-regex>
    local label="$1" re="$2" hits
    hits="$(grep -rnE "${re}" "${DIRS[@]}" 2>/dev/null)"
    if [ -n "${hits}" ]; then
        echo "FAIL: ${label}"
        echo "${hits}" | head -20
        fail=1
    fi
}

# Same as report(), but drops doc/comment lines first: the bare word "Real"
# is a legit physics term in prose ("Real-space Coulomb", "Real → complex FFT")
# — only in CODE is a bare `Real` the deleted catch-all type.
report_code() { # <label> <extended-regex>
    local label="$1" re="$2" hits
    hits="$(grep -rnE "${re}" "${DIRS[@]}" 2>/dev/null \
            | grep -vE ':[0-9]+:[[:space:]]*(//|\*|/\*)')"
    if [ -n "${hits}" ]; then
        echo "FAIL: ${label}"
        echo "${hits}" | head -20
        fail=1
    fi
}

report      "deleted core/real.h include"   'core/real\.h'
report_code "deleted catch-all Real type"   '\bReal\b'
report "renamed unit constants"             '\b(kFtm2v|kMvv2e|kBoltzmann|kPI|kCoulomb)\b'
report "host cpu::Frame span .data() access" '\.atoms\.[A-Za-z_]+\.data\(\)'
report "host cpu::Frame .atoms.size()"      '\.atoms\.size\(\)'
report "host cpu::Frame topology method"    '\.(bond_atom|angle_atom|dihedral_atom)\('
report "host cpu::Frame .scalars value"     '\.scalars\.[A-Za-z_]'

if [ ${fail} -eq 0 ]; then
    echo "PASS: cuda backend free of deleted host Frame / precision / unit API"
else
    exit 1
fi
