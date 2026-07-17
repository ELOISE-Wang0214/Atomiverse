#!/usr/bin/env bash
# Grep gate: no free-function factories; no intermediate-type construction
# in user-facing trees (examples/app).
#
# Public API: concrete methods in methods.h (RHF / UHF / RKS / UKS) + Frame.
# Kernel unit tests may use scf/internal/ — not gated here.
#
# Usage: no_factory_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: no_factory_scrub.sh <repo-root>}"

FACTORY_PATS=(
    'build_integrals[[:space:]]*\('
    'make_atom_list[[:space:]]*\('
    'build_shell_list[[:space:]]*\('
    'build_integral_derivatives[[:space:]]*\('
    'build_am1_system[[:space:]]*\('
    'build_nddo_system[[:space:]]*\('
    'rhf_config[[:space:]]*\('
    'uhf_config[[:space:]]*\('
    'uks_config[[:space:]]*\('
    'default_config[[:space:]]*\('
    'b3lyp_config[[:space:]]*\('
    'pbe0_config[[:space:]]*\('
    'baseline_config[[:space:]]*\('
)

INTERMEDIATE_PATS=(
    'IntegralPack[[:space:]]*\{'
    'AtomList[[:space:]]*\{'
    'ShellList[[:space:]]*\{'
    'IntegralDerivPack[[:space:]]*\{'
    'SCFConfig[[:space:]]*\{'
    'SCFSolver[[:space:]]+[A-Za-z_]'
    'NDDOSystem[[:space:]]+[A-Za-z_]'
)

status=0

for dir in "include/atomiverse" "src" "examples" "tests" "regressions" "app"; do
    target="${ROOT}/${dir}"
    [ -d "${target}" ] || continue
    for pat in "${FACTORY_PATS[@]}"; do
        matches="$(grep -rnE "${pat}" "${target}" 2>/dev/null \
            | grep -v 'no_factory_scrub.sh' \
            | grep -v 'build_integrals_libcint' \
            | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*)' \
            | grep -vE 'Wrong:|banned|Never `' \
            || true)"
        if [ -n "${matches}" ]; then
            echo "FAIL: banned free-function factory /${pat}/ in ${dir}:"
            echo "${matches}"
            echo "  -> use concrete methods (RHF/UHF/RKS/UKS) with Frame."
            status=1
        fi
    done
done

# Intermediate construction: examples + app only (no internal-include gate).
for dir in "examples" "app"; do
    target="${ROOT}/${dir}"
    [ -d "${target}" ] || continue
    for pat in "${INTERMEDIATE_PATS[@]}"; do
        matches="$(grep -rnE "${pat}" "${target}" 2>/dev/null \
            | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*)' \
            || true)"
        if [ -n "${matches}" ]; then
            echo "FAIL: intermediate type in user-facing ${dir}/ /${pat}/:"
            echo "${matches}"
            echo "  -> Frame + RHF/UHF/RKS/UKS only."
            status=1
        fi
    done
done

if [ "${status}" -eq 0 ]; then
    echo "PASS: no banned factories; examples/app free of intermediate types."
fi
exit "${status}"
