#!/usr/bin/env bash
# Grep gate (scf-review-03-posthf-dedup): the two post-HF hoists each live in
# exactly ONE place, and the free `build_imds` factory is gone.
#
#   (1) The engine headers exist and declare their anchor symbols:
#         include/atomiverse/cpu/scf/internal/eom_ccsd_common.h  (make_tau,
#             struct EomImds)
#         include/atomiverse/cpu/scf/internal/ci_engine.h        (ci_element,
#             dense_ci_solve, struct CiSolution)
#   (2) NO free function `build_imds` anywhere under src/ or include/ — the
#       intermediates build is the computing constructor
#       EomImds(const CCSDSolution&) (this clears the 3 pre-existing
#       no_factory_scrub hits eom_ee:70 / eom_ea:66 / eom_ip:67).
#   (3) make_tau + `struct EomImds` are DEFINED only in eom_ccsd_common.{h,cpp};
#       the determinant-CI machinery (choose / enumerate_strings / occ_list /
#       jw_single_sign / diagonal_element / single_element /
#       same_spin_double_element / opposite_spin_double_element / ci_element /
#       dense_ci_solve / make_rdm12 / struct CiSolution) is DEFINED only in
#       ci_engine.{h,cpp}. No local copy survives in fci.cpp / casscf.cpp /
#       eom_{ee,ea,ip}_ccsd.cpp / ccsd.cpp.
#
# Definition vs call: each symbol's DEFINITION/DECLARATION carries its typed
# first parameter (e.g. `ci_element(std::uint64_t`, `make_tau(const`,
# `choose(int`), which no call site spells. That parameter-type anchor is how
# this pure-grep gate tells a definition from a call in every build config.
#
# Usage: posthf_dedup_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: posthf_dedup_scrub.sh <repo-root>}"

EOM_HDR="${ROOT}/include/atomiverse/cpu/scf/internal/eom_ccsd_common.h"
EOM_SRC="${ROOT}/src/cpu/scf/internal/eom_ccsd_common.cpp"
CI_HDR="${ROOT}/include/atomiverse/cpu/scf/internal/ci_engine.h"
CI_SRC="${ROOT}/src/cpu/scf/internal/ci_engine.cpp"
SCF_SRC="${ROOT}/src/cpu/scf"
SCF_INC="${ROOT}/include/atomiverse/cpu/scf"

# Drop comment-only match lines ("path:line:   // ...") so the gate tracks real
# code, mirroring the no_ao2mo_public_scrub.sh / no_factory_scrub.sh convention.
strip_comments() {
    grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*|#)' || true
}

# Drop the two sanctioned engine files: their definitions are the home.
drop_engine_files() {
    grep -vE '(eom_ccsd_common|ci_engine)\.(h|cpp):' || true
}

status=0

# ── (1) engine headers exist and declare their anchor symbols ─────────────
if [ ! -f "${EOM_HDR}" ]; then
    echo "FAIL: missing ${EOM_HDR}"
    echo "  -> declare make_tau + struct EomImds there; impl in ${EOM_SRC}."
    status=1
else
    for sym in 'make_tau' 'struct EomImds'; do
        if ! grep -qE "${sym}" "${EOM_HDR}"; then
            echo "FAIL: eom_ccsd_common.h does not declare '${sym}'"
            status=1
        fi
    done
fi
if [ ! -f "${CI_HDR}" ]; then
    echo "FAIL: missing ${CI_HDR}"
    echo "  -> declare the CI engine there; impl in ${CI_SRC}."
    status=1
else
    for sym in 'ci_element' 'dense_ci_solve' 'struct CiSolution'; do
        if ! grep -qE "${sym}" "${CI_HDR}"; then
            echo "FAIL: ci_engine.h does not declare '${sym}'"
            status=1
        fi
    done
fi

# ── (2) no free build_imds anywhere under src/ or include/ ────────────────
matches="$(grep -rnE '\<build_imds' "${ROOT}/src" "${ROOT}/include" 2>/dev/null \
    | strip_comments)"
if [ -n "${matches}" ]; then
    echo "FAIL: free function build_imds still present (must become the"
    echo "      computing constructor EomImds(const CCSDSolution&)):"
    echo "${matches}"
    status=1
fi

# ── (3a) make_tau / struct EomImds defined only in eom_ccsd_common.* ──────
matches="$(grep -rnE '\<make_tau[[:space:]]*\(const|struct EomImds' \
    "${SCF_SRC}" "${SCF_INC}" 2>/dev/null | drop_engine_files | strip_comments)"
if [ -n "${matches}" ]; then
    echo "FAIL: make_tau / struct EomImds defined outside eom_ccsd_common:"
    echo "${matches}"
    echo "  -> hoist into internal/eom_ccsd_common.{h,cpp}."
    status=1
fi

# ── (3b) CI-engine symbols defined only in ci_engine.* ────────────────────
CI_DEF_PAT='\<choose[[:space:]]*\(int'
CI_DEF_PAT="${CI_DEF_PAT}|\<enumerate_strings[[:space:]]*\(int"
CI_DEF_PAT="${CI_DEF_PAT}|\<occ_list[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<jw_single_sign[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<diagonal_element[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<single_element[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<same_spin_double_element[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<opposite_spin_double_element[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<ci_element[[:space:]]*\(std::uint64_t"
CI_DEF_PAT="${CI_DEF_PAT}|\<dense_ci_solve[[:space:]]*\(const"
CI_DEF_PAT="${CI_DEF_PAT}|\<make_rdm12[[:space:]]*\(const"
CI_DEF_PAT="${CI_DEF_PAT}|struct CiSolution"
matches="$(grep -rnE "${CI_DEF_PAT}" "${SCF_SRC}" "${SCF_INC}" 2>/dev/null \
    | drop_engine_files | strip_comments)"
if [ -n "${matches}" ]; then
    echo "FAIL: determinant-CI symbol defined outside ci_engine (fci/casscf fork):"
    echo "${matches}"
    echo "  -> hoist into internal/ci_engine.{h,cpp}."
    status=1
fi

if [ "${status}" -eq 0 ]; then
    echo "PASS: post-HF EOM intermediates + determinant-CI engine each hosted once."
fi
exit "${status}"
