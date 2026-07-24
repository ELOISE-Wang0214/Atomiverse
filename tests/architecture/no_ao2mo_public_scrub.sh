#!/usr/bin/env bash
# Grep gate (scf-review-02-eri-hoists): the AO→MO transforms live in the
# private internal header, not the public MP2 header, and the spin-orbital
# antisymmetrization helpers exist in exactly ONE place.
#
#   (a) include/atomiverse/cpu/scf/mp2.h must NOT declare
#       ao2mo_ovov / ao2mo_full / mo_hcore (they move to internal/ao2mo.h).
#   (b) include/atomiverse/cpu/scf/internal/ao2mo.h must EXIST and declare all
#       three transforms.
#   (c) src/cpu/scf must contain no local `struct SO` / `so_eri` copies
#       (the ccsd.cpp / cisd.cpp spin-orbital ERI duplicates).
#   (d) src/cpu/scf must DEFINE no local spat_ / spin_ / phys_ / antisym_
#       helpers (the mp3.cpp / mp4.cpp duplicates); calls to the hoisted
#       `antisym_eri` are allowed and do not match the definition pattern.
#
# Regression gate: keeps the ao2mo transforms off the public mp2.h surface and
# the spin-orbital helpers hosted solely in internal/spin_orbital_eri.h.
#
# Usage: no_ao2mo_public_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: no_ao2mo_public_scrub.sh <repo-root>}"

MP2_HDR="${ROOT}/include/atomiverse/cpu/scf/mp2.h"
AO2MO_HDR="${ROOT}/include/atomiverse/cpu/scf/internal/ao2mo.h"
SCF_SRC="${ROOT}/src/cpu/scf"

# Drop comment-only match lines ("path:line:   // ...") so the gate tracks real
# code, mirroring the no_factory_scrub.sh convention.
strip_comments() {
    grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*|/\*|#)' || true
}

status=0

# ── (a) public mp2.h must not re-declare the ao2mo transforms ─────────────
if [ -f "${MP2_HDR}" ]; then
    matches="$(grep -nE 'ao2mo_(ovov|full)|mo_hcore' "${MP2_HDR}" \
        | strip_comments)"
    if [ -n "${matches}" ]; then
        echo "FAIL: public mp2.h still declares AO→MO transforms:"
        echo "${matches}"
        echo "  -> move ao2mo_ovov/ao2mo_full/mo_hcore to internal/ao2mo.h."
        status=1
    fi
fi

# ── (b) internal/ao2mo.h must exist and declare all three transforms ──────
if [ ! -f "${AO2MO_HDR}" ]; then
    echo "FAIL: missing internal header ${AO2MO_HDR}"
    echo "  -> create it with ao2mo_ovov/ao2mo_full/mo_hcore declarations."
    status=1
else
    for sym in ao2mo_ovov ao2mo_full mo_hcore; do
        if ! grep -qE "\<${sym}\>" "${AO2MO_HDR}"; then
            echo "FAIL: internal/ao2mo.h does not declare ${sym}"
            status=1
        fi
    done
fi

# ── (c) no local struct SO / so_eri copies in src/cpu/scf ─────────────────
if [ -d "${SCF_SRC}" ]; then
    matches="$(grep -rnE 'struct SO\>|so_eri' "${SCF_SRC}" 2>/dev/null \
        | strip_comments)"
    if [ -n "${matches}" ]; then
        echo "FAIL: local spin-orbital ERI duplicate (struct SO / so_eri):"
        echo "${matches}"
        echo "  -> use hosted SpinOrbital / antisym_eri from"
        echo "     internal/spin_orbital_eri.h."
        status=1
    fi

    # ── (d) no local spat_/spin_/phys_/antisym_ DEFINITIONS in src/cpu/scf ─
    # Match only definitions ("inline int spat_(" / "inline double antisym_(");
    # the trailing \> boundary excludes calls to the hoisted `antisym_eri`.
    matches="$(grep -rnE 'inline (int|double) +(spat_|spin_|phys_|antisym_)\>' \
        "${SCF_SRC}" 2>/dev/null | strip_comments)"
    if [ -n "${matches}" ]; then
        echo "FAIL: local spin-orbital helper definition (spat_/spin_/phys_/antisym_):"
        echo "${matches}"
        echo "  -> use hosted SpinOrbital / antisym_eri from"
        echo "     internal/spin_orbital_eri.h."
        status=1
    fi
fi

if [ "${status}" -eq 0 ]; then
    echo "PASS: ao2mo hoisted to internal/ao2mo.h; spin-orbital ERI unified."
fi
exit "${status}"
