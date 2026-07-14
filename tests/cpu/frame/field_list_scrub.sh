#!/usr/bin/env bash
# Grep gate (frame-block-01, ac-001): the registry is THE single field
# enumeration. The functions that used to hand-copy the per-atom field list —
# Atom::set_capacity, Atom::copy_row, comm::pack_exchange,
# comm::unpack_exchange, comm::max_exchange_doubles — must contain NO member
# names; they iterate for_each_column / collect_columns instead. (grow()
# delegates to set_capacity, and the probe-4 grow data-destruction bug was
# exactly a second hand copy of this list falling out of sync.)
#
# Usage: field_list_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: field_list_scrub.sh <repo-root>}"

# Distinctive member tokens: a hand-written field list can't avoid these.
TOKENS='vy|vz|fy|fz|charge|mass|bond_atom|angle_atom|dihedral_atom|nspecial|special'

status=0

# Extract one function body: from the signature line to the first column-0 '}'.
body_of() { # file, signature-regex
    awk -v sig="$2" '
        $0 ~ sig { inside = 1 }
        inside { print }
        inside && /^}/ { exit }
    ' "$1"
}

check() { # file, signature-regex, label
    local hits
    hits=$(body_of "$1" "$2" | sed 's://.*::' | grep -nE "(^|[^a-zA-Z0-9_])(${TOKENS})([^a-zA-Z0-9_]|$)" | grep -v for_each_column || true)
    if [ -n "$hits" ]; then
        echo "FAIL: hand-written field names inside $3 (must derive from the registry):"
        echo "$hits"
        status=1
    fi
}

ATOM="$ROOT/src/cpu/frame/atom.cpp"
PACK="$ROOT/src/cpu/comm/atom_pack.cpp"

check "$ATOM" '^void Atom::set_capacity'        'Atom::set_capacity'
check "$ATOM" '^void Atom::copy_row'            'Atom::copy_row'
check "$PACK" '^int pack_exchange'              'comm::pack_exchange'
check "$PACK" '^int unpack_exchange'            'comm::unpack_exchange'
check "$PACK" '^int max_exchange_doubles'       'comm::max_exchange_doubles'

if [ "$status" -eq 0 ]; then
    echo "PASS: registry functions carry no hand-written field lists"
fi
exit "$status"
