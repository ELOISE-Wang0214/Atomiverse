#!/usr/bin/env bash
# Grep gate (frame-block-02, ac-001): src/cpu/ and include/atomiverse/cpu/
# must not access per-atom schema fields through the direct `atoms.<field>[`
# / `atoms.<field>(` column syntax anymore. Every consumer resolves a typed
# `Frame::span(kAtom_*)` ONCE at the launch boundary and loops over the raw
# pointer instead (the 13.7x rule -- cpu/frame/frame.h, core/frame_schema.h).
#
# Field set: the per-atom schema columns that carry a plain member-style
# accessor on cpu::Atom -- x y z vx vy vz fx fy fz mass charge type tag
# image. Topology arrays (num_bond/bond_atom/angle_atom/.../special) and
# bookkeeping (nlocal/nghost/nmax/map) are NOT migrated by this step -- see
# the frame-block-02 spec -- so they are deliberately absent from the list.
#
# Exemptions (each is the registry's OWN binding, or a boundary this step
# deliberately leaves alone -- every one is a conscious call, not an
# oversight):
#   - src/cpu/frame/atom.cpp, cpu/frame/atom.h  -- the member declarations
#     and the for_each_column registry binding THEMSELVES; this is the one
#     place `atoms.x` etc. legitimately name the backing storage.
#   - src/cpu/comm/atom_pack.cpp                -- already registry-driven
#     (iterates for_each_column / ColumnRef; field_list_scrub.sh separately
#     gates its two entry points against hand-written field names).
#   - src/cpu/frame/materialize.cpp             -- the molrs FFI boundary: a
#     structural block <-> member conversion that runs once per `run()`,
#     never in a step loop.
#   - include/atomiverse/cpu/frame/frame.h      -- Frame::span/has/insert ARE
#     the mapping facade; they resolve through Atom::for_each_column, never
#     hand-write a member name.
#   - cuda/ and torch/ subtrees                 -- out of scope for this
#     step; their own gate needs a GPU / LibTorch build to even run. See
#     .claude/specs/frame-block-02-cpu-consumers.md, "Out of scope".
#   - comment lines (`//`, `///`, or a `*` continuation line) -- several
#     docstrings quote the banned pattern to explain this very rule.
#
# Usage: direct_access_scrub.sh <repo-root>
set -u

ROOT="${1:?usage: direct_access_scrub.sh <repo-root>}"

# Distinctive schema-field tokens, immediately after `atoms.` and immediately
# before a column-access `[` or `(` -- the direct-access shape this gate bans.
FIELDS='x|y|z|vx|vy|vz|fx|fy|fz|mass|charge|type|tag|image'
PATTERN="atoms\.(${FIELDS})[[(]"
# Same pattern, doubled backslash: awk's `-v` assignment runs its own
# string-escape pass before the dynamic regex is compiled, so a single
# `\.` would be silently un-escaped (gawk warns "escape sequence treated as
# plain '.'") and the regex would loosely match ANY character there instead
# of a literal dot. `\\.` survives that pass as the literal `\.` we want.
AWK_PATTERN="atoms\\\\.(${FIELDS})[[(]"

status=0
for dir in "src/cpu" "include/atomiverse/cpu"; do
    target="${ROOT}/${dir}"
    [ -d "${target}" ] || continue

    # Two comment filters, since a hit can be a comment TWO different ways:
    #   1. the whole line is a comment (`///`, `//`, or a `*` continuation) --
    #      dropped outright.
    #   2. the match only lives in a TRAILING `//`/`///<` comment after real
    #      code on the same line (e.g. `int type_i; ///< ... atoms.type[i] ...`)
    #      -- re-check the pattern against the code portion alone (everything
    #      before the first `//`) and drop the hit if it no longer matches.
    matches="$(grep -rnE "${PATTERN}" "${target}" 2>/dev/null \
        | grep -vE '^[^:]+:[0-9]+:[[:space:]]*(///|//|\*)' \
        | grep -vE '/(cuda|torch)/' \
        | grep -vE '(^|/)src/cpu/frame/atom\.cpp:' \
        | grep -vE '(^|/)src/cpu/comm/atom_pack\.cpp:' \
        | grep -vE '(^|/)src/cpu/frame/materialize\.cpp:' \
        | grep -vE '(^|/)include/atomiverse/cpu/frame/frame\.h:' \
        | grep -vE '(^|/)include/atomiverse/cpu/frame/atom\.h:' \
        | awk -F: -v pat="${AWK_PATTERN}" '
            {
                content = $0
                sub(/^[^:]*:[^:]*:/, "", content)
                code = content
                sub(/\/\/.*/, "", code)
                if (code ~ pat) print $0
            }' || true)"

    if [ -n "${matches}" ]; then
        echo "FAIL: direct atoms.<field>[/( access remains in ${dir}:"
        echo "${matches}"
        echo "  -> resolve via Frame::span(kAtom_*) once at the launch boundary"
        echo "     instead (cpu/frame/frame.h, core/frame_schema.h)."
        status=1
    fi
done

if [ "${status}" -eq 0 ]; then
    echo "PASS: no direct atoms.<field>[/( access outside the registry mapping facade."
fi
exit "${status}"
