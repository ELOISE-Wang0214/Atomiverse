#!/usr/bin/env bash
# tools/audit_lammps_deck.sh — deck-fidelity audit for the MSU-PE Phase D.1
# reference run (.claude/specs/msu-pe-phase-d-01-reference-runs.md, ac-001).
#
# Asserts that the *physics* of benchmarks/md/reference/polyethylene/lammps.in
# is identical to the immutable upstream paper deck lammps_msu_paper.in. The
# two decks differ only in non-physics output hooks (fix print / dump / RDF /
# log / write_restart / reset_timestep / thermo_style ...); every fix,
# integrator, ensemble, timestep, run step-count and force-field coefficient
# must match.
#
# Usage:
#   ./tools/audit_lammps_deck.sh [DECK] [PAPER_DECK]
#
#   DECK        derived deck to audit  (default: the in-tree lammps.in)
#   PAPER_DECK  immutable paper deck   (default: the in-tree lammps_msu_paper.in)
#
# Exit codes:
#   0  — physics lines are set-equal: lammps.in physics ≡ lammps_msu_paper.in.
#   1  — a physics deviation was found; the offending lines are printed.
#   2  — a deck file is missing or a usage error.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [DECK] [PAPER_DECK]

  DECK        derived deck to audit (default: in-tree lammps.in)
  PAPER_DECK  immutable paper deck  (default: in-tree lammps_msu_paper.in)

Exit 0 when the physics lines are set-equal, 1 on a deviation, 2 on a
missing file. Non-physics output hooks (fix print / dump / compute rdf /
log / write_restart / reset_timestep / thermo_style ...) are stripped
before the comparison.
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" ]]; then
  usage
  exit 0
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$REPO_ROOT/benchmarks/md/reference/polyethylene"
DECK="${1:-$REF_DIR/lammps.in}"
PAPER="${2:-$REF_DIR/lammps_msu_paper.in}"

for f in "$DECK" "$PAPER"; do
  if [[ ! -f "$f" ]]; then
    echo "error: deck file not found: $f" >&2
    exit 2
  fi
done

# extract_physics — emit a deck's physics lines, one per line, normalized.
#
# Step 1 — coalesce `&`-continued LAMMPS statements into one logical line
#          (the canonical-CSV `fix print` in lammps.in spans 4 lines; LAMMPS
#          uses a trailing `&` for line continuation, not a backslash).
# Step 2 — drop `#`-comment lines and blank lines.
# Step 3 — drop non-physics output hooks: fix print / fix thermolog /
#          fix ss5 / fix rdfavg / fix def1 / fix def2 / fix ave/time /
#          dump / undump / compute *rdf / variable p* / variable t* /
#          variable pstep|strain|sigxx|lx0|L0|tmp|simname|fname /
#          thermo_style / thermo / thermo_modify / log / write_restart /
#          reset_timestep / bare `run 0` / `print` / and the `unfix` of any
#          of those output-hook fixes (unfix thermolog|ss5|rdfavg|def1|def2).
# Step 4 — normalize: collapse whitespace, trim ends, and reduce
#          `read_data <file>` to a bare `read_data` (the two decks read the
#          byte-identical data file under different names — msu_pe.data is a
#          renamed copy of PE_nc10_cl1000.dat).
extract_physics() {
  local deck="$1"
  awk '
    # --- step 1: coalesce `&` line continuations ---
    {
      cur = $0
      while (sub(/[ \t]*&[ \t]*$/, " ", cur) > 0) {
        if ((getline nxt) <= 0) break
        cur = cur nxt
      }
      $0 = cur
    }

    # --- step 2: strip comments and blank lines ---
    { sub(/[ \t]*#.*$/, "") }
    /^[ \t]*$/ { next }

    # --- step 3: drop non-physics output hooks (either deck) ---
    /^[ \t]*fix[ \t]+[A-Za-z0-9_]+[ \t]+all[ \t]+print[ \t]/ { next }
    /^[ \t]*fix[ \t]+thermolog[ \t]/ { next }
    /^[ \t]*fix[ \t]+ss5[ \t]/       { next }
    /^[ \t]*fix[ \t]+rdfavg[ \t]/    { next }
    /^[ \t]*fix[ \t]+def1[ \t]/      { next }
    /^[ \t]*fix[ \t]+def2[ \t]/      { next }
    /^[ \t]*fix[ \t]+[A-Za-z0-9_]+[ \t]+all[ \t]+ave\/time[ \t]/ { next }
    /^[ \t]*unfix[ \t]+(thermolog|ss5|rdfavg|def1|def2)[ \t]*$/ { next }
    /^[ \t]*dump[ \t]/        { next }
    /^[ \t]*undump[ \t]/      { next }
    /^[ \t]*compute[ \t].*[ \t]rdf[ \t]/ { next }
    /^[ \t]*variable[ \t]+p[0-9A-Za-z_]*[ \t]/ { next }
    /^[ \t]*variable[ \t]+t[0-9A-Za-z_]*[ \t]/ { next }
    /^[ \t]*variable[ \t]+pstep[ \t]/  { next }
    /^[ \t]*variable[ \t]+strain[ \t]/ { next }
    /^[ \t]*variable[ \t]+sigxx[ \t]/  { next }
    /^[ \t]*variable[ \t]+lx0[ \t]/    { next }
    /^[ \t]*variable[ \t]+L0[ \t]/     { next }
    /^[ \t]*variable[ \t]+tmp[ \t]/    { next }
    /^[ \t]*variable[ \t]+simname[ \t]/ { next }
    /^[ \t]*variable[ \t]+fname[ \t]/   { next }
    /^[ \t]*thermo_style[ \t]/  { next }
    /^[ \t]*thermo_modify[ \t]/ { next }
    /^[ \t]*thermo[ \t]/        { next }
    /^[ \t]*log[ \t]/           { next }
    /^[ \t]*write_restart[ \t]/ { next }
    /^[ \t]*reset_timestep[ \t]/ { next }
    /^[ \t]*print[ \t]/         { next }
    # bare `run 0` is a non-physics state-flush (the paper deck issues it
    # before the stage-5 variable block); a real physics run has a count.
    /^[ \t]*run[ \t]+0[ \t]*$/  { next }

    # --- step 4: physics line — normalize whitespace, normalize read_data ---
    {
      line = $0
      gsub(/[ \t]+/, " ", line)
      sub(/^ /, "", line); sub(/ $/, "", line)
      if (line ~ /^read_data /) { line = "read_data" }
      print line
    }
  ' "$deck" | sort
}

export LC_ALL=C

PHYS_DECK="$(extract_physics "$DECK")"
PHYS_PAPER="$(extract_physics "$PAPER")"

ONLY_DECK="$(comm -23 <(printf '%s\n' "$PHYS_DECK") <(printf '%s\n' "$PHYS_PAPER"))"
ONLY_PAPER="$(comm -13 <(printf '%s\n' "$PHYS_DECK") <(printf '%s\n' "$PHYS_PAPER"))"

if [[ -z "$ONLY_DECK" && -z "$ONLY_PAPER" ]]; then
  N_PHYS="$(printf '%s\n' "$PHYS_DECK" | grep -c . || true)"
  echo "[audit_lammps_deck] PASS — physics lines set-equal ($N_PHYS lines)"
  echo "  deck       : $DECK"
  echo "  paper deck : $PAPER"
  exit 0
fi

echo "[audit_lammps_deck] FAIL — physics deviation between decks" >&2
echo "  deck       : $DECK" >&2
echo "  paper deck : $PAPER" >&2
if [[ -n "$ONLY_DECK" ]]; then
  echo "physics lines only in $(basename "$DECK"):" >&2
  printf '  + %s\n' "$ONLY_DECK" | sed '/^  + $/d' >&2
fi
if [[ -n "$ONLY_PAPER" ]]; then
  echo "physics lines only in $(basename "$PAPER"):" >&2
  printf '  - %s\n' "$ONLY_PAPER" | sed '/^  - $/d' >&2
fi
exit 1
