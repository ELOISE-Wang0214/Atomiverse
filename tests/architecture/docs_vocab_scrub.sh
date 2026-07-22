#!/usr/bin/env bash
# backend-mirror-07-docs gate: the Architecture prose in CLAUDE.md and AGENTS.md
# must not resurrect the pre-restructure vocabulary. Bans the class/path names
# that no longer exist; the bare words "backend", "Backend", "kspace" stay legit.
set -u
ROOT="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
fail=0
# Build the ban list without embedding the literal tokens contiguously, so this
# script never matches itself. Paths use a trailing slash; class names are exact.
bans=(
  "fa""cade"
  "CPU""Driver"
  "CUDA""Driver"
  "Loop""Handler"
  "add_""handler"
  "cuda/""kspace"
  "cpu::""Backend"
  "cpu/""md/" "cpu/""comm/" "cpu/""async/" "cpu/""pw/" "cpu/""mc/"
  "cpu/""minimizer/" "cpu/""vib/" "cpu/""observer/" "cpu/""recorder/"
)
for f in "${ROOT}/CLAUDE.md" "${ROOT}/AGENTS.md"; do
    for b in "${bans[@]}"; do
        if grep -qF "$b" "$f" 2>/dev/null; then
            echo "FAIL: $(basename "$f") still contains legacy token '$b'"; fail=1
        fi
    done
done
if [ ${fail} -eq 0 ]; then echo "PASS: CLAUDE/AGENTS free of legacy layout vocabulary"; else exit 1; fi
