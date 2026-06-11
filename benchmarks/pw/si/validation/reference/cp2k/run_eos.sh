#!/bin/bash
# Generate the CP2K Si EOS reference: lattice-constant scan matching
# bench_pw_si_cp2k_eos (5.25 … 5.55 Å, 2-atom fcc primitive, LDA).
#
# Usage: ./run_eos.sh /path/to/cp2k.ssmp [workdir]
# Output: eos_reference.csv next to this script (a_angstrom,energy [Ha]).
#
# A CP2K binary is NOT bundled. Either load a site module, or download
# the official static build (any recent version works, e.g.
#   https://github.com/cp2k/cp2k/releases — cp2k-<ver>-Linux-gnu-x86_64.ssmp
# ) and pass its path. The deck inlines basis + pseudopotential, so no
# CP2K data directory is needed.
set -euo pipefail

CP2K="${1:?usage: run_eos.sh /path/to/cp2k_binary [workdir]}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="${2:-$HERE/eos_work}"
CUTOFF=400   # density-grid cutoff [Ry]
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"

mkdir -p "$WORK"
CSV="$HERE/eos_reference.csv"
{
  echo "# CP2K $($CP2K --version 2>/dev/null | head -1 | tr ',' ';')"
  echo "# Si fcc 2-atom primitive, DZVP-GTH-q4 / GTH-PADE-q4, PADE LDA,"
  echo "# MP 4x4x4, CUTOFF ${CUTOFF} Ry / REL_CUTOFF 60, EPS_SCF 1e-7"
  echo "# energy_unit: hartree"
  echo "# a_angstrom,energy"
} > "$CSV"

for A in 5.25 5.30 5.35 5.40 5.45 5.50 5.55; do
  H=$(awk "BEGIN{printf \"%.6f\", $A/2}")
  Q=$(awk "BEGIN{printf \"%.6f\", $A/4}")
  INP="$WORK/si_eos_${A}.inp"
  sed -e "s/@A@/$A/g" -e "s/@H@/$H/g" -e "s/@Q@/$Q/g" \
      -e "s/@CUTOFF@/$CUTOFF/g" -e "s/@PROJECT@/si_eos_${A}/g" \
      "$HERE/si_bulk.inp.template" > "$INP"
  echo "=== a = $A Å ==="
  (cd "$WORK" && "$CP2K" -i "$INP" -o "si_eos_${A}.out")
  E=$(grep "ENERGY| Total FORCE_EVAL" "$WORK/si_eos_${A}.out" | tail -1 \
      | awk '{print $NF}')
  echo "$A,$E" >> "$CSV"
  echo "    E = $E Ha"
done

echo "Wrote $CSV"
