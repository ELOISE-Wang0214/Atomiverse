#!/bin/bash
# CP2K leg of the plane-wave SCF performance benchmark: single-point
# bulk Si (a = 5.43 Å) at increasing density-grid cutoffs, appended to
# the same CSV schema bench_pw_scf_perf emits, backend=cp2k.
#
# Usage: ./run_perf.sh /path/to/cp2k.ssmp [csv_out] [workdir]
# Note: CP2K's CUTOFF (density grid, Ry) is not the same quantity as
# the plane-wave wavefunction cutoff in the atomiverse rows — rows are
# comparable as time-to-converged-energy for the same physical system,
# not at equal cutoff values. See ../README.md.
set -euo pipefail

CP2K="${1:?usage: run_perf.sh /path/to/cp2k_binary [csv_out] [workdir]}"
HERE="$(cd "$(dirname "$0")" && pwd)"
CSV="${2:-$HERE/pw_scf_perf.csv}"
WORK="${3:-$HERE/perf_work}"
TEMPLATE="$HERE/../../validation/reference/cp2k/si_bulk.inp.template"
A=5.43
H=$(awk "BEGIN{printf \"%.6f\", $A/2}")
Q=$(awk "BEGIN{printf \"%.6f\", $A/4}")
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"

mkdir -p "$WORK"
if [ ! -f "$CSV" ]; then
  echo "backend,n_atoms,a_angstrom,ecut_ry,k_mesh,threads,iterations,converged,energy_ry,wall_s,s_per_iter" > "$CSV"
fi

for CUTOFF in 200 300 400; do
  INP="$WORK/si_perf_${CUTOFF}.inp"
  OUT="$WORK/si_perf_${CUTOFF}.out"
  sed -e "s/@A@/$A/g" -e "s/@H@/$H/g" -e "s/@Q@/$Q/g" \
      -e "s/@CUTOFF@/$CUTOFF/g" -e "s/@PROJECT@/si_perf_${CUTOFF}/g" \
      "$TEMPLATE" > "$INP"
  echo "=== CP2K CUTOFF = $CUTOFF Ry (OMP_NUM_THREADS=$OMP_NUM_THREADS) ==="
  T0=$(date +%s.%N)
  (cd "$WORK" && "$CP2K" -i "$INP" -o "$OUT")
  T1=$(date +%s.%N)
  WALL=$(awk "BEGIN{printf \"%.3f\", $T1-$T0}")
  E_HA=$(grep "ENERGY| Total FORCE_EVAL" "$OUT" | tail -1 | awk '{print $NF}')
  E_RY=$(awk "BEGIN{printf \"%.8f\", $E_HA*2}")
  ITER=$(grep -c "^ *[0-9]\+ .*OT\|^ *[0-9]\+ .*Diag\." "$OUT" || true)
  CONV=$(grep -q "SCF run converged" "$OUT" && echo 1 || echo 0)
  SPI=$(awk "BEGIN{if ($ITER>0) printf \"%.3f\", $WALL/$ITER; else print 0}")
  echo "cp2k,2,$A,$CUTOFF,4,$OMP_NUM_THREADS,$ITER,$CONV,$E_RY,$WALL,$SPI" >> "$CSV"
  echo "    E = $E_HA Ha, wall = ${WALL}s, scf_steps = $ITER"
done

echo "CSV appended to $CSV"
