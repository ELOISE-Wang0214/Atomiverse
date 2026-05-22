#!/usr/bin/env bash
# tools/run_msu_atomiverse.sh — standardized launcher for the Atomiverse
# MSU amorphous-PE tension benchmark.
#
# Usage: ./tools/run_msu_atomiverse.sh [RANKS] [-- <extra polyethylene_tension args>]
#
# Examples:
#   ./tools/run_msu_atomiverse.sh             # 8 ranks, default scale
#   ./tools/run_msu_atomiverse.sh 1
#   ./tools/run_msu_atomiverse.sh 8 -- --scale 0.1 --quiet
#
# Honours OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES from the caller.
# Writes per-component CSV output under
# benchmarks/md/reference/polyethylene/runs/atomiverse_<RANKS>r_<timestamp>/
# and prints both the final per-step wall time and the run directory path.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [RANKS] [-- <extra polyethylene_tension args>]

  RANKS     MPI rank count (default: 8). Positional argument.
  --        Forward subsequent args verbatim to polyethylene_tension.

Environment:
  OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES — propagated as-is.
  ATV_BUILD_DIR  — override build directory (default: ./build).
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" ]]; then
  usage
  exit 0
fi

RANKS="${1-8}"
shift || true
if [[ "${1-}" == "--" ]]; then
  shift
fi
EXTRA_ARGS=("$@")

if ! [[ "$RANKS" =~ ^[0-9]+$ ]]; then
  echo "error: RANKS must be a positive integer, got: $RANKS" >&2
  usage >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ATV_BUILD_DIR:-$REPO_ROOT/build}"
BIN="$BUILD_DIR/benchmarks/polyethylene_tension"
REF_DIR="$REPO_ROOT/benchmarks/md/reference/polyethylene"

if [[ ! -x "$BIN" ]]; then
  echo "error: $BIN not found or not executable." >&2
  echo "       build with: cmake -B $BUILD_DIR -DATV_BUILD_BENCHMARKS=ON" >&2
  echo "                   cmake --build $BUILD_DIR -j" >&2
  exit 3
fi

TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$REPO_ROOT/benchmarks/md/reference/polyethylene/runs/atomiverse_${RANKS}r_${TS}"
mkdir -p "$RUN_DIR"

# Record env / git tag so runs are auditable later.
{
  echo "RANKS=$RANKS"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS-unset}"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND-unset}"
  echo "OMP_PLACES=${OMP_PLACES-unset}"
  echo "git_describe=$(git -C "$REPO_ROOT" describe --always --dirty 2>/dev/null || echo unknown)"
  echo "host=$(uname -a)"
  echo "args=${EXTRA_ARGS[*]-(none)}"
} > "$RUN_DIR/run_env.txt"

# Pick mpirun / srun if available; otherwise serial fallback for RANKS=1.
LAUNCH=()
if [[ "$RANKS" == "1" ]]; then
  if command -v mpirun >/dev/null 2>&1; then
    LAUNCH=(mpirun -n 1)
  fi
else
  if command -v mpirun >/dev/null 2>&1; then
    LAUNCH=(mpirun -n "$RANKS")
  elif command -v srun >/dev/null 2>&1; then
    LAUNCH=(srun -n "$RANKS")
  else
    echo "error: RANKS=$RANKS but neither mpirun nor srun is on PATH." >&2
    exit 4
  fi
fi

START="$(date +%s.%N)"
# Paper-scale defaults: emit a trajectory frame every 10000 steps and a
# strided RDF. Anything in EXTRA_ARGS comes last and overrides these
# (the benchmark parser takes the final occurrence of a repeated flag).
"${LAUNCH[@]}" "$BIN" --out-dir "$RUN_DIR" \
  --traj-every 10000 --rdf 10.5,200,1000 "${EXTRA_ARGS[@]}"
END="$(date +%s.%N)"

WALL="$(awk -v s="$START" -v e="$END" 'BEGIN { printf "%.3f", e - s }')"

# Crude per-step extraction: count lines in thermo.csv (minus header) and
# divide wall time by that count multiplied by --print-every (default 100).
THERMO="$RUN_DIR/thermo.csv"
PER_STEP_MS="n/a"
if [[ -s "$THERMO" ]]; then
  N_ROWS="$(awk 'NR>1' "$THERMO" | wc -l | tr -d ' ')"
  if [[ "$N_ROWS" -gt 0 ]]; then
    # Steps recorded = rows × (last-step - first-step) / (N_ROWS-1) approx.
    # Conservative: use the last "step" column directly.
    LAST_STEP="$(awk -F, 'NR>1 {s=$1} END{print s}' "$THERMO")"
    if [[ -n "$LAST_STEP" && "$LAST_STEP" -gt 0 ]]; then
      PER_STEP_MS="$(awk -v w="$WALL" -v n="$LAST_STEP" 'BEGIN { printf "%.4f", 1000.0 * w / n }')"
    fi
  fi
fi

# Publish canonical-schema reference artifacts (engine + rank suffixed)
# under the reference directory so Phase D.2/D.3/D.4 reference them by a
# stable path. Canonicalization is delegated to the shared, pytest-covered
# tools/canonicalize_msu_output.py so LAMMPS and Atomiverse use one tested
# converter (no inline awk to drift).
#
# Phase G.2 TODO: the benchmark's thermo.csv carries an isotropic pressure
# column only; the canonical schema's Pxx/Pyy/Pzz are filled with that
# isotropic value as a stand-in. The real per-axis (Pxx/Pyy/Pzz) pressure
# path is umbrella-phase msu-pe-validation.md Phase G.2 — until it lands,
# downstream consumers of per-axis Atomiverse stress must treat these three
# columns as an isotropic placeholder. canonicalize_msu_output.py preserves
# this stand-in deliberately; do NOT "fix" it in this Phase-D.1 link.
CANON="$REPO_ROOT/tools/canonicalize_msu_output.py"
THERMO_REF="$REF_DIR/thermo_atomiverse_${RANKS}r.csv"
STRESS_REF="$REF_DIR/stress_strain_atomiverse_${RANKS}r.csv"
if [[ -s "$THERMO" ]]; then
  python3 "$CANON" atomiverse thermo "$THERMO" "$THERMO_REF"
fi
STRESS="$RUN_DIR/stress_strain.csv"
if [[ -s "$STRESS" ]]; then
  python3 "$CANON" atomiverse stress-strain "$STRESS" "$STRESS_REF"
fi

echo
echo "[run_msu_atomiverse] complete"
echo "  ranks       : $RANKS"
echo "  wall time   : ${WALL}s"
echo "  per step    : ${PER_STEP_MS} ms"
echo "  run dir     : $RUN_DIR"
echo "  env         : $RUN_DIR/run_env.txt"
echo "  thermo      : $RUN_DIR/thermo.csv"
echo "  stress/strain: $RUN_DIR/stress_strain.csv"
echo "  thermo ref  : $THERMO_REF"
echo "  stress ref  : $STRESS_REF"
