#!/usr/bin/env bash
# tools/run_msu_lammps.sh — standardized launcher for the LAMMPS MSU
# amorphous-PE tension reference run.
#
# Usage: ./tools/run_msu_lammps.sh [RANKS] [-- <extra lmp args>]
#
# Examples:
#   ./tools/run_msu_lammps.sh           # 8 ranks
#   ./tools/run_msu_lammps.sh 1
#
# Honours OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES from the caller.
# Writes log.lammps + reference CSVs into a timestamped subdirectory of
# benchmarks/md/reference/polyethylene/runs/.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [RANKS] [-- <extra lmp args>]

  RANKS     MPI rank count (default: 8). Positional argument.
  --        Forward subsequent args verbatim to lmp.

Environment:
  OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES — propagated as-is.
  LMP_BIN   — override LAMMPS binary (default: lmp on PATH).
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
REF_DIR="$REPO_ROOT/benchmarks/md/reference/polyethylene"
INPUT="$REF_DIR/lammps.in"
DATA="$REF_DIR/msu_pe.data"

if [[ ! -f "$INPUT" ]]; then
  echo "error: $INPUT not found." >&2
  exit 3
fi
if [[ ! -f "$DATA" ]]; then
  echo "error: $DATA not present. See $REF_DIR/README.md for generation steps." >&2
  exit 3
fi

LMP="${LMP_BIN:-lmp}"
if ! command -v "$LMP" >/dev/null 2>&1; then
  echo "error: LAMMPS binary '$LMP' not on PATH (override with LMP_BIN=...)." >&2
  exit 4
fi

TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="$REF_DIR/runs/lammps_${RANKS}r_${TS}"
mkdir -p "$RUN_DIR"

cp "$INPUT" "$RUN_DIR/lammps.in"
cp "$DATA"  "$RUN_DIR/msu_pe.data"

{
  echo "RANKS=$RANKS"
  echo "LMP_BIN=$LMP"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS-unset}"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND-unset}"
  echo "OMP_PLACES=${OMP_PLACES-unset}"
  echo "lmp_version=$("$LMP" -h 2>&1 | head -2 | tr '\n' ' ')"
  echo "host=$(uname -a)"
  echo "args=${EXTRA_ARGS[*]-(none)}"
} > "$RUN_DIR/run_env.txt"

LAUNCH=()
if [[ "$RANKS" == "1" ]]; then
  if command -v mpirun >/dev/null 2>&1; then LAUNCH=(mpirun -n 1); fi
else
  if   command -v mpirun >/dev/null 2>&1; then LAUNCH=(mpirun -n "$RANKS")
  elif command -v srun   >/dev/null 2>&1; then LAUNCH=(srun -n "$RANKS")
  else
    echo "error: RANKS=$RANKS but neither mpirun nor srun is on PATH." >&2
    exit 4
  fi
fi

START="$(date +%s.%N)"
( cd "$RUN_DIR" && "${LAUNCH[@]}" "$LMP" -in lammps.in -log log.lammps "${EXTRA_ARGS[@]}" )
END="$(date +%s.%N)"

WALL="$(awk -v s="$START" -v e="$END" 'BEGIN { printf "%.3f", e - s }')"

echo
echo "[run_msu_lammps] complete"
echo "  ranks      : $RANKS"
echo "  wall time  : ${WALL}s"
echo "  run dir    : $RUN_DIR"
echo "  env        : $RUN_DIR/run_env.txt"
echo "  log        : $RUN_DIR/log.lammps"
echo "  thermo     : $RUN_DIR/thermo_reference.csv"
echo "  stress/str : $RUN_DIR/stress_strain_reference.csv"
