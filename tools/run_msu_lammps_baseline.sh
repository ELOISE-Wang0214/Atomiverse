#!/usr/bin/env bash
# tools/run_msu_lammps_baseline.sh — single-rank-count baseline launcher.
#
# Runs the lammps_baseline/deck.in (physics-verbatim MSU PE protocol with
# extended thermo on every stage) at one MPI rank count, then converts the
# log.lammps into thermo.csv + stage_summary.csv via the molpy log reader.
#
# Usage: ./tools/run_msu_lammps_baseline.sh <RANKS>
#
# Writes everything into:
#   benchmarks/md/reference/polyethylene/lammps_baseline/runs/n<RANKS>r/
#
# Environment:
#   LMP_BIN    override LAMMPS binary (default: lmp from PATH)
#   MOLPY_VENV salamander venv with molpy + molrs
#              (default: ~/salamander_2025/250501/.venv)

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") <RANKS>

  RANKS   MPI rank count (positive integer)

Environment:
  LMP_BIN     override LAMMPS binary (default: 'lmp')
  MOLPY_VENV  override molpy venv path
              (default: ~/salamander_2025/250501/.venv)
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" || $# -lt 1 ]]; then
  usage
  exit 0
fi

RANKS="$1"
if ! [[ "$RANKS" =~ ^[0-9]+$ ]] || [[ "$RANKS" -lt 1 ]]; then
  echo "error: RANKS must be a positive integer, got: $RANKS" >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASELINE_DIR="$REPO_ROOT/benchmarks/md/reference/polyethylene/lammps_baseline"
DECK="$BASELINE_DIR/deck.in"
DATA="$BASELINE_DIR/msu_pe.data"

if [[ ! -f "$DECK" ]]; then
  echo "error: $DECK not found." >&2; exit 3
fi
if [[ ! -e "$DATA" ]]; then
  echo "error: $DATA not found (expected symlink to ../msu_pe.data)." >&2; exit 3
fi

LMP="${LMP_BIN:-lmp}"
if ! command -v "$LMP" >/dev/null 2>&1; then
  echo "error: LAMMPS binary '$LMP' not on PATH (override with LMP_BIN=...)." >&2
  exit 4
fi

RUN_DIR="$BASELINE_DIR/runs/n${RANKS}r"
if [[ -d "$RUN_DIR" ]]; then
  echo "error: $RUN_DIR already exists -- remove it or rename before re-running." >&2
  exit 5
fi
mkdir -p "$RUN_DIR"

# Stage the deck + data into the run dir so the run is self-contained
# and re-runnable from RUN_DIR alone.
cp "$DECK" "$RUN_DIR/deck.in"
cp -L "$DATA" "$RUN_DIR/msu_pe.data"

{
  echo "RANKS=$RANKS"
  echo "LMP_BIN=$LMP"
  echo "lmp_version=$("$LMP" -h 2>&1 | head -2 | tr '\n' ' ')"
  echo "host=$(uname -a)"
  echo "started=$(date -Iseconds)"
} > "$RUN_DIR/run_env.txt"

LAUNCH=()
if   command -v srun   >/dev/null 2>&1; then LAUNCH=(srun -n "$RANKS")
elif command -v mpirun >/dev/null 2>&1; then LAUNCH=(mpirun -n "$RANKS")
else
  echo "error: neither srun nor mpirun on PATH; cannot launch ${RANKS}-rank run." >&2
  exit 4
fi

START="$(date +%s.%N)"
( cd "$RUN_DIR" && "${LAUNCH[@]}" "$LMP" -in deck.in -log log.lammps )
END="$(date +%s.%N)"
WALL="$(awk -v s="$START" -v e="$END" 'BEGIN { printf "%.3f", e - s }')"
echo "ended=$(date -Iseconds)" >> "$RUN_DIR/run_env.txt"
echo "wall_seconds=${WALL}"    >> "$RUN_DIR/run_env.txt"

# ── log → CSV ─────────────────────────────────────────────────────────
MOLPY_VENV="${MOLPY_VENV:-$HOME/salamander_2025/250501/.venv}"
if [[ -f "$MOLPY_VENV/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$MOLPY_VENV/bin/activate"
else
  echo "[run_msu_lammps_baseline] warning: MOLPY_VENV $MOLPY_VENV missing; trying system python." >&2
fi
python3 "$REPO_ROOT/tools/lammps_baseline_log_to_csv.py" --run-dir "$RUN_DIR"

echo
echo "[run_msu_lammps_baseline] complete"
echo "  ranks      : $RANKS"
echo "  wall time  : ${WALL}s"
echo "  run dir    : $RUN_DIR"
echo "  thermo csv : $RUN_DIR/thermo.csv"
echo "  summary csv: $RUN_DIR/stage_summary.csv"
echo "  log        : $RUN_DIR/log.lammps"
