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
Usage: $(basename "$0") [RANKS] [--mode cpu|gpu] [-- <extra lmp args>]

  RANKS         MPI rank count (default: 8). Positional argument.
  --mode cpu    CPU mode (default): no GPU-package args; lmp_cuda runs as
                a plain MPI/OpenMP CPU build.
  --mode gpu    GPU mode (msu-pe-cuda-components-08): adds
                -pk gpu 1 -sf gpu to the lmp command line, wiring the
                GPU-package style with one GPU per rank. Run directory
                is tagged lammps_<RANKS>r_gpu_<TS>/ so the
                canonicalize-msu-output backend-autodetect picks "cuda".
  --            Forward subsequent args verbatim to lmp.

Environment:
  OMP_NUM_THREADS, OMP_PROC_BIND, OMP_PLACES — propagated as-is.
  LMP_BIN   — override LAMMPS binary (default: lmp_cuda).
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" ]]; then
  usage
  exit 0
fi

# Default mode is cpu; --mode {cpu,gpu} may appear before, between, or
# after the RANKS positional (but before `--`).
MODE="cpu"
POSITIONAL=()
while (( $# )); do
  case "$1" in
    --mode)
      shift
      MODE="${1-}"
      if [[ "$MODE" != "cpu" && "$MODE" != "gpu" ]]; then
        echo "error: --mode must be cpu or gpu, got: $MODE" >&2
        usage >&2
        exit 2
      fi
      shift || true
      ;;
    --)
      shift
      break
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done
EXTRA_ARGS=("$@")
RANKS="${POSITIONAL[0]-8}"

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

# The reference LAMMPS binary (lmp_cuda) is dynamically linked against the
# CUDA 12.x runtime (libcudart.so.12) even when run in CPU mode. Source the
# project build-environment file — it puts the CUDA runtime + OpenBLAS on
# LD_LIBRARY_PATH. Guarded: skipped if absent. Override the path with
# ATV_ENV_FILE; point it at an empty file to opt out.
ATV_ENV_FILE="${ATV_ENV_FILE:-$HOME/atv-build-env/env.sh}"
if [[ -f "$ATV_ENV_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$ATV_ENV_FILE"
fi

# Default to the user's single lmp_cuda binary run in CPU mode: the CPU arm
# passes NO GPU package args (no `-sf gpu` / `-pk gpu`), so lmp_cuda runs as
# a plain MPI/OpenMP CPU build. LMP_BIN still overrides the binary.
LMP="${LMP_BIN:-lmp_cuda}"
if ! command -v "$LMP" >/dev/null 2>&1; then
  echo "error: LAMMPS binary '$LMP' not on PATH (override with LMP_BIN=...)." >&2
  exit 4
fi

TS="$(date +%Y%m%d_%H%M%S)"
# In GPU mode, tag the run-dir name with `_gpu_` so
# canonicalize_msu_output.py's --backend auto picks `cuda` from the path
# (per msu-pe-cuda-components-08 § Design).
if [[ "$MODE" == "gpu" ]]; then
  RUN_DIR="$REF_DIR/runs/lammps_${RANKS}r_gpu_${TS}"
else
  RUN_DIR="$REF_DIR/runs/lammps_${RANKS}r_${TS}"
fi
mkdir -p "$RUN_DIR"

cp "$INPUT" "$RUN_DIR/lammps.in"
cp "$DATA"  "$RUN_DIR/msu_pe.data"

{
  echo "RANKS=$RANKS"
  echo "MODE=$MODE"
  echo "LMP_BIN=$LMP"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS-unset}"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND-unset}"
  echo "OMP_PLACES=${OMP_PLACES-unset}"
  echo "lmp_version=$("$LMP" -h 2>&1 | head -2 | tr '\n' ' ')"
  echo "host=$(uname -a)"
  echo "args=${EXTRA_ARGS[*]-(none)}"
} > "$RUN_DIR/run_env.txt"

# msu-pe-cuda-components-08 § Design: GPU-package args are prepended to
# the lmp command line so they apply to every loaded fix / pair style
# that has a /gpu variant. `-sf gpu` is the suffix flag (LAMMPS will try
# the /gpu variant before the base style); `-pk gpu 1` tells the GPU
# package to use one GPU per rank.
MODE_ARGS=()
if [[ "$MODE" == "gpu" ]]; then
  MODE_ARGS=(-pk gpu 1 -sf gpu)
fi

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
( cd "$RUN_DIR" && "${LAUNCH[@]}" "$LMP" "${MODE_ARGS[@]}" -in lammps.in -log log.lammps "${EXTRA_ARGS[@]}" )
END="$(date +%s.%N)"

WALL="$(awk -v s="$START" -v e="$END" 'BEGIN { printf "%.3f", e - s }')"

# Publish the deck's output as committed canonical-schema reference
# artifacts, suffixed by engine and rank count, directly under the
# reference directory so Phase D.2/D.3/D.4 reference them by stable path.
# The deck's `fix thermolog` / `fix ss5` already stream the canonical
# schema; routing through canonicalize_msu_output.py (rather than a bare
# cp) audits the header / column order with the shared, tested converter
# so LAMMPS and Atomiverse use one canonicalization path.
CANON="$REPO_ROOT/tools/canonicalize_msu_output.py"
# msu-pe-cuda-components-08: GPU runs emit backend-tagged canonical
# names so the link-08 compare_msu.py can resolve them by stable path.
# CPU runs keep the legacy untagged names for back-compat with link-01.
if [[ "$MODE" == "gpu" ]]; then
  THERMO_REF="$REF_DIR/thermo_lammps_cuda_${RANKS}r.csv"
  STRESS_REF="$REF_DIR/stress_strain_lammps_cuda_${RANKS}r.csv"
else
  THERMO_REF="$REF_DIR/thermo_lammps_${RANKS}r.csv"
  STRESS_REF="$REF_DIR/stress_strain_lammps_${RANKS}r.csv"
fi
if [[ -s "$RUN_DIR/thermo.lammps.csv" ]]; then
  python3 "$CANON" lammps thermo "$RUN_DIR/thermo.lammps.csv" "$THERMO_REF"
fi
if [[ -s "$RUN_DIR/stress_strain.lammps.csv" ]]; then
  python3 "$CANON" lammps stress-strain \
    "$RUN_DIR/stress_strain.lammps.csv" "$STRESS_REF"
fi

echo
echo "[run_msu_lammps] complete"
echo "  ranks      : $RANKS"
echo "  mode       : $MODE"
echo "  wall time  : ${WALL}s"
echo "  run dir    : $RUN_DIR"
echo "  env        : $RUN_DIR/run_env.txt"
echo "  log        : $RUN_DIR/log.lammps"
echo "  thermo ref : $THERMO_REF"
echo "  stress ref : $STRESS_REF"
