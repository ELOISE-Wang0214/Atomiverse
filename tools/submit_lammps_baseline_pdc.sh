#!/usr/bin/env bash
# tools/submit_lammps_baseline_pdc.sh — submit the lammps_baseline reference
# runs to PDC Dardel `main` partition, one sbatch job per rank count.
#
# Default rank-count list: 1, 32, 64, 128 (one node holds 128 cores on
# Dardel `main` so every job fits on a single node).
#
# Usage: ./tools/submit_lammps_baseline_pdc.sh -A <account> [options]
#
# Account: -A is REQUIRED (no default). PDC CPU project is naiss2025-1-25.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -A <account> [options]

  -A <account>   PDC project account (REQUIRED). CPU: naiss2025-1-25.
  -r <ranks>     Comma-separated rank-count list (default: 1,32,64,128).
  -p <part>      SLURM partition (default: main).
  -t <wall>      Wall-time HH:MM:SS (default: 24:00:00).
  --dry-run      Print job scripts instead of submitting.
  -h, --help     Show this help.

For each rank count R in the list, this submits ONE sbatch job that
runs tools/run_msu_lammps_baseline.sh R on a single allocated node
under srun -n R. Logs land in lammps_baseline/runs/nRr/.
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" || $# -lt 1 ]]; then
  usage
  exit 0
fi

ACCOUNT=""
RANKS_LIST="1,32,64,128"
PARTITION="main"
WALLTIME="24:00:00"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -A) ACCOUNT="${2-}"; shift 2 || { echo "error: -A needs a value" >&2; exit 2; } ;;
    -r) RANKS_LIST="${2-}"; shift 2 || { echo "error: -r needs a value" >&2; exit 2; } ;;
    -p) PARTITION="${2-}"; shift 2 || { echo "error: -p needs a value" >&2; exit 2; } ;;
    -t) WALLTIME="${2-}"; shift 2 || { echo "error: -t needs a value" >&2; exit 2; } ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) echo "error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$ACCOUNT" ]]; then
  echo "error: -A <account> is REQUIRED (no default). CPU project: naiss2025-1-25." >&2
  usage >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_RUNNER="$REPO_ROOT/tools/run_msu_lammps_baseline.sh"
if [[ ! -x "$LOCAL_RUNNER" ]]; then
  echo "error: $LOCAL_RUNNER missing or not executable." >&2
  exit 3
fi

IFS=',' read -r -a RANKS_ARRAY <<< "$RANKS_LIST"
for R in "${RANKS_ARRAY[@]}"; do
  if ! [[ "$R" =~ ^[0-9]+$ ]] || [[ "$R" -lt 1 ]]; then
    echo "error: invalid rank count in list: '$R'" >&2
    exit 2
  fi
done

TS="$(date +%Y%m%d_%H%M%S)"

for R in "${RANKS_ARRAY[@]}"; do
  JOB_NAME="lmp_baseline_${R}r"
  # Dardel `main` nodes have 128 CPU cores; pack each rank count into a
  # single node so MPI traffic stays on-node.
  NODES=1
  TASKS="$R"

  JOB_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/${JOB_NAME}_${TS}.XXXXXX.sbatch")"
  {
    echo "#!/usr/bin/env bash"
    echo "#SBATCH --job-name=${JOB_NAME}"
    echo "#SBATCH --account=${ACCOUNT}"
    echo "#SBATCH --partition=${PARTITION}"
    echo "#SBATCH --time=${WALLTIME}"
    echo "#SBATCH --nodes=${NODES}"
    echo "#SBATCH --ntasks=${TASKS}"
    echo "#SBATCH --output=${REPO_ROOT}/benchmarks/md/reference/polyethylene/lammps_baseline/runs/slurm_${JOB_NAME}_%j.out"
    echo "#SBATCH --error=${REPO_ROOT}/benchmarks/md/reference/polyethylene/lammps_baseline/runs/slurm_${JOB_NAME}_%j.err"
    echo
    echo "set -euo pipefail"
    echo "# Atomiverse lammps_baseline reference run — submitted by"
    echo "# tools/submit_lammps_baseline_pdc.sh on ${TS}."
    echo "# The local runner drives its own 'srun -n ${TASKS}' across this"
    echo "# allocation. We invoke it DIRECTLY -- never wrap it in another srun."
    echo "${LOCAL_RUNNER} ${R}"
  } > "$JOB_SCRIPT"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "── ${JOB_NAME} (dry-run) ──────────────────────────────────────"
    cat "$JOB_SCRIPT"
    echo "────────────────────────────────────────────────────────────────"
    continue
  fi

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "error: sbatch not on PATH; cannot submit. Use --dry-run to inspect." >&2
    exit 4
  fi
  echo "[submit] ${JOB_NAME}  nodes=${NODES} ntasks=${TASKS} time=${WALLTIME} account=${ACCOUNT}"
  sbatch "$JOB_SCRIPT"
done
