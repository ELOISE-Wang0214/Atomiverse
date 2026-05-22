#!/usr/bin/env bash
# tools/run_msu_pdc.sh — PDC Dardel SLURM submission wrapper for the MSU-PE
# Phase D.1 reference runs (.claude/specs/msu-pe-phase-d-01-reference-runs.md,
# ac-005).
#
# Generates a SLURM batch script that runs one MSU-PE reference engine under
# `srun -n <RANKS>` and submits it with `sbatch`. The job body invokes the
# matching local launcher (tools/run_msu_lammps.sh / tools/run_msu_atomiverse.sh)
# which itself performs the run and the canonicalization.
#
# This wrapper only orchestrates submission — it never reads or writes physics
# artifacts.
#
# Usage:
#   ./tools/run_msu_pdc.sh <engine> <ranks> -A <account> [options]
#
#   engine        lammps | atomiverse                       (required)
#   ranks         MPI rank count, positive integer          (required)
#   -A <account>  PDC project account                       (REQUIRED, no default)
#   -p <part>     SLURM partition          (default: main — the CPU arm)
#   -t <wall>     wall-time HH:MM:SS       (default: 24:00:00 — the PDC cap)
#   -n <nodes>    explicit node count      (default: SLURM picks from -n ranks)
#   --            forward subsequent args verbatim to the run script
#
# PDC Dardel accounts (for reference; the -A flag stays REQUIRED, no hardcoded
# default):
#   naiss2025-1-25       CPU  — partition `main`   (this CPU arm)
#   naiss2025-1-25-gpu   GPU  — partition `gpugh`  (umbrella Phase I, not here)
#
# When `sbatch` is absent (e.g. on a non-SLURM login shell) the generated job
# script is printed to stdout with a clear note and nothing is submitted, so
# the wrapper is testable without a scheduler.
#
# Exit codes:
#   0  — job submitted, or job script printed because sbatch is absent.
#   2  — a usage error (missing engine/ranks, missing -A account, bad value).

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") <engine> <ranks> -A <account> [options]

  engine        lammps | atomiverse                  (required)
  ranks         MPI rank count, positive integer     (required)
  -A <account>  PDC project account                  (REQUIRED, no default)
  -p <part>     SLURM partition   (default: main)
  -t <wall>     wall-time HH:MM:SS (default: 24:00:00)
  -n <nodes>    explicit node count (default: scheduler-chosen)
  --            forward subsequent args verbatim to the run script

PDC Dardel accounts (the -A flag stays required — these are reference only):
  naiss2025-1-25       CPU  — partition main   (this CPU arm)
  naiss2025-1-25-gpu   GPU  — partition gpugh  (umbrella Phase I)

Submits via sbatch a SLURM job whose body runs run_msu_<engine>.sh under
srun -n <ranks>. When sbatch is absent the job script is printed instead.
EOF
}

if [[ "${1-}" == "--help" || "${1-}" == "-h" ]]; then
  usage
  exit 0
fi

# ── positional: engine + ranks ───────────────────────────────────────────────

ENGINE="${1-}"
RANKS="${2-}"
shift 2 2>/dev/null || true

if [[ -z "$ENGINE" || -z "$RANKS" ]]; then
  echo "error: engine and ranks are required positional arguments." >&2
  usage >&2
  exit 2
fi
if [[ "$ENGINE" != "lammps" && "$ENGINE" != "atomiverse" ]]; then
  echo "error: engine must be 'lammps' or 'atomiverse', got: $ENGINE" >&2
  exit 2
fi
if ! [[ "$RANKS" =~ ^[0-9]+$ ]] || [[ "$RANKS" -lt 1 ]]; then
  echo "error: ranks must be a positive integer, got: $RANKS" >&2
  exit 2
fi

# ── options ──────────────────────────────────────────────────────────────────

ACCOUNT=""
PARTITION="main"
WALLTIME="24:00:00"
NODES=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -A) ACCOUNT="${2-}"; shift 2 || { echo "error: -A needs a value" >&2; exit 2; } ;;
    -p) PARTITION="${2-}"; shift 2 || { echo "error: -p needs a value" >&2; exit 2; } ;;
    -t) WALLTIME="${2-}"; shift 2 || { echo "error: -t needs a value" >&2; exit 2; } ;;
    -n) NODES="${2-}"; shift 2 || { echo "error: -n needs a value" >&2; exit 2; } ;;
    --) shift; EXTRA_ARGS=("$@"); break ;;
    *)  echo "error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# The PDC account is mandatory: no fallback, no default — per project
# convention an unspecified account is a hard error, not a guess.
if [[ -z "$ACCOUNT" ]]; then
  echo "error: -A <account> is REQUIRED (no default PDC account)." >&2
  echo "       CPU arm uses naiss2025-1-25 on the 'main' partition." >&2
  usage >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_SCRIPT="$REPO_ROOT/tools/run_msu_${ENGINE}.sh"
if [[ ! -f "$RUN_SCRIPT" ]]; then
  echo "error: run script not found: $RUN_SCRIPT" >&2
  exit 2
fi

TS="$(date +%Y%m%d_%H%M%S)"
JOB_NAME="msu_pe_${ENGINE}_${RANKS}r"
JOB_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/${JOB_NAME}_${TS}.XXXXXX.sbatch")"

# ── generate the SLURM job script ────────────────────────────────────────────

{
  echo "#!/usr/bin/env bash"
  echo "#SBATCH --job-name=${JOB_NAME}"
  echo "#SBATCH --account=${ACCOUNT}"
  echo "#SBATCH --partition=${PARTITION}"
  echo "#SBATCH --time=${WALLTIME}"
  echo "#SBATCH --ntasks=${RANKS}"
  if [[ -n "$NODES" ]]; then
    echo "#SBATCH --nodes=${NODES}"
  fi
  echo "#SBATCH --output=${JOB_NAME}_%j.out"
  echo "#SBATCH --error=${JOB_NAME}_%j.err"
  echo
  echo "set -euo pipefail"
  echo
  echo "# MSU-PE Phase D.1 reference run — generated by tools/run_msu_pdc.sh."
  echo "# The local launcher does the run + canonicalization; this wrapper"
  echo "# only drives srun on the allocated PDC Dardel ${PARTITION} nodes."
  echo "srun -n ${RANKS} ${RUN_SCRIPT} ${RANKS}${EXTRA_ARGS[*]+ -- ${EXTRA_ARGS[*]}}"
} > "$JOB_SCRIPT"

# ── submit, or print when sbatch is absent ───────────────────────────────────

if command -v sbatch >/dev/null 2>&1; then
  echo "[run_msu_pdc] submitting ${JOB_NAME} (${RANKS} ranks, ${PARTITION}, ${WALLTIME})"
  echo "  account    : ${ACCOUNT}"
  echo "  job script : ${JOB_SCRIPT}"
  sbatch "$JOB_SCRIPT"
else
  echo "[run_msu_pdc] sbatch not on PATH — NOT submitting." >&2
  echo "[run_msu_pdc] generated SLURM job script (${JOB_SCRIPT}):" >&2
  echo "----------------------------------------------------------------------"
  cat "$JOB_SCRIPT"
  echo "----------------------------------------------------------------------"
  echo "[run_msu_pdc] run this on a PDC Dardel login node to submit it." >&2
fi
