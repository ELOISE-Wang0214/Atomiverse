#!/usr/bin/env bash
# tools/submit_atomiverse_baseline_pdc.sh — submit the Atomiverse reference
# runs for the MSU amorphous-PE protocol to PDC Dardel `main` partition,
# mirroring `submit_lammps_baseline_pdc.sh` (one sbatch per rank count).
#
# A single x86 build is produced inside a dedicated "build" sbatch job;
# the four benchmark jobs are submitted with `--dependency=afterok:<build_id>`
# so they only fire after a successful build. Each benchmark job writes to
#   benchmarks/md/reference/polyethylene/atomiverse_baseline/runs/n<R>r/
# parallel to the LAMMPS baseline layout, so the comparison drivers can
# point at matching paths.
#
# Default rank-count list: 1, 32, 64, 128 (one Dardel main node = 128 cores).
#
# Usage: ./tools/submit_atomiverse_baseline_pdc.sh -A <account> [options]
#
# Account: -A is REQUIRED. PDC CPU project is naiss2025-1-25.

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -A <account> [options]

  -A <account>   PDC project account (REQUIRED). CPU: naiss2025-1-25.
  -r <ranks>     Comma-separated rank-count list (default: 1,32,64,128).
  -p <part>      SLURM partition (default: main).
  -t <wall>      Wall-time HH:MM:SS for benchmark jobs (default: 24:00:00).
  --build-time   Wall-time HH:MM:SS for the build job (default: 01:00:00).
  --build-dir    Out-of-tree build directory name (default: build-main).
  --skip-build   Reuse an existing --build-dir; do NOT submit a build job.
  --dry-run      Print job scripts instead of submitting.
  -h, --help     Show this help.

The build job runs cmake + make on an allocated main node (x86, 128 cores)
producing <repo>/<build-dir>/benchmarks/polyethylene_tension. Each subsequent
benchmark job allocates 1 node + <R> ntasks, runs polyethylene_tension under
srun -n <R> with --out-dir set to atomiverse_baseline/runs/n<R>r/, then
canonicalizes thermo.csv via tools/canonicalize_msu_output.py.
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
BUILD_TIME="01:00:00"
BUILD_DIR="build-main"
SKIP_BUILD=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -A) ACCOUNT="${2-}"; shift 2 || { echo "error: -A needs a value" >&2; exit 2; } ;;
    -r) RANKS_LIST="${2-}"; shift 2 || { echo "error: -r needs a value" >&2; exit 2; } ;;
    -p) PARTITION="${2-}"; shift 2 || { echo "error: -p needs a value" >&2; exit 2; } ;;
    -t) WALLTIME="${2-}"; shift 2 || { echo "error: -t needs a value" >&2; exit 2; } ;;
    --build-time) BUILD_TIME="${2-}"; shift 2 || { echo "error: --build-time needs a value" >&2; exit 2; } ;;
    --build-dir)  BUILD_DIR="${2-}";  shift 2 || { echo "error: --build-dir needs a value"  >&2; exit 2; } ;;
    --skip-build) SKIP_BUILD=1; shift ;;
    --dry-run)    DRY_RUN=1; shift ;;
    *) echo "error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$ACCOUNT" ]]; then
  echo "error: -A <account> is REQUIRED (no default). CPU project: naiss2025-1-25." >&2
  usage >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNS_PARENT="$REPO_ROOT/benchmarks/md/reference/polyethylene/atomiverse_baseline/runs"
DATA_FILE="$REPO_ROOT/benchmarks/md/reference/polyethylene/msu_pe.data"
CANON_TOOL="$REPO_ROOT/tools/canonicalize_msu_output.py"
BIN_REL="$BUILD_DIR/benchmarks/polyethylene_tension"
BIN_PATH="$REPO_ROOT/$BIN_REL"
MOLPY_VENV_DEFAULT="$HOME/salamander_2025/250501/.venv"

if [[ ! -f "$DATA_FILE" ]]; then
  echo "error: data file missing: $DATA_FILE" >&2; exit 3
fi

IFS=',' read -r -a RANKS_ARRAY <<< "$RANKS_LIST"
for R in "${RANKS_ARRAY[@]}"; do
  if ! [[ "$R" =~ ^[0-9]+$ ]] || [[ "$R" -lt 1 ]]; then
    echo "error: invalid rank count in list: '$R'" >&2
    exit 2
  fi
done

mkdir -p "$RUNS_PARENT"

TS="$(date +%Y%m%d_%H%M%S)"

# ── env preamble used by build + run jobs ───────────────────────────────────
# Kept inline in every emitted sbatch so the job is self-describing on disk
# and a re-submission survives changes to the orchestrator. The two paths
# baked in are:
#   - /opt/cray/pe/gcc-native/13/bin/g++  (host C++ compiler — same as the
#     aarch64 GH200 build, identical major version on x86 Dardel main)
#   - $HOME/atv-build-env/blas/lib/libopenblas.so  (x86 OpenBLAS symlink to
#     /usr/lib64/libopenblas.so.0; arch-portable since libopenblas-devel is
#     installed on both Dardel arms)
read -r -d '' COMMON_ENV <<'EOF' || true
set -euo pipefail

# Module setup for x86 main partition. The Cray default environment is
# PrgEnv-cray + cray-mpich + cray-libsci. CMake's find_package(MPI) on
# Cray needs the cc/CC compiler wrappers (NOT mpicc/mpicxx — those are
# blocked under any PrgEnv with "CrayPE is loaded, use cc/CC/ftn").
# Swap to PrgEnv-gnu so CC delegates to gcc-13 + auto-MPI (verified
# 2026-05-25 via build-main-smoketest configure).
module unload PrgEnv-cray 2>/dev/null || true
module load PrgEnv-gnu
module load cmake/4.1.2 rust/1.93.0

# Cargo cache lives on scratch — the home fileset is inode-limited and
# the cargo registry/target trees blow past 100k inodes. The user's own
# env.sh routes here for the aarch64 build; same paths reused for x86.
ATV_RUST_CACHE="/cfs/klemming/scratch/j/jichenli/salamander_2025/rust_cache"
export RUSTUP_HOME="$ATV_RUST_CACHE/rustup-home"
export CARGO_HOME="$ATV_RUST_CACHE/cargo-home"
export CARGO_TARGET_DIR="$ATV_RUST_CACHE/cargo-target"
export PATH="$CARGO_HOME/bin:$PATH"

# Force CMake to use the Cray cc/CC wrappers so MPI auto-links.
export CC=cc
export CXX=CC
EOF

# ── 1. build job ────────────────────────────────────────────────────────────

BUILD_JOB_ID=""
if [[ "$SKIP_BUILD" -eq 0 ]]; then
  BUILD_JOB_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/atv_build_${TS}.XXXXXX.sbatch")"
  {
    echo "#!/usr/bin/env bash"
    echo "#SBATCH --job-name=atv_build_main"
    echo "#SBATCH --account=${ACCOUNT}"
    echo "#SBATCH --partition=${PARTITION}"
    echo "#SBATCH --time=${BUILD_TIME}"
    echo "#SBATCH --nodes=1"
    echo "#SBATCH --ntasks=1"
    echo "#SBATCH --cpus-per-task=128"
    echo "#SBATCH --output=${RUNS_PARENT}/slurm_atv_build_main_%j.out"
    echo "#SBATCH --error=${RUNS_PARENT}/slurm_atv_build_main_%j.err"
    echo
    echo "${COMMON_ENV}"
    echo
    echo "cd '${REPO_ROOT}'"
    echo "BLAS_LIB='${HOME}/atv-build-env/blas/lib/libopenblas.so'"
    echo "if [[ ! -e \"\$BLAS_LIB\" ]]; then"
    echo "  echo 'error: BLAS lib missing at '\"\$BLAS_LIB\" >&2"
    echo "  exit 4"
    echo "fi"
    echo
    echo "cmake -B '${BUILD_DIR}' \\"
    echo "  -DCMAKE_BUILD_TYPE=Release \\"
    echo "  -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \\"
    echo "  -DATV_BUILD_CUDA=OFF \\"
    echo "  -DATV_BUILD_RHF=OFF \\"
    echo "  -DATV_BUILD_TESTS=OFF \\"
    echo "  -DATV_BUILD_BENCHMARKS=ON \\"
    echo "  -DATV_USE_MPI=ON \\"
    echo "  -DBLAS_LIBRARIES=\"\$BLAS_LIB\""
    echo
    echo "cmake --build '${BUILD_DIR}' -j 128"
    echo
    echo "test -x '${BIN_REL}' || { echo 'error: ${BIN_REL} not built' >&2; exit 5; }"
    echo "file '${BIN_REL}'"
  } > "$BUILD_JOB_SCRIPT"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "── atv_build_main (dry-run) ────────────────────────────────────"
    cat "$BUILD_JOB_SCRIPT"
    echo "──────────────────────────────────────────────────────────────────"
  else
    if ! command -v sbatch >/dev/null 2>&1; then
      echo "error: sbatch not on PATH; cannot submit. Use --dry-run to inspect." >&2
      exit 4
    fi
    echo "[submit] atv_build_main  partition=${PARTITION} time=${BUILD_TIME} account=${ACCOUNT}"
    BUILD_SUBMIT_OUT="$(sbatch "$BUILD_JOB_SCRIPT")"
    echo "  $BUILD_SUBMIT_OUT"
    BUILD_JOB_ID="$(awk '{print $NF}' <<<"$BUILD_SUBMIT_OUT")"
    if ! [[ "$BUILD_JOB_ID" =~ ^[0-9]+$ ]]; then
      echo "error: could not parse build job id from sbatch output: $BUILD_SUBMIT_OUT" >&2
      exit 4
    fi
  fi
fi

# ── 2. one benchmark job per rank count ────────────────────────────────────

for R in "${RANKS_ARRAY[@]}"; do
  JOB_NAME="atv_baseline_${R}r"
  RUN_DIR="$RUNS_PARENT/n${R}r"
  # If the slot already has content, refuse — preserves prior runs.
  if [[ -e "$RUN_DIR/thermo.csv" || -e "$RUN_DIR/stress_strain.csv" ]]; then
    echo "error: $RUN_DIR already has outputs -- remove or rename before re-running." >&2
    exit 5
  fi

  JOB_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/${JOB_NAME}_${TS}.XXXXXX.sbatch")"
  {
    echo "#!/usr/bin/env bash"
    echo "#SBATCH --job-name=${JOB_NAME}"
    echo "#SBATCH --account=${ACCOUNT}"
    echo "#SBATCH --partition=${PARTITION}"
    echo "#SBATCH --time=${WALLTIME}"
    echo "#SBATCH --nodes=1"
    echo "#SBATCH --ntasks=${R}"
    echo "#SBATCH --cpus-per-task=1"
    echo "#SBATCH --output=${RUNS_PARENT}/slurm_${JOB_NAME}_%j.out"
    echo "#SBATCH --error=${RUNS_PARENT}/slurm_${JOB_NAME}_%j.err"
    echo
    echo "${COMMON_ENV}"
    echo
    echo "cd '${REPO_ROOT}'"
    echo "test -x '${BIN_REL}' || { echo 'error: ${BIN_REL} missing' >&2; exit 3; }"
    echo
    echo "mkdir -p '${RUN_DIR}'"
    echo "{"
    echo "  echo \"RANKS=${R}\""
    echo "  echo \"OMP_NUM_THREADS=\${OMP_NUM_THREADS-unset}\""
    echo "  echo \"git_describe=\$(git describe --always --dirty 2>/dev/null || echo unknown)\""
    echo "  echo \"binary=${BIN_REL}\""
    echo "  echo \"binary_arch=\$(file -b '${BIN_REL}' | head -c 80)\""
    echo "  echo \"host=\$(uname -a)\""
    echo "  echo \"started=\$(date -Iseconds)\""
    echo "} > '${RUN_DIR}/run_env.txt'"
    echo
    echo "START=\$(date +%s.%N)"
    echo "srun -n ${R} '${BIN_REL}' \\"
    echo "  --out-dir '${RUN_DIR}' \\"
    echo "  --data '${DATA_FILE}' \\"
    echo "  --traj-every 10000 --rdf 10.5,200,1000"
    echo "END=\$(date +%s.%N)"
    echo "WALL=\$(awk -v s=\"\$START\" -v e=\"\$END\" 'BEGIN { printf \"%.3f\", e - s }')"
    echo "echo \"ended=\$(date -Iseconds)\" >> '${RUN_DIR}/run_env.txt'"
    echo "echo \"wall_seconds=\${WALL}\" >> '${RUN_DIR}/run_env.txt'"
    echo
    echo "# Canonicalize for downstream comparison drivers."
    echo "if [[ -f '${MOLPY_VENV_DEFAULT}/bin/activate' ]]; then"
    echo "  source '${MOLPY_VENV_DEFAULT}/bin/activate' || true"
    echo "fi"
    echo "python3 '${CANON_TOOL}' atomiverse thermo \\"
    echo "  '${RUN_DIR}/thermo.csv' \\"
    echo "  '${RUN_DIR}/thermo_canonical.csv' || true"
    echo "if [[ -s '${RUN_DIR}/stress_strain.csv' ]]; then"
    echo "  python3 '${CANON_TOOL}' atomiverse stress-strain \\"
    echo "    '${RUN_DIR}/stress_strain.csv' \\"
    echo "    '${RUN_DIR}/stress_strain_canonical.csv' || true"
    echo "fi"
  } > "$JOB_SCRIPT"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "── ${JOB_NAME} (dry-run) ─────────────────────────────────────────"
    cat "$JOB_SCRIPT"
    echo "──────────────────────────────────────────────────────────────────"
    continue
  fi

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "error: sbatch not on PATH; cannot submit. Use --dry-run to inspect." >&2
    exit 4
  fi

  DEP_FLAG=()
  if [[ -n "$BUILD_JOB_ID" ]]; then
    DEP_FLAG=(--dependency="afterok:${BUILD_JOB_ID}")
  fi
  echo "[submit] ${JOB_NAME}  nodes=1 ntasks=${R} time=${WALLTIME} account=${ACCOUNT}${BUILD_JOB_ID:+  dep=afterok:${BUILD_JOB_ID}}"
  sbatch "${DEP_FLAG[@]}" "$JOB_SCRIPT"
done
