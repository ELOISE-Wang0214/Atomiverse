#!/bin/bash -l
# CPU MD weak+strong scaling sweep: Atomiverse vs (optional) LAMMPS, 1..N MPI ranks.
#
#   sbatch benchmarks/md/common/bench_scaling_cpu.sbatch          # LJ, AV only
#   LMP=/path/to/lmp sbatch benchmarks/md/common/bench_scaling_cpu.sbatch        # + LAMMPS
#   SYS=pe CORES="1 8 32" sbatch benchmarks/md/common/bench_scaling_cpu.sbatch   # PE, fewer ranks
#   SYS=min sbatch benchmarks/md/common/bench_scaling_cpu.sbatch  # LJ minimizer (CG; MIN_OPT=lbfgs)
#
# Everything is an env-overridable variable in the CONFIG block below — no
# hardcoded repo / LAMMPS / output paths, no meta.json, no REPS loop.
# Writes plot_scaling.py-compatible results to
#   benchmarks/md/<sys>/performance/runs/<sha>_<date>/results_cpu.csv
# (schema engine,system,mode,ng,atoms,matom_step_s,vram_gb[=RSS],rep) so the
# figure is one command away — the script prints it at the end.
#SBATCH -J av_scan_cpu
#SBATCH -p cpu
#SBATCH -A naiss2026-3-381-cpu
#SBATCH -N 4
#SBATCH --exclusive
#SBATCH -t 2:00:00
# Override account/partition on the command line, e.g.:
#   sbatch -A other-proj -p other-part benchmarks/md/bench_scaling_cpu.sbatch
set -e

# ── CONFIG (override via env) ──────────────────────────────────────────────
SYS="${SYS:-lj}"                            # lj | pe | min (LJ minimizer)
CORES="${CORES:-1 8 32 128 256 512 1024}"   # MPI rank counts to sweep
STEPS="${STEPS:-50}"                        # timed MD steps
WEAK_PER_RANK="${WEAK_PER_RANK:-200000}"    # weak scaling: atoms per rank
STRONG_TOTAL="${STRONG_TOTAL:-8000000}"     # strong scaling: fixed total atoms
MIN_ITERS="${MIN_ITERS:-100}"               # SYS=min: optimizer iterations
MIN_OPT="${MIN_OPT:-cg}"                    # SYS=min: cg | lbfgs
LMP="${LMP:-}"                              # path to LAMMPS 'lmp'; empty ⇒ AV only
MOLRS="${MOLRS:-$HOME/work/molcrafts/molrs}" # molrs checkout (../molrs breaks from a worktree)

# ── ENV (auto-derived; portable) ───────────────────────────────────────────
# Under sbatch the script is copied to a spool dir, so $0 is NOT the repo path.
# SLURM_SUBMIT_DIR is the directory sbatch was invoked from — resolve the repo
# root from there (works for a normal checkout or a git worktree). Submit from
# inside the repo:  cd <repo> && sbatch benchmarks/md/common/bench_scaling_cpu.sbatch
REPO="$(cd "${SLURM_SUBMIT_DIR:-$PWD}" 2>/dev/null && git rev-parse --show-toplevel 2>/dev/null)"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$REPO"
module load buildenv-gcc/2026.03 2>/dev/null || true
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
export OMP_NUM_THREADS=1
AV="$REPO/build-cpu/benchmarks/scaling"
AVMIN="$REPO/build-cpu/benchmarks/minimize_scaling"
SR="srun --mpi=pmi2 -c 1"   # account/partition inherited from the SBATCH allocation
SYSU=$(echo "$SYS" | tr a-z A-Z)
# plot_scaling.py reads <sys>/runs/<id>/results_cpu.csv and picks the latest id;
# keep the <sha>_<date> run dir so the figure regenerates with one command.
# SYS=min reuses the lj system dir but writes results_cpu_min.csv so the MD
# plotter never pools minimizer rows with MD rows (it groups by mode,engine).
SYSDIR="$SYS"; CSV="results_cpu.csv"
[ "$SYS" = min ] && { SYSDIR=lj; CSV="results_cpu_min.csv"; }
RUN_ID="$(git rev-parse --short HEAD 2>/dev/null || echo local)_$(date +%F)"
OUT="$REPO/benchmarks/md/$SYSDIR/performance/runs/$RUN_ID/$CSV"
mkdir -p "$(dirname "$OUT")"
echo "engine,system,mode,ng,atoms,matom_step_s,vram_gb,rep" | tee "$OUT"

# ── BUILD (once, if missing) ───────────────────────────────────────────────
if [ ! -x "$AV" ]; then
  # molrs (Rust dep) needs the x86_64 toolchain; the login/compute default is
  # wrongly aarch64. Harmless if already correct.
  export RUSTUP_TOOLCHAIN="${RUSTUP_TOOLCHAIN:-stable-x86_64-unknown-linux-gnu}"
  cmake -S . -B build-cpu -DATV_BUILD_CUDA=OFF -DATV_USE_MPI=ON \
        -DATV_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release \
        -DATV_MOLRS_PATH="$MOLRS" >/dev/null
  cmake --build build-cpu -j --target scaling minimize_scaling polyethylene_tension >/dev/null
fi

# ── helpers ────────────────────────────────────────────────────────────────
cbrt(){ awk -v n="$1" 'BEGIN{printf "%d", int(n^(1/3)+0.5)}'; }
# Sum of per-rank peak VmHWM [GB] — wrap each LAMMPS rank in `/usr/bin/time -v`
# (AV self-reports the same metric in its CSV column 8).
lmp_rss(){ local e; e=$(mktemp); "$@" >/dev/null 2>"$e"
  awk -F: '/Maximum resident set size/{v=$NF;gsub(/[^0-9]/,"",v);s+=v}END{printf "%.2f",s/1048576}' "$e"; rm -f "$e"; }
# matom-step/s + atoms from a LAMMPS log's production "Loop time" line.
lmp_tput(){ awk '/^Loop time of/{t=$4;s=$9;a=$12}END{if(t>0)printf "%.1f %d",a*s/t/1e6,a}' "$1"; }

# ── AV runners ─────────────────────────────────────────────────────────────
av_lj(){ local mode=$1 p=$2 flag=$3 val=$4 row
  row=$($SR -n $p $AV $flag $val --steps $STEPS --warmup 10 --quiet | grep -E '^[0-9]+,' | tail -1) || return 0
  [ -n "$row" ] || return 0
  echo "AV,$SYSU,$mode,$p,$(echo $row|cut -d, -f2),$(awk -F, -v r="$row" 'BEGIN{split(r,a,",");printf "%.1f",a[2]*1000/a[5]/1e6}'),$(awk -F, -v r="$row" 'BEGIN{split(r,a,",");printf "%.2f",a[8]/1048576}'),1" | tee -a "$OUT"; }

# Minimizer (SYS=min): CSV row is ranks,atoms,iters,total_s,per_iter_ms,
# per_rank,mode,opt,rss_kb — same throughput formula as av_lj (col2 atoms,
# col5 per-iter ms), RSS in col 9. Matom-iter/s == Matom-feval/s (fixed-alpha
# optimizers do exactly 1 force evaluation per iteration).
av_min(){ local mode=$1 p=$2 flag=$3 val=$4 row
  row=$($SR -n $p $AVMIN $flag $val --iters $MIN_ITERS --opt $MIN_OPT --quiet | grep -E '^[0-9]+,' | tail -1) || return 0
  [ -n "$row" ] || return 0
  echo "AV,$SYSU,$mode,$p,$(echo $row|cut -d, -f2),$(awk -F, -v r="$row" 'BEGIN{split(r,a,",");printf "%.1f",a[2]*1000/a[5]/1e6}'),$(awk -F, -v r="$row" 'BEGIN{split(r,a,",");printf "%.2f",a[9]/1048576}'),1" | tee -a "$OUT"; }

av_pe(){ local mode=$1 p=$2 grid=$3 atoms=$4 rd; rd=$(mktemp -d)
  local rss; rss=$(lmp_rss $SR -n $p /usr/bin/time -v "$REPO/build-cpu/benchmarks/polyethylene_tension" \
        --backend cpu --pe-melt-grid "$grid" --mode stage:1 --scale "$(awk -v s=$STEPS 'BEGIN{printf "%.4f",s/10000}')" --out-dir "$rd" --quiet)
  local f="$rd/stage_component_timing.csv"; [ -f "$f" ] || { rm -rf "$rd"; return 0; }
  local ts st; ts=$(tail -1 "$f"|awk -F, '{print $4}'); st=$(tail -1 "$f"|awk -F, '{print $3}')
  echo "AV,$SYSU,$mode,$p,$atoms,$(awk -v a=$atoms -v s=$st -v t=$ts 'BEGIN{if(t>0)printf "%.1f",a*s/t/1e6}'),$rss,1" | tee -a "$OUT"; rm -rf "$rd"; }

# ── LAMMPS runners (only when $LMP set) ────────────────────────────────────
lmp_lj(){ [ -n "$LMP" ] || return 0; local mode=$1 p=$2 cells=$3 log; log=$(mktemp)
  local rss; rss=$(lmp_rss $SR -n $p /usr/bin/time -v "$LMP" -in "$REPO/benchmarks/md/lj/performance/in.lj_kk" -v cells $cells -v steps $STEPS -log "$log")
  local r; r=$(lmp_tput "$log"); rm -f "$log"; [ -n "$r" ] || return 0
  echo "LAMMPS,$SYSU,$mode,$p,$(echo $r|cut -d' ' -f2),$(echo $r|cut -d' ' -f1),$rss,1" | tee -a "$OUT"; }

# LAMMPS minimize: work unit is FORCE EVALUATIONS (min_style cg line-searches
# ~2 fevals/iter), parsed from the minimize stats block; Matom-feval/s is the
# number comparable to AV's Matom-iter/s.
lmp_min(){ [ -n "$LMP" ] || return 0; local mode=$1 p=$2 cells=$3 log; log=$(mktemp)
  local rss; rss=$(lmp_rss $SR -n $p /usr/bin/time -v "$LMP" -in "$REPO/benchmarks/md/lj/performance/in.lj_min_scale_kk" -v cells $cells -v iters $MIN_ITERS -log "$log")
  local r; r=$(awk '/Iterations, force evaluations/{fe=$NF} /^Loop time of/{t=$4;a=$12} END{if(t>0&&fe>0)printf "%.1f %d",a*fe/t/1e6,a}' "$log"); rm -f "$log"; [ -n "$r" ] || return 0
  echo "LAMMPS,$SYSU,$mode,$p,$(echo $r|cut -d' ' -f2),$(echo $r|cut -d' ' -f1),$rss,1" | tee -a "$OUT"; }

lmp_pe(){ [ -n "$LMP" ] || return 0; local mode=$1 p=$2 nx=$3 ny=$4 nz=$5 atoms=$6 log; log=$(mktemp)
  local rss; rss=$(cd "$REPO/benchmarks/md/pe/performance" && lmp_rss $SR -n $p /usr/bin/time -v "$LMP" -in in.pe_kk -v nx $nx -v ny $ny -v nz $nz -v steps $STEPS -log "$log")
  local r; r=$(lmp_tput "$log"); rm -f "$log"; [ -n "$r" ] || return 0
  echo "LAMMPS,$SYSU,$mode,$p,$atoms,$(echo $r|cut -d' ' -f1),$rss,1" | tee -a "$OUT"; }

# PE weak grid per rank count (atoms = nx*ny*nz*18). Strong = fixed 38^3 chains.
declare -A PE_GRID=( [1]="16 16 16" [8]="32 32 32" [32]="32 64 64" [128]="64 64 128" \
                     [256]="64 128 128" [512]="128 128 128" [1024]="128 128 256" )
declare -A PE_AT=(   [1]=73728 [8]=589824 [32]=2359296 [128]=9437184 \
                     [256]=18874368 [512]=37748736 [1024]=75497472 )

# ── SWEEP ──────────────────────────────────────────────────────────────────
for p in $CORES; do
  if [ "$SYS" = lj ]; then
    av_lj  weak   $p --atoms-per-rank $WEAK_PER_RANK; lmp_lj weak   $p "$(cbrt $((WEAK_PER_RANK*p)))"
    av_lj  strong $p --atoms          $STRONG_TOTAL ; lmp_lj strong $p "$(cbrt $STRONG_TOTAL)"
  elif [ "$SYS" = min ]; then
    av_min weak   $p --atoms-per-rank $WEAK_PER_RANK; lmp_min weak   $p "$(cbrt $((WEAK_PER_RANK*p)))"
    av_min strong $p --atoms          $STRONG_TOTAL ; lmp_min strong $p "$(cbrt $STRONG_TOTAL)"
  else
    g=${PE_GRID[$p]}; a=${PE_AT[$p]}; [ -n "$g" ] || continue; set -- $g
    av_pe  weak   $p "$1,$2,$3" $a;                   lmp_pe weak   $p $1 $2 $3 $a
    av_pe  strong $p "38,38,38" 987696;               lmp_pe strong $p 38 38 38 987696
  fi
done
echo "=== wrote $OUT ==="
# Always plot (per project convention: throughput on left-y, MEMORY on the
# right-y, TOTAL ATOMS on the top-x). Best-effort — needs python+matplotlib;
# if absent, run the printed command on the login node.
PLOT="$REPO/benchmarks/md/$SYSDIR/scaling_cpu.png"
if [ "$SYS" = min ]; then
  echo "=== SYS=min: no auto-plot (results_cpu_min.csv is analysed directly) ==="
elif python3 -c "import matplotlib" 2>/dev/null; then
  python3 "$REPO/benchmarks/md/common/plot_scaling.py" --backend cpu --id "$RUN_ID" --out "$PLOT" \
    && echo "=== wrote figure $PLOT ===" || echo "(plot failed; run the command below on the login node)"
fi
echo "=== plot:  python3 benchmarks/md/common/plot_scaling.py --backend cpu --id $RUN_ID --out scaling_cpu.png ==="
