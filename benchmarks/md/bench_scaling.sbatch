#!/bin/bash
# ONE canonical multi-GPU scaling scan (AV vs LAMMPS-KOKKOS, GH200). Run per
# system; reuse this exact script every time so runs are comparable.
#
#   sbatch benchmarks/md/bench_scaling.sbatch lj      # or: pe
#
# Builds the CURRENT committed Atomiverse, then writes
#   benchmarks/md/<system>/runs/<git-sha>_<date>/{results.csv, meta.json, logs/}
# results.csv schema: engine,system,mode,ng,atoms,matom_step_s,vram_gb
# Plot with plot_scaling.py. See README.md.
#SBATCH -J av_scan
#SBATCH -p gpu
#SBATCH -A naiss2026-4-715-gpu
#SBATCH --gres=gpu:4
#SBATCH -N 1
#SBATCH -t 0:45:00
#SBATCH -o /home/jicli594/work/molcrafts/Atomiverse/benchmarks/md/.scanout/scan_%j.out
set +e
SYS="${1:?usage: sbatch bench_scaling.sbatch <lj|pe>}"
source /etc/profile 2>/dev/null || true
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
export MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=generic
REPO=/home/jicli594/work/molcrafts/Atomiverse
BD="$REPO/benchmarks/md"; cd "$REPO"
source "$BD/lib_scaling.sh"
source "$BD/$SYS/config.sh"

SHA=$(git rev-parse --short HEAD); DATE=$(date +%F)
RUN="$BD/$SYS/runs/${SHA}_${DATE}"; L="$RUN/logs"; mkdir -p "$L"
CSV="$RUN/results_gpu.csv"; echo "engine,system,mode,ng,atoms,matom_step_s,vram_gb" > "$CSV"

echo "== build committed HEAD ($SHA) =="
cmake --build build-cuda -j 32 --target bench_multi_gpu_scaling polyethylene_tension 2>&1 | tail -2
LJ=build-cuda/benchmarks/bench_multi_gpu_scaling
PE=build-cuda/benchmarks/polyethylene_tension
LMP=/home/jicli594/work/lammps/build-kk-mpi-mixed/lmp

# provenance — makes each run auditable
cat > "$RUN/meta.json" <<META
{
  "system": "$SYS",
  "av_git_sha": "$SHA",
  "av_git_subject": "$(git log -1 --format=%s | sed 's/"/\\"/g')",
  "date": "$DATE",
  "host": "$(hostname)",
  "gpu": "$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)",
  "module": "GPU/buildenv-gcccuda/2026.03-cu13.0",
  "lammps": "$LMP",
  "ngs": "$NGS", "steps": $STEPS, "reneigh": $RENEIGH
}
META

vis(){ seq -s, 0 $(($1-1)); }

if [ "$SYS" = lj ]; then
  avlj(){ local mode=$1 ng=$2 N=$3 l; l=$(CUDA_VISIBLE_DEVICES=$(vis $ng) $LJ $N $ng $STEPS $RENEIGH 2>/dev/null|grep '^CSV,')
    echo "AV,LJ,$mode,$ng,$(echo $l|cut -d, -f2),$(echo $l|cut -d, -f7),$(echo $l|awk -F, '{printf "%.2f",$9/1024}')">>"$CSV"
    echo "AV LJ $mode ng=$ng -> $(echo $l|cut -d, -f7) Matom/s"; }
  lmplj(){ local mode=$1 ng=$2 cells=$3; local log="$L/lmp_${mode}_${ng}.log"
    local vr; vr=$(pv srun --mpi=pmi2 -n $ng $LMP -k on g $ng -sf kk -pk kokkos gpu/aware on -in $BD/lj/in.lj_kk -v cells $cells -v steps $STEPS -log $log)
    local r; r=$(lmp_throughput "$log") || { echo "LAMMPS LJ $mode ng=$ng FAILED"; tail -3 "$log"; return; }
    echo "LAMMPS,LJ,$mode,$ng,$(echo $r|cut -d' ' -f2),$(echo $r|cut -d' ' -f1),$vr">>"$CSV"
    echo "LAMMPS LJ $mode ng=$ng -> $(echo $r|cut -d' ' -f1) Matom/s ${vr}GB"; }
  for ng in $NGS; do avlj weak $ng $((WEAK_PER_GPU*ng)); lmplj weak $ng $(cbrt $((WEAK_PER_GPU*ng))); done
  for ng in $NGS; do avlj strong $ng $STRONG_TOTAL; lmplj strong $ng $(cbrt $STRONG_TOTAL); done

elif [ "$SYS" = pe ]; then
  avpe(){ local mode=$1 ng=$2 grid=$3 atoms=$4; local rd="$L/av_${mode}_${ng}"; mkdir -p "$rd"
    local vr; vr=$(CUDA_VISIBLE_DEVICES=$(vis $ng) pv srun --mpi=pmi2 -n 1 $PE --backend cuda --gpus $ng --pe-melt-grid "$grid" --mode stage:1 --scale $PESCALE --out-dir "$rd" --quiet)
    local sps; sps=$(tail -1 "$rd/stage_timing.csv" 2>/dev/null|awk -F, '{print $5}')
    echo "AV,PE,$mode,$ng,$atoms,$(awk -v a=$atoms -v s=$sps 'BEGIN{printf "%.1f",a*s/1e6}'),$vr">>"$CSV"
    echo "AV PE $mode ng=$ng -> $(awk -v a=$atoms -v s=$sps 'BEGIN{printf "%.1f",a*s/1e6}') Matom/s ${vr}GB"; }
  lmppe(){ local mode=$1 ng=$2 nx=$3 ny=$4 nz=$5 atoms=$6; local log="$L/lmp_${mode}_${ng}.log"
    local vr; vr=$( cd "$BD/pe" && pv srun --mpi=pmi2 -n $ng $LMP -k on g $ng -sf kk -pk kokkos gpu/aware on -in in.pe_kk -v nx $nx -v ny $ny -v nz $nz -v steps $STEPS -log $log )
    local r; r=$(lmp_throughput "$log") || { echo "LAMMPS PE $mode ng=$ng FAILED"; tail -3 "$log"; return; }
    echo "LAMMPS,PE,$mode,$ng,$atoms,$(echo $r|cut -d' ' -f1),$vr">>"$CSV"
    echo "LAMMPS PE $mode ng=$ng -> $(echo $r|cut -d' ' -f1) Matom/s ${vr}GB"; }
  for ng in $NGS; do g=PE_WEAK_GRID_$ng; a=PE_WEAK_ATOMS_$ng; set -- ${!g}
    avpe weak $ng "$1,$2,$3" ${!a}; lmppe weak $ng $1 $2 $3 ${!a}; done
  set -- $PE_STRONG_GRID
  for ng in $NGS; do avpe strong $ng "$1,$2,$3" $PE_STRONG_ATOMS; lmppe strong $ng $1 $2 $3 $PE_STRONG_ATOMS; done
fi

echo "=== $RUN/results.csv ==="; cat "$CSV"
echo "===== SCAN DONE ($SYS, $SHA) ====="
