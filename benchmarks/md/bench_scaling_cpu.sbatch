#!/bin/bash
# CPU counterpart of bench_scaling.sbatch: AV vs LAMMPS over MPI ranks (cores) on
# one exclusive cpu node. Same lib + per-system config.sh; writes
#   benchmarks/md/<system>/runs/<git-sha>_<date>/{results_cpu.csv, meta.json, logs/}
# results_cpu.csv schema: engine,system,mode,ng,atoms,matom_step_s,vram_gb
#   (ng = MPI ranks/cores; vram_gb column = peak RSS [GB]).
#   sbatch benchmarks/md/bench_scaling_cpu.sbatch lj      # or: pe
#SBATCH -J av_scan_cpu
#SBATCH -p cpu
#SBATCH -A naiss2026-3-381-cpu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 3:00:00
#SBATCH -o /home/jicli594/work/molcrafts/Atomiverse/benchmarks/md/.scanout/scancpu_%j.out
set +e
SYS="${1:?usage: sbatch bench_scaling_cpu.sbatch <lj|pe>}"
source /etc/profile 2>/dev/null || true
module load buildenv-gcc/2026.03 2>/dev/null || true
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
export OMP_NUM_THREADS=1
REPO=/home/jicli594/work/molcrafts/Atomiverse
BD="$REPO/benchmarks/md"; cd "$REPO"
source "$BD/lib_scaling.sh"; source "$BD/$SYS/config.sh"
CPU_STEPS=50

SHA=$(git rev-parse --short HEAD); DATE=$(date +%F)
RUN="$BD/$SYS/runs/${SHA}_${DATE}"; L="$RUN/logs"; mkdir -p "$L"
CSV="$RUN/results_cpu.csv"; echo "engine,system,mode,ng,atoms,matom_step_s,vram_gb" > "$CSV"

echo "== build-cpu ($SHA) =="
[ -f build-cpu/benchmarks/scaling ] || cmake -S . -B build-cpu -DATV_BUILD_CUDA=OFF -DATV_USE_MPI=ON -DATV_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release 2>&1|tail -3
cmake --build build-cpu -j 64 --target scaling polyethylene_tension 2>&1 | tail -2
AVLJ=build-cpu/benchmarks/scaling
AVPE=build-cpu/benchmarks/polyethylene_tension
LMP=/home/jicli594/work/lammps/build/lmp
cat > "$RUN/meta_cpu.json" <<META
{ "system":"$SYS","av_git_sha":"$SHA","date":"$DATE","host":"$(hostname)",
  "backend":"cpu","cores":"$CPU_CORES","steps":$CPU_STEPS,"lammps":"$LMP" }
META

if [ "$SYS" = lj ]; then
  avlj(){ local mode=$1 p=$2 flag=$3 val=$4 l
    l=$(srun --mpi=pmi2 -n $p -c 1 $AVLJ $flag $val --steps $CPU_STEPS --warmup 10 --quiet 2>/dev/null|grep -E '^[0-9]+,[0-9]+,'|tail -1)
    [ -z "$l" ] && { echo "AV LJ $mode p=$p FAILED"; return; }
    local at=$(echo $l|cut -d, -f2) tps=$(echo $l|cut -d, -f5) rss=$(echo $l|cut -d, -f8)
    local mas=$(awk -v a=$at -v t=$tps 'BEGIN{if(t>0)printf "%.1f",a*1000/t/1e6}')
    echo "AV,LJ,$mode,$p,$at,$mas,$(awk -v r=$rss 'BEGIN{printf "%.2f",r/1048576}')">>"$CSV"
    echo "AV LJ $mode p=$p -> $mas Matom/s"; }
  lmplj(){ local mode=$1 p=$2 cells=$3; local log="$L/cpu_lmp_${mode}_${p}.log"
    local rss; rss=$(pv_cpu srun --mpi=pmi2 -n $p -c 1 $LMP -in $BD/lj/in.lj_kk -v cells $cells -v steps $CPU_STEPS -log $log)
    local r; r=$(lmp_throughput "$log") || { echo "LAMMPS LJ $mode p=$p FAILED"; return; }
    echo "LAMMPS,LJ,$mode,$p,$(echo $r|cut -d' ' -f2),$(echo $r|cut -d' ' -f1),$rss">>"$CSV"
    echo "LAMMPS LJ $mode p=$p -> $(echo $r|cut -d' ' -f1) Matom/s ${rss}GB"; }
  for p in $CPU_CORES; do avlj weak $p --atoms-per-rank $CPU_WEAK_PER_CORE; lmplj weak $p $(cbrt $((CPU_WEAK_PER_CORE*p))); done
  for p in $CPU_CORES; do avlj strong $p --atoms $CPU_STRONG_TOTAL; lmplj strong $p $(cbrt $CPU_STRONG_TOTAL); done

elif [ "$SYS" = pe ]; then
  PESC=$(awk -v s=$CPU_STEPS 'BEGIN{printf "%.4f", s/10000}')
  avpe(){ local mode=$1 p=$2 grid=$3 atoms=$4; local rd="$L/cpu_av_${mode}_${p}"; mkdir -p "$rd"
    local rss; rss=$(pv_cpu srun --mpi=pmi2 -n $p -c 1 $AVPE --backend cpu --pe-melt-grid "$grid" --mode stage:1 --scale $PESC --out-dir "$rd" --quiet)
    local f="$rd/stage_component_timing.csv"; [ -f "$f" ] || { echo "AV PE $mode p=$p FAILED"; return; }
    local ts=$(tail -1 "$f"|awk -F, '{print $4}') st=$(tail -1 "$f"|awk -F, '{print $3}')
    local mas=$(awk -v a=$atoms -v s=$st -v t=$ts 'BEGIN{if(t>0)printf "%.1f",a*s/t/1e6}')
    echo "AV,PE,$mode,$p,$atoms,$mas,$rss">>"$CSV"; echo "AV PE $mode p=$p -> $mas Matom/s ${rss}GB"; }
  lmppe(){ local mode=$1 p=$2 nx=$3 ny=$4 nz=$5 atoms=$6; local log="$L/cpu_lmp_${mode}_${p}.log"
    local rss; rss=$( cd "$BD/pe" && pv_cpu srun --mpi=pmi2 -n $p -c 1 $LMP -in in.pe_kk -v nx $nx -v ny $ny -v nz $nz -v steps $CPU_STEPS -log $log )
    local r; r=$(lmp_throughput "$log") || { echo "LAMMPS PE $mode p=$p FAILED"; return; }
    echo "LAMMPS,PE,$mode,$p,$atoms,$(echo $r|cut -d' ' -f1),$rss">>"$CSV"
    echo "LAMMPS PE $mode p=$p -> $(echo $r|cut -d' ' -f1) Matom/s ${rss}GB"; }
  for p in $CPU_CORES; do g=${PE_CPU_WEAK_GRID[$p]}; a=${PE_CPU_WEAK_ATOMS[$p]}; set -- $g
    avpe weak $p "$1,$2,$3" $a; lmppe weak $p $1 $2 $3 $a; done
  set -- $PE_CPU_STRONG_GRID
  for p in $CPU_CORES; do avpe strong $p "$1,$2,$3" $PE_CPU_STRONG_ATOMS; lmppe strong $p $1 $2 $3 $PE_CPU_STRONG_ATOMS; done
fi
echo "=== $CSV ==="; cat "$CSV"; echo "===== CPU SCAN DONE ($SYS) ====="
