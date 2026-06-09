#!/bin/bash -l
# MPI=ON rebuild (the --pe-melt-grid generator is #ifdef ATV_HAS_MPI), full
# ctest, LJ+PE throughput, and AV_MG_PROFILE diagnostics. Output: jobs/out3/.
set -x
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
export MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=generic
WT=/home/jicli594/work/molcrafts/Atomiverse/.claude/worktrees/cuda-mg-overlap
OUT=$WT/jobs/out3; mkdir -p $OUT
cd $WT

LAPACK="$OPENBLAS_ROOT/lib/libopenblas.so"
[ -f "$OPENBLAS_ROOT/lib/liblapack.so" ] && LAPACK="$OPENBLAS_ROOT/lib/liblapack.so"
cmake -B build-cuda -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON \
  -DATV_BUILD_BENCHMARKS=ON -DATV_USE_MPI=ON \
  -DATV_MOLRS_PATH=$HOME/work/molcrafts/molrs \
  -DBLAS_LIBRARIES="$OPENBLAS_ROOT/lib/libopenblas.so" \
  -DLAPACK_LIBRARIES="$LAPACK" > $OUT/configure.log 2>&1 || { tail -20 $OUT/configure.log; exit 1; }

if ! cmake --build build-cuda -j 48 > $OUT/build.log 2>&1; then
  echo "===== BUILD FAILED ====="
  grep -E "error|Error" $OUT/build.log | head -40
  exit 1
fi
echo "===== BUILD OK ====="

echo "===== ctest ====="
ctest --test-dir build-cuda --output-on-failure -j 4 > $OUT/ctest.log 2>&1
echo "ctest exit=$?"
grep -E "tests passed|tests failed" $OUT/ctest.log
grep -A5 "The following tests FAILED" $OUT/ctest.log || true

B=$WT/build-cuda/benchmarks
vis(){ seq -s, 0 $(($1-1)); }

echo "===== LJ weak 1-4 + strong (no profile) ====="
for ng in 1 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling $((8000000*ng)) $ng 200 20 2>$OUT/lj_weak_$ng.err | grep '^CSV' | tee -a $OUT/lj_weak.csv
done
for ng in 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling 32000000 $ng 200 20 2>$OUT/lj_strong_$ng.err | grep '^CSV' | tee -a $OUT/lj_strong.csv
done

echo "===== LJ phase profile ng=3 vs ng=4 (weak) ====="
for ng in 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) AV_MG_PROFILE=1 $B/bench_multi_gpu_scaling $((8000000*ng)) $ng 100 20 2> $OUT/lj_prof_$ng.err >/dev/null
  echo "-- LJ ng=$ng profile:"; grep -E "MG_ROWS|_vv_|_forward|_reneighbor|_force|R_|TOTAL" $OUT/lj_prof_$ng.err | head -25
done

echo "===== PE throughput (srun, MPI singleton) ====="
pe(){ local ng=$1 grid=$2; local rd=$OUT/pe_$ng; mkdir -p $rd
  CUDA_VISIBLE_DEVICES=$(vis $ng) srun --mpi=pmi2 -n 1 $B/polyethylene_tension \
    --backend cuda --gpus $ng --pe-melt-grid "$grid" --mode stage:1 \
    --scale 0.05 --out-dir $rd --quiet 2> $OUT/pe_$ng.err
  echo "-- PE ng=$ng:"; tail -1 $rd/stage_timing.csv; }
pe 1 "60,60,60"
pe 2 "76,76,76"
pe 3 "87,87,88"
pe 4 "96,96,96"

echo "===== PE phase profile ng=4 ====="
rm -rf $OUT/pe_4p; mkdir -p $OUT/pe_4p
CUDA_VISIBLE_DEVICES=0,1,2,3 srun --mpi=pmi2 -n 1 --export=ALL,AV_MG_PROFILE=1 \
  $B/polyethylene_tension --backend cuda --gpus 4 --pe-melt-grid "96,96,96" \
  --mode stage:1 --scale 0.05 --out-dir $OUT/pe_4p --quiet 2> $OUT/pe_prof_4.err
grep -E "MG_ROWS|phase profile" -A 22 $OUT/pe_prof_4.err | head -35
echo "===== JOB3 DONE ====="
