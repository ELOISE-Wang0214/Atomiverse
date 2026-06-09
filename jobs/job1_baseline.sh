#!/bin/bash -l
# Baseline: build current worktree HEAD on GH200, reproduce LJ 1-4 GPU
# throughput (the 3->4 collapse), PE phase profile (AV_MG_PROFILE), and the
# multi-GPU parity test. All output under jobs/out1/.
set -x
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
WT=/home/jicli594/work/molcrafts/Atomiverse/.claude/worktrees/cuda-mg-overlap
OUT=$WT/jobs/out1; mkdir -p $OUT
cd $WT

LAPACK="$OPENBLAS_ROOT/lib/libopenblas.so"
[ -f "$OPENBLAS_ROOT/lib/liblapack.so" ] && LAPACK="$OPENBLAS_ROOT/lib/liblapack.so"
cmake -B build-cuda -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON -DATV_USE_MPI=OFF \
  -DATV_MOLRS_PATH=$HOME/work/molcrafts/molrs \
  -DBLAS_LIBRARIES="$OPENBLAS_ROOT/lib/libopenblas.so" \
  -DLAPACK_LIBRARIES="$LAPACK" \
  2>&1 | tail -5 || exit 1
cmake --build build-cuda -j 48 --target bench_multi_gpu_scaling polyethylene_tension test_cuda_multi_gpu 2>&1 | tail -5 || exit 1

B=$WT/build-cuda/benchmarks
vis(){ seq -s, 0 $(($1-1)); }

echo "===== multi-GPU parity test (baseline) ====="
$WT/build-cuda/tests/test_cuda_multi_gpu > $OUT/test_multi_gpu.log 2>&1
echo "test_cuda_multi_gpu exit=$? (log: $OUT/test_multi_gpu.log)"
tail -3 $OUT/test_multi_gpu.log

echo "===== LJ weak 1-4 GPU (8M/card, 200 steps, reneigh 20) ====="
for ng in 1 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling $((8000000*ng)) $ng 200 20 2>$OUT/lj_weak_$ng.err | grep '^CSV' | tee -a $OUT/lj_weak.csv
done

echo "===== PE phase profile (AV_MG_PROFILE=1, stage:1 500 steps) ====="
pe(){ local ng=$1 grid=$2; local rd=$OUT/pe_$ng; mkdir -p $rd
  CUDA_VISIBLE_DEVICES=$(vis $ng) AV_MG_PROFILE=1 $B/polyethylene_tension \
    --backend cuda --gpus $ng --pe-melt-grid "$grid" --mode stage:1 \
    --scale 0.05 --out-dir $rd --quiet 2> $OUT/pe_prof_$ng.err
  echo "-- PE ng=$ng stage_timing:"; tail -2 $rd/stage_timing.csv
  echo "-- PE ng=$ng phase profile:"; grep -A 20 'phase profile' $OUT/pe_prof_$ng.err; }
pe 1 "60,60,60"
pe 2 "76,76,76"
pe 4 "96,96,96"
echo "===== BASELINE DONE ====="
