#!/bin/bash -l
# Validate the O(migrants) migrate rewrite: rebuild, CUDA-test subset + the
# MPI-free correctness set, LJ 1-4 GPU + PE 1-4 GPU throughput, LJ ng=4
# profile re-check. Output: jobs/out4/.
set -x
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
export MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=generic
WT=/home/jicli594/work/molcrafts/Atomiverse/.claude/worktrees/cuda-mg-overlap
OUT=$WT/jobs/out4; mkdir -p $OUT
cd $WT

if ! cmake --build build-cuda -j 48 > $OUT/build.log 2>&1; then
  echo "===== BUILD FAILED ====="
  grep -E "error|Error" $OUT/build.log | head -40
  exit 1
fi
echo "===== BUILD OK ====="

echo "===== ctest (CUDA subset) ====="
ctest --test-dir build-cuda --output-on-failure -R '^test_cuda' > $OUT/ctest_cuda.log 2>&1
echo "ctest_cuda exit=$?"
grep -E "tests passed|tests failed" $OUT/ctest_cuda.log
grep -A6 "The following tests FAILED" $OUT/ctest_cuda.log || true

B=$WT/build-cuda/benchmarks
vis(){ seq -s, 0 $(($1-1)); }

echo "===== LJ weak 1-4 + strong 2-4 ====="
for ng in 1 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling $((8000000*ng)) $ng 200 20 2>$OUT/lj_weak_$ng.err | grep '^CSV' | tee -a $OUT/lj_weak.csv
done
for ng in 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling 32000000 $ng 200 20 2>$OUT/lj_strong_$ng.err | grep '^CSV' | tee -a $OUT/lj_strong.csv
done

echo "===== LJ ng=4 phase profile (post-migrate-fix) ====="
CUDA_VISIBLE_DEVICES=0,1,2,3 AV_MG_PROFILE=1 $B/bench_multi_gpu_scaling 32000000 4 100 20 2> $OUT/lj_prof_4.err >/dev/null
grep -E "^  [0-9A-Za-z_]+ +[0-9.]+ s|TOTAL" $OUT/lj_prof_4.err | head -14

echo "===== PE 1-4 GPU ====="
pe(){ local ng=$1 grid=$2; local rd=$OUT/pe_$ng; mkdir -p $rd
  CUDA_VISIBLE_DEVICES=$(vis $ng) srun --mpi=pmi2 -n 1 $B/polyethylene_tension \
    --backend cuda --gpus $ng --pe-melt-grid "$grid" --mode stage:1 \
    --scale 0.05 --out-dir $rd --quiet 2> $OUT/pe_$ng.err
  echo "-- PE ng=$ng:"; tail -1 $rd/stage_timing.csv; }
pe 1 "60,60,60"
pe 2 "76,76,76"
pe 3 "87,87,88"
pe 4 "96,96,96"
echo "===== JOB4 DONE ====="
