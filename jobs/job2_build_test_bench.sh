#!/bin/bash -l
# Build the async-halo worktree, run the CUDA test suite, then LJ/PE
# multi-GPU throughput + PE phase profile. Output under jobs/out2/.
set -x
module load GPU/buildenv-gcccuda/2026.03-cu13.0
export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
WT=/home/jicli594/work/molcrafts/Atomiverse/.claude/worktrees/cuda-mg-overlap
OUT=$WT/jobs/out2; mkdir -p $OUT
cd $WT

LAPACK="$OPENBLAS_ROOT/lib/libopenblas.so"
[ -f "$OPENBLAS_ROOT/lib/liblapack.so" ] && LAPACK="$OPENBLAS_ROOT/lib/liblapack.so"
cmake -B build-cuda -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON \
  -DATV_BUILD_BENCHMARKS=ON -DATV_USE_MPI=OFF \
  -DATV_MOLRS_PATH=$HOME/work/molcrafts/molrs \
  -DBLAS_LIBRARIES="$OPENBLAS_ROOT/lib/libopenblas.so" \
  -DLAPACK_LIBRARIES="$LAPACK" > $OUT/configure.log 2>&1 || { tail -20 $OUT/configure.log; exit 1; }

if ! cmake --build build-cuda -j 48 > $OUT/build.log 2>&1; then
  echo "===== BUILD FAILED ====="
  grep -E "error|Error" $OUT/build.log | head -40
  exit 1
fi
echo "===== BUILD OK ====="

echo "===== ctest (full suite) ====="
ctest --test-dir build-cuda --output-on-failure -j 4 > $OUT/ctest.log 2>&1
CT=$?
tail -15 $OUT/ctest.log
echo "ctest exit=$CT"

B=$WT/build-cuda/benchmarks
vis(){ seq -s, 0 $(($1-1)); }

echo "===== LJ weak 1-4 GPU (8M/card, 200 steps, reneigh 20) ====="
for ng in 1 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling $((8000000*ng)) $ng 200 20 2>$OUT/lj_weak_$ng.err | grep '^CSV' | tee -a $OUT/lj_weak.csv
done
echo "===== LJ strong 32M total, 2-4 GPU ====="
for ng in 2 3 4; do
  CUDA_VISIBLE_DEVICES=$(vis $ng) $B/bench_multi_gpu_scaling 32000000 $ng 200 20 2>$OUT/lj_strong_$ng.err | grep '^CSV' | tee -a $OUT/lj_strong.csv
done

echo "===== PE weak (stage:1, 500 steps) ====="
pe(){ local ng=$1 grid=$2 prof=$3; local rd=$OUT/pe_$ng; mkdir -p $rd
  CUDA_VISIBLE_DEVICES=$(vis $ng) AV_MG_PROFILE=$prof $B/polyethylene_tension \
    --backend cuda --gpus $ng --pe-melt-grid "$grid" --mode stage:1 \
    --scale 0.05 --out-dir $rd --quiet 2> $OUT/pe_prof_$ng.err
  echo "-- PE ng=$ng stage_timing:"; tail -2 $rd/stage_timing.csv
  grep -A 20 'phase profile' $OUT/pe_prof_$ng.err || true; }
pe 1 "60,60,60" ""
pe 2 "76,76,76" ""
pe 4 "96,96,96" ""
echo "===== PE 4-GPU phase profile (separate run, AV_MG_PROFILE distorts overlap) ====="
rm -rf $OUT/pe_4p; mkdir -p $OUT/pe_4p
CUDA_VISIBLE_DEVICES=0,1,2,3 AV_MG_PROFILE=1 $B/polyethylene_tension \
  --backend cuda --gpus 4 --pe-melt-grid "96,96,96" --mode stage:1 \
  --scale 0.05 --out-dir $OUT/pe_4p --quiet 2> $OUT/pe_prof_4p.err
grep -A 20 'phase profile' $OUT/pe_prof_4p.err || true
echo "===== JOB2 DONE ====="
