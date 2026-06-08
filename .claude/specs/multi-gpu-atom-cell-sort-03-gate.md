---
title: Multi-GPU atom cell-sort — beneficial gate + LAMMPS-beating benchmark
status: code-complete
created: 2026-06-08
chain: multi-gpu-atom-cell-sort
link: 3
---

# Multi-GPU atom cell-sort — gate + benchmark (link 3/3)

## Summary
链条第 3 环(收尾 + 性能合同)。加入 beneficial-or-skip 启发式(只在 build/force 收益大于 ~2× 置换访存代价时
排序),确认 PE 生产负载不回归,并在 GH200 上证明:多卡 LJ build even-count 悬崖消除,且 AV 多卡 LJ 32M
强扩展在 **每个** 卡数 1/2/3/4 上都**超过 LAMMPS-KOKKOS**。

## Domain basis
性能优化,力不变(link-02 已逐字节 parity)。本环只衡量吞吐与代价权衡。

## Design
- **启发式**:当每卡子域较小(奇数卡已有 ~66% L1 命中、不需要排序)或排序的 ~2× per-atom 数据搬运超过
  build/force 节省时,跳过排序。门控开启后 parity 仍须通过。
- **性能合同**:`bench_multi_gpu_scaling.cu` @ 32M,对照 LAMMPS `in.lj_kk`(`srun --mpi=pmi2`)。
  基线(link-02 前,float4-only):AV ng1=2688 / ng2=2041 / ng3=5146 / ng4=3217;LAMMPS ng1=1547 /
  ng2=2933 / ng3=4014 / ng4=5442 Matom-step/s。目标:AV 在 1/2/3/4 全面 > LAMMPS。
- 非回归:`polyethylene_tension --backend cuda --gpus N`(N=2,4)吞吐不低于变更前。

## Implementation notes (post-impl)
- Also fixed `classify_owner` atomic contention (per-block shared-memory partial counts) — a migrate win that lifts ng2/ng3 further and speeds the PE migrate too; parity (NveMigrationParity) green.
- PE production multi-GPU is FASTER with the sort: +42% (2 GPU) / +32% (4 GPU).

## Files
- `src/cuda/md/multi_gpu.cu`:beneficial-or-skip 启发式。
- `benchmarks/bench_multi_gpu_scaling.cu`:性能证明(已存在)。
- `benchmarks/md/reference/lj_scaling/in.lj_kk`:LAMMPS 对照(已存在)。

## Tasks
- [x] Implement beneficial-or-skip heuristic (skip single-GPU + per-GPU<2M owned; AV_CELL_SORT_FORCE override) + reused per-device buffers.
- [x] GREEN: parity passes BOTH gate-skip (18/18) and forced-sort/reused-buffers (18/18 + 9/9 md) — job 11187.
- [x] Bench GH200: build_verlet_direct_k 71.8ms→28.8ms @16M (even-count cliff largely removed).
- [~] Bench GH200 32M (sort + classify-atomic fix): AV beats LAMMPS ng=1(1.73x)/2(1.17x)/3(1.36x); ng=4 4550 vs 5442 (0.84x) — build (R_nlist 55%) now dominates, needs the out-of-scope cluster rewrite to cross.
- [x] Bench GH200: PE multi-GPU NOT regressed — actually +42% (2 GPU) / +32% (4 GPU) with sort (job 11203).
- [x] Hygiene: simplify ran clean.

## Testing
GH200 sbatch(partition gpu). parity gtest + `bench_multi_gpu_scaling` + LAMMPS `lmp -k on g N -sf kk` via
`srun --mpi=pmi2`(4 卡 `gpu/aware on` 需 `MPIR_CVAR_CH4_IPC_GPU_HANDLE_CACHE=generic`)。

## Out of scope
GROMACS 式 8-原子 cluster list(更大重写,作为备选/未来工作记录,不在本链)。
