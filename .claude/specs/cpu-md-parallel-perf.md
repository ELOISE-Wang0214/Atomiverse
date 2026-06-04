---
title: CPU MD multi-rank (MPI) performance — close the parallel gap to LAMMPS
status: code-complete
created: 2026-06-04
layer: cpu/md (+ cpu/md/mpi, cpu/comm)
---

# CPU MD 并行（多 rank）性能优化 —— 追平/超过 LAMMPS

## 综述 (Summary)

单 rank 的 CPU MD 已经比 LAMMPS 快约 14%（见 `cpu-md-lammps-perf-assessment`
的结论 + 已落地的 half-stencil 邻居表、LJ shift、AVX-512 向量化累加）。但这些
胜利**大多被门控在单 rank**，没有传导到并行：实测全量 5 阶段 MSU-PE（native
fabric，1 线程/rank）**NEW-AV 32r=321s / 64r=225s vs LAMMPS 164s/86s**（仍慢
~2–2.6×；64r 因向量化累加改善了 33%，32r 基本持平）。本规范把单 rank 的内核
优化扩展到多 rank，并削减通信开销，目标是把并行差距追平并力争反超 LAMMPS。

## 背景数据（已测，stage5 占全程 52%）

| ranks | OLD-AV | NEW-AV | LAMMPS |
|------:|------:|------:|------:|
| 32    | 324s  | 321s  | 164s  |
| 64    | 336s  | 225s  | 86s   |

多 rank 瓶颈归因（stage5）：**Neigh 42–47%**（邻居构建在 MPI 下走的是**未优化的
full-stencil 路径**——half-stencil + per-cell-shift 门控为 `!tag_dedup`=单 rank）
+ **Comm 22–28%**。

**关键否定结论：负载不均不是瓶颈。** RCB 负载均衡实测：把 64r 的
nlocal 49/193 改善到 127/177，但 64r 总时间**零加速**、32r 反而因再均衡开销
（每 500 步 histogram+migration+邻居重建）**变慢到 502s**。已从 benchmark 回退。
⇒ 多 rank 步长由 Comm + Neigh + 全局同步决定，不由原子数倾斜决定。

## Domain basis (science.required)

所有优化保持物理不变：与 `lammps_msu_paper.in` 一致的 5 阶段协议，结果落在
`msu-pe-validation` 容差内（能量/密度/应力-应变）。参考：LAMMPS `Comm`/
`Neighbor` 多 rank 设计（newton-on 半邻居表 + 反向力归约 + ghost 镜像）、
Plimpton 1995（comm ∝ surface），本仓 `.claude/notes/cpu_md_perf_review_2026-05-24.md`
（HIGH/MED 通信项）+ `.claude/notes/cpu_md_perf_assessment_2026-06.md`。

## Design — 三个杠杆，按收益排序

**杠杆 1（最大）：把 half-stencil + per-cell-shift 邻居构建扩展到多 rank.**
`cell_list.cpp` 的 `refresh_half` 串行快路径目前 `use_half = !tag_dedup && ...`，
多 rank（tag_dedup=true）退回 full-stencil（27 cell）+ while-loop min-image。
多 rank 用 ghost halo（`axis_extended_`，ghost 在真实位置），所以：
- 让 half-stencil 在多 rank 也生效（home cell `j>i` + forward cell 全 j），
  保留 ghost 的 `tag_ptr[j]<tag_ptr[i]` 去重；处理 ghost/extended 轴的边界
  （extended 轴 `continue` 出界，不 wrap，shift=0）。
- 验证 small-grid/aliasing 保护（≥3 cell/轴）在多 rank 仍成立；否则保留
  full-stencil 回退。
目标：Neigh 42–47% → ~半。

**杠杆 2：多 rank 也用存好的 per-pair shift（去掉 LJ 内核 min-image）.**
`pair_lj_cut.cpp` 的 `k_has_shifts` 现门控为单 rank。多 rank 下 ghost 在真实
位置、min-image 对 ghost 对是 no-op，存的 shift 对这些对恰为 0（非分解轴的
wrap 对则非 0，shift 已正确编码）。把门控从 `!(dd && mpi_size>1)` 放开为仅
依赖 `has_shifts`（serial build 已填充），多 rank Pair 即跳过 min-image。需
8-rank 正确性验证（ghost 对 dx 一致）。

**杠杆 3：通信路径收紧**（perf-review MED 项，按需）：
- `comm/halo_exchange.cpp` send-list 扫描 OMP/SoA 化；recv-buf 不再每步
  `assign(0.0)` 清零；
- `comm/atom_migration.cpp` 复用 member-state 缓冲（去掉每次 `vector<vector<int>>`
  分配）；
- 把残余的独立 `MPI_Allreduce` 合并进 `ReduceEnergy` 18-double 批；
- `compute_bonded_extension_local` 结果按 `neighbor_generation` 缓存。

## Files

- `src/cpu/md/cell_list.cpp`（杠杆 1：多 rank half-stencil + shift；最关键）。
- `src/cpu/potentials/pair_lj_cut.cpp`（杠杆 2：放开 `k_has_shifts` 多 rank 门控）。
- `src/cpu/comm/halo_exchange.cpp`、`src/cpu/comm/atom_migration.cpp`、
  `src/cpu/md/mpi/reduce.cpp`（杠杆 3，按测量取舍）。
- 测量：`benchmarks/md/.../perf_assessment_2026-06-04/` 的 32r/64r 全量运行 +
  `stage_component_timing.csv`（已有计时基建 `core/timer.h` 等）。

## Tasks

- [x] Extend half-stencil + per-cell-shift neighbor build to the multi-rank path (ghost cells, tag_dedup), with small-grid fallback
- [x] Ungate the LJ per-pair-shift fast path for multi-rank (drop the single-rank guard; rely on has_shifts)
- [ ] Comm-path tightening per measurement (halo send-list, recv-buf no-zero, migration buffers, batched reduce) — DEFERRED: spec-gated on measurement; the comm bottleneck is only re-attributed by the operator-action ac-004 SLURM run, which has not yet happened. Levers 1+2 already target the dominant Neigh (42–47%) cost; tightening working MPI comm code blind, for unmeasured gain, is out of the safe diff budget. Revisit after ac-004 shows whether Comm is still the gate.
- [ ] Re-measure 32r/64r full-run + per-component attribution; confirm Neigh fraction drops and total beats/approaches LAMMPS — operator-action (ac-004), SLURM compute-node timing; cannot run from this harness.
- [x] Verify 8-rank correctness (test_pe_msu_pipeline green) and physics within msu-pe-validation tolerances — 8-rank correctness GREEN (test_mpi_pe_msu_pipeline + new test_mpi_lj_energy_bruteforce 1r/8r brute-force equality + lj_nve/pe_full/pe_chains_8rank all pass). Physics within msu-pe-validation tolerances (ac-005) is the operator/bench cluster gate.

## Testing

8-rank `tests/mpi/test_pe_msu_pipeline` (bonded + LJ exclusions + halo + migration
+ NVT + barostat + deform) must stay green; 1r-vs-8r bonded/LJ energy diff within
tolerance; operator-action 32r/64r SLURM full runs for timing (native fabric,
`srun --mpi=pmi2`/mpiexec, exclude flaky nodes); formal `msu-pe-validation`
energy/density/stress cluster gate before commit.

## Out of scope

- RCB / dynamic load balancing (measured: not the gate for this workload; reverted).
- Single-rank kernels (already beat LAMMPS; this spec is multi-rank).
- CUDA backend; algorithm/physics changes.
