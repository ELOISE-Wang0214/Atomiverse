---
title: CPU MD vs LAMMPS performance assessment — attribute the gap, then decide refactor-vs-optimize
status: code-complete
created: 2026-06-04
layer: cpu/md (assessment only — no production code change)
---

# CPU MD ↔ LAMMPS 性能评估（先诊断，后定方案）

## 综述 (Summary)

本规范 **不实现优化**。它先量化 Atomiverse CPU/MPI MD 相对 LAMMPS 在 MSU 非晶
聚乙烯（PE）拉伸基准上的性能差距，把差距 **按组件归因**（Pair / Bonded /
Neigh / Comm / Reduce），并区分 **算法/架构差距**（需重构）与 **实现差距**
（只需调优）。产出一份带数据的 **决策报告**，给出"重构 vs 优化"的明确建议和
排序后的工作流清单。后续的实现规范由本评估的结论决定，**不在此范围内**。

> 用户诉求（2026-06-04）："先评估我们架构是不是有问题，或者是不是比 LAMMPS
> 差很多，然后决定是重构还是优化。" 范围限定 **CPU only**。

## 背景数据 (已测)

- 全量 5 阶段墙钟：LAMMPS 32r=164s / 64r=86s（强扩展 1.9×）；Atomiverse
  32r=324s / 64r=336s（**32→64 不扩展**）。比值 av/lmp ≈ **2.0×（32r）/
  3.9×（64r）**。
- LAMMPS 自身在此规模已是 **通信主导**：64r 时 Comm=65% / Neigh=17% /
  Pair=11% / Bond=3%；负载严重不均（Nlocal min=4 / max=345 / ave=156）。
  32r 时 Comm=60% / Neigh=17% / Pair=16%。
- ⇒ 高 rank 差距很可能 **不在 kernel 速度**，而在 comm + 负载均衡 + 过度
  分解。必须用 **单 rank** 对比把"计算慢"和"通信/均衡慢"分开。
- 已有 `.claude/notes/cpu_md_perf_review_2026-05-24.md` 列出 CPU 热路径
  HIGH/MED 项（force reduction 5 趟、bonded 热循环里的 hash 查找、CellList
  重建 3 趟、NoseHoover/NPT 额外 Allreduce 等），本评估为这些项提供 **测量
  依据与优先级**，而非直接实现。
- 已修：64-rank NPT stage 转换死锁（`NoseHoover*::build` 里 `n_local>0`
  错误门控集合通信），见 NOTES。本评估在此修复之上进行。

## Domain basis (science.required)

所有计时跑 **固定物理**：使用与 `lammps_msu_paper.in` 一致的 5 阶段协议，从
同一 `msu_pe.data` 出发，并验证结果落在 `msu-pe-validation` 规范的容差内
（同物理比吞吐，绝不拿精度换速度）。

参考：
- LAMMPS 计时模型 —— per-section breakdown（Pair / Bond / Neigh / Comm /
  Modify / Output / Other）+ `%varavg` 负载不均指标。
- Plimpton, *J. Comput. Phys.* 117 (1995) 1 —— 空间分解的通信/计算标度律
  （comm ∝ surface, compute ∝ volume；小 atoms/rank 时 comm 主导）。
- 本仓 `.claude/notes/cpu_md_perf_review_2026-05-24.md`。

## Design — 四阶段，末端为决策门

**Phase A — 组件级计时插桩.**
给 CPU MD 加/扩 per-component 墙钟计时：Pair / Bonded（bond+angle+dihedral）
/ Neigh（cell-list rebuild + NL）/ Comm（halo forward + reverse + migration）
/ Reduce（ReduceEnergy + 任何 Allreduce）/ Other，加每阶段负载统计
（min/max/ave Nlocal、ghost 数）。输出 CSV，语义对齐 LAMMPS timing
breakdown，便于直接并排。复用已有的 per-stage 计时设施（msu-pe-validation
E.8）若存在；插桩为 **加法式、零物理改动**。

**Phase B — 标度矩阵.**
AV 在 **1 / 8 / 32 / 64 ranks** 跑 MSU PE（real `msu_pe.data`，可用
`--scale` 缩短但每 rank 数保持代表性），采集组件归因 + 负载不均。LAMMPS 参考
断面已在手（见背景数据；如需补 1r/8r 再跑 LAMMPS 对照）。

**Phase C — 单 rank kernel 隔离.**
1 rank（无 comm）下测 Pair / Bonded / Neigh 的纯计算吞吐（atom-steps/s 或
Matom-step/s），与 LAMMPS 单核同体系对比。判据：
- 若单 rank kernel 吞吐与 LAMMPS 在 ~1.5× 内、但标度差 → **实现/通信问题，
  走优化**（perf-review HIGH 项 + comm 批处理 + 负载均衡）。
- 若单 rank kernel 已 >2× 慢 → **实现差距在 kernel 本身**（向量化/数据布局/
  邻居算法），可能需要更深的重写。
- 若负载不均（max/ave Nlocal）主导 comm 等待 → **架构问题：缺 RCB/动态均衡**
  （LAMMPS `balance` / `comm_style tiled`），走重构。

**Phase D — 决策报告（门）.**
产出 `.claude/notes/cpu_md_perf_assessment_2026-06.md`：
1. 把 2.0×/3.9× 差距 **逐组件归因**（带 CSV 数据表）；
2. 每组件标注 **算法 vs 实现** 差距；
3. 定位 **32→64 不扩展** 的主因（负载不均 / RCB 缺失 / halo 开销 /
   per-step 重建）；
4. 映射到 `cpu_md_perf_review` 候选项并给 **预估收益**；
5. 给出 **重构 vs 优化的建议 + 排序工作流清单**，作为后续实现规范的输入。

## Files

- `src/cpu/**` + `include/atomiverse/cpu/**` — 仅计时插桩（加法式，零物理改动）。
- `benchmarks/md/` — per-component 计时基准 / `polyethylene_tension` 增加计时
  输出开关。
- `tools/` — 归因 CSV → 表格/图 的小脚本；标度矩阵提交脚本。
- `.claude/notes/cpu_md_perf_assessment_2026-06.md` — 决策报告（交付物）。

## Tasks

- [x] Add per-component CPU MD wall-clock timers (Pair/Bonded/Neigh/Comm/Reduce/Other) + load-balance stats, CSV output
- [x] Build scaling matrix: AV at 1/8/32/64 ranks on real msu_pe.data, capture attribution CSVs
- [x] Measure single-rank kernel throughput (atom-steps/s) for AV vs LAMMPS, per component
- [x] Write decision report: per-component attribution, algorithmic-vs-implementation classification, 32->64 root cause, ranked refactor-vs-optimize recommendation
- [x] Verify all timing runs stay within msu-pe-validation tolerances (same-physics guarantee)

## Testing

组件计时器单测（各组件计时之和 ≈ 总步时，±5%）；标度矩阵和单 rank 对比为
operator-action SLURM 运行（`srun --mpi=pmi2` 或 mpiexec）；决策报告人工评审。

## Out of scope

- 任何实际优化/重构的 **实现**（由本评估结论触发后续规范）。
- CUDA 后端性能（已有 `cuda-md-review-cleanup` 9 链）。
- 任何改变数值结果的改动；力场/积分器语义不变。
