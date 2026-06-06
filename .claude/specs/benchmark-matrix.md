---
title: Complete AV-vs-LAMMPS benchmark matrix (3 reps, mean±std)
status: approved
created: 2026-06-06
---

# Complete AV-vs-LAMMPS benchmark matrix (3 reps, mean±std)

## Summary
完整 benchmark 矩阵:CPU/GPU × LJ/PE × 单卡多核/多卡 × 弱/强/体系扩展,**每点 3 次重复**,保存每次原始数据,
统一算 **mean ± std**,与 LAMMPS(CPU MPI / GPU KOKKOS-CUDA)对比。多卡 PE(G5/G6)依赖 `multigpu-bonded` 链完成。
`dihedral multi/harmonic/kk` 确认存在 → LAMMPS GPU-PE 公平。

## Domain basis
统计可靠性:已观测节点争用产生离群(CPU 512r 85 s、GPU weak 4-card 54 s)。3 次重复 + std 让离群可见可剔除
(`wall_s > 1.5×median` 标记并排除,记 n_dropped),报告明示。体系参数沿用既有 LJ(ρ*=0.712,rc=2σ,skin 0.6,
reneigh 20)与 PE(gen_pe Dreiding-UA,20-chain,bond/angle/dihedral+lj/cut)。

## Design
矩阵:C1 LJ-CPU-weak、C2 LJ-CPU-strong、C3 PE-CPU-weak、C4 PE-CPU-strong、G1 LJ-GPU-single-size、
G2 LJ-GPU-multi-weak、G3 LJ-GPU-multi-strong、G4 PE-GPU-single-size、G5 PE-GPU-multi-weak、G6 PE-GPU-multi-strong。
每点 AV + LAMMPS × 3 reps。raw.csv:`code,system,platform,scaling,point,atoms,gpus_or_cores,rep,steps,wall_s,matom_step_s,mem_mb`。
`aggregate.py` 分组 → summary.csv:`*_mean,*_std,n_rep,n_dropped`,离群规则 `wall>1.5×median`。图复用统一风格,误差棒=std。
内存:CPU `/usr/bin/time -v` max RSS/rank(两码);GPU 峰值显存。

## Files to create or modify
- `benchmarks/md/results/matrix_<date>/` — sbatch 全套(CPU+GPU,3 reps)、`raw.csv`、`aggregate.py`、`summary.csv`、图、README。

## Tasks
- [ ] **PB.1** raw-data schema + aggregate.py (mean/std + wall>1.5×median outlier rule, n_dropped)
- [ ] **PB.2** CPU matrix C1–C4, 3 reps each (sbatch, EPYC)
- [ ] **PB.3** GPU matrix G1–G4, 3 reps each (GH200)
- [ ] **PB.4** multi-GPU PE matrix G5–G6, 3 reps (depends on multigpu-bonded chain)
- [ ] **PB.5** unified-style figure (error bars = std) + README documenting method

## Testing
- aggregate.py unit-checked on synthetic raw rows (mean/std/outlier); steady-state std/mean < ~5% or contention noted.

## Out of scope
- Engine changes (the multigpu-bonded chain owns multi-GPU PE); new potentials; cross-node multi-GPU.
