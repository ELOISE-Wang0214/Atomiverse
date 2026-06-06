---
title: Multi-GPU bonded potentials + complete AV-vs-LAMMPS benchmark matrix
status: draft
created: 2026-06-06
---

# Multi-GPU bonded potentials + complete AV-vs-LAMMPS benchmark matrix

## Summary
两件事合为一个 spec，因为第二件依赖第一件：

1. **引擎(Phase 1):实现并优化跨卡(multi-GPU)版本的各类势** —— 当前 `MultiGpu` /
   `bench_multi_gpu_scaling` 只支持 LJ;`multi_gpu.cu` 完全没有 bonded 力。要测多卡
   PE,必须先让 `BondHarmonic / AngleHarmonic / DihedralMultiHarmonic`(以及
   special-bonds 排除)在单进程多卡域分解下正确运行:键合 ghost 拓扑随光环过去、各卡在
   owned+ghost 上算键力、ghost 力经 `cudaMemcpyPeer` reverse 回 owner,跨域 special
   排除生效。并对其做性能优化(减少跨卡传输、合并打包、与 LJ 力核共用 halo)。

2. **基准(Phase 2):完整测试矩阵** —— CPU/GPU × LJ/PE × 单卡多核/多卡 × 弱/强/体系
   扩展,**每个数据点跑 3 次**,保存每次的**原始数据**,再统一算 **mean ± std**。对手
   全程为 LAMMPS(CPU MPI / GPU KOKKOS-CUDA)。`dihedral multi/harmonic` 确认有
   `/kk` GPU 核,故 LAMMPS GPU-PE 是公平对比(此前关于二面角在 host 的担心作废)。

## Domain basis
不引入新物理(science.required=true,声明:无新方程)。两条可证伪基础:

- **跨卡键力正确性**:多卡 PE 的能量分解(E_pair / E_bond / E_angle / E_dihedral)与
  KE/virial,在相同初值下必须与单卡 PE 一致(浮点容差内),且与 CPU 参考一致。键力是
  N 体局域力,跨卡只改变"谁持有 owner、谁是 ghost",不改变物理 —— 与 CPU 多 rank 的
  `ResolveView`/`HaloRecv` 同构(见 `force_scratch.h` / `halo_exchange.cpp`)。
- **统计可靠性**:每点 3 次重复,报告 mean ± std。我们已经观测到节点争用会产生离群值
  (CPU 512r LAMMPS 85 s、GPU weak 4-card 54 s 的 straggler),单次测量不可信;3 次重复
  + std 让离群可见、可剔除(报告中说明剔除规则)。

## Design

### Phase 1 — 跨卡势(engine)
- **键合 ghost 拓扑**:扩展多卡 halo,使 ghost 原子携带其 bond/angle/dihedral 成员的
  全局 tag(类比 CPU `compute_bonded_extension` + tag→local 解析)。halo 截断取
  `max(pair_cutoff+skin, bonded_extension)`。
- **跨卡键力核**:每卡在 `[0, n_local+n_ghost)` 上跑 `cuda::Bond/Angle/Dihedral`,
  owner 过滤规则与单卡一致(bond=min-tag owner、angle=center、dihedral=2nd atom);
  partner 解析到本卡 ghost 行,解析失败计数(drop tally,门控诊断)。
- **reverse 力通信**:ghost 行上累加的键力经 `cudaMemcpyPeer`(NVLink)送回 owner 卡,
  与 LJ 的 cluster 力 reverse 合并到同一次 halo,避免二次同步(吸取 CPU "双 halo" 教训:
  reverse 恰好跑一次)。
- **special-bonds 跨域排除**:LJ 排除表的 partner tag 必须能解析到 ghost;排除表随
  halo/migrate 重建(类比 `pair_lj_cut.cpp::refresh_exclusion_map` 的 multi-rank 修复)。
- **优化**:bonded ghost 打包/解包合并为一次 P2P;键力核 `__launch_bounds__`;owner-grouped
  atomicAdd;与 cuda-md-review-cleanup 链的已落地优化(box by value、reductions)对齐。
- **bench 入口**:`bench_multi_gpu_scaling.cu` 增加 `--system pe` 路径(或新增
  `bench_multi_gpu_pe.cu`),用 `make_pe_melt_local` 等价的 GPU 端构建 + gen_pe 力场,
  发出 CSV(throughput + peak VRAM/卡)。

### Phase 2 — 完整测试矩阵(benchmark)
全部:每点 **3 次重复**,原始每次一行落 CSV,聚合脚本算 mean/std。

| # | 平台 | 体系 | 缩放 | 扫描点 | 备注 |
|---|------|------|------|--------|------|
| C1 | CPU 1–1024 核 | LJ | weak | 1,2,4,…,1024 核(16k/核) | 已有引擎 |
| C2 | CPU 1–1024 核 | LJ | strong | 固定 1M 原子 | 已有 |
| C3 | CPU 1–1024 核 | PE | weak | 1…1024 核(8k/核) | 已有 |
| C4 | CPU 1–1024 核 | PE | strong | 固定 ~1M 原子 | 已有 |
| G1 | 单 GPU | LJ | size | 32k→64M 原子 | 已有 |
| G2 | 1/2/4 GPU | LJ | weak | 16M/卡(总 16M/32M/64M) | 已有 |
| G3 | 1/2/4 GPU | LJ | strong | 固定 64M 总 | 已有 |
| **G4** | **单 GPU** | **PE** | size | 256k→8M(VRAM 上限) | `run_cuda_pipeline` 已有 |
| **G5** | **1/2/4 GPU** | **PE** | weak | 固定原子/卡 | **依赖 Phase 1** |
| **G6** | **1/2/4 GPU** | **PE** | strong | 固定总原子 | **依赖 Phase 1** |

- 每点两套:Atomiverse + LAMMPS,共享体系参数(§2 的 LJ / PE 力场,与 CPU 现状一致)。
- LAMMPS:CPU=MPI(`lj.in` / `lammps_bench.in`),GPU=KOKKOS-CUDA(`-k on g N -sf kk`),
  PE 全部 bonded 风格均有 `/kk`(已确认含 `dihedral multi/harmonic/kk`)。
- 步数:CPU 1000 步、GPU 100–500 步(warmup 后),reneigh 与现状一致(LJ 20 / PE 10)。

### 原始数据 + 统计
- **原始**:`results/<campaign>/raw.csv`,列:
  `code,system,platform,scaling,point,atoms,gpus_or_cores,rep,steps,wall_s,matom_step_s,mem_mb`
  (每点 3 行,`rep∈{1,2,3}`)。
- **聚合**:`aggregate.py` 读 raw.csv,按 (code,system,platform,scaling,point) 分组,
  输出 `summary.csv`:`…,matom_step_s_mean,matom_step_s_std,mem_mb_mean,mem_mb_std,n_rep,n_dropped`。
  离群剔除规则:wall_s > 1.5×median 的 rep 标记并排除(记 `n_dropped`),在报告中明示。
- 图:复用 `summary_2026-06-06/plot_summary.py` 的统一风格,误差棒 = std。

## Files to create or modify
- `src/cuda/md/multi_gpu.cu`, `include/atomiverse/cuda/md/multi_gpu.cuh` — 接入跨卡 bonded 力引擎 + 键合 ghost 拓扑 + reverse 力通信。
- `src/cuda/md/multigpu_bonded.cu` (可能新增) — 跨卡键力打包/解包/解析。
- `src/cuda/md/domain_decomposition.*`(多卡 DD)— ghost 携带键合 tag、跨域 tag→local。
- `src/cuda/potentials/{bond,angle,dihedral}_*.cu` — 若需 owned+ghost 范围参数化(对齐 CPU `ResolveView`)。
- `benchmarks/bench_multi_gpu_scaling.cu`(或新 `bench_multi_gpu_pe.cu`)— PE 多卡路径 + VRAM。
- `tests/cuda/md/test_cuda_multi_gpu_bonded.cu` (new) — 跨卡 PE 能量/力与单卡一致性。
- `benchmarks/md/results/<campaign>/` — `raw.csv`、`aggregate.py`、`summary.csv`、图。
- 全套 sbatch 脚本(CPU + GPU,3 reps)。

## Progress
- 2026-06-06: **reverse force halo 落地**(`MultiGpu::reverse_comm()` + `gather_ghost_forces`/`scatter_add_forces` 核 + `HaloPlan` 增 dst 侧 reverse buffer)。这是跨卡键力的前置:键力核在 ghost 行 atomicAdd 后,reverse_comm 按 z→y→x 逆 cascade 把 ghost 力 scatter-add 回 owner。LJ 路径不调用它(owned 全力,ghost 力恒 0),故无回归。待 GH200 编译 + `test_cuda_multi_gpu` 基线验证(job 9926)。

## Tasks
- [x] **P1.0** reverse force halo(`reverse_comm`)— foundational,待 GH200 编译验证
- [ ] **P1.1** 多卡 DD ghost 携带键合拓扑 + 跨域 tag→local 解析(RED: 跨卡键对解析单测)
- [ ] **P1.2** 跨卡 `BondHarmonic` 力 + reverse 通信;单卡↔双卡能量一致(test_cuda_multi_gpu_bonded)
- [ ] **P1.3** `AngleHarmonic`、`DihedralMultiHarmonic` 同上(逐个加,逐个验证)
- [ ] **P1.4** special-bonds 跨域排除;与 LJ 力核共享同一次 halo(reverse 恰跑一次)
- [ ] **P1.5** 优化:合并 P2P、`__launch_bounds__`、owner-grouped atomicAdd;profile 前后对比
- [ ] **P1.6** bench PE 多卡入口 + VRAM 报告
- [ ] **P2.1** 原始数据 schema + `aggregate.py`(mean/std/离群剔除)
- [ ] **P2.2** CPU 矩阵 C1–C4 全跑(3 reps,sbatch)
- [ ] **P2.3** GPU 矩阵 G1–G4 全跑(3 reps)
- [ ] **P2.4** GPU 矩阵 G5–G6(多卡 PE,依赖 P1)全跑(3 reps)
- [ ] **P2.5** 统一风格图(误差棒=std)+ README + summary.csv

## Testing
- **跨卡键力等价性**(新):单卡 PE vs 2/4 卡 PE,能量分解 + |F| 在浮点容差内一致;
  并与 `--backend cpu` 单 rank 参考一致 —— bonded 版的 `test_cuda_multi_gpu`。
- **reverse-once 不变量**:Σ F = 0(动量守恒)在多卡 PE 每步成立(防"双 halo"回归)。
- 已有 `test_cuda_multi_gpu`(LJ 13/13)保持绿。
- 统计:每点 3 reps,std/mean 在稳态点应 < ~5%(否则报告争用并重测)。

## Out of scope
- 新势(Coulomb/PME、EAM 等)的跨卡版本 —— 仅 bond/angle/dihedral + lj/cut。
- 非正交盒子的多卡 PE。
- CPU 引擎改动(CPU 多卡 PE 已工作)。
- 跨节点多 GPU(本 spec 限单节点 NVLink P2P,与现有多卡 LJ 一致)。
