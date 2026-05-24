---
title: CUDA MD Review Cleanup — 总索引 (reference)
slug: cuda-md-review-cleanup
status: reference
created: 2026-05-24
axes: arch, perf
findings_count: 47
severity_distribution: 4 🚨 / 16 🔴 / 21 🟡 / 6 🟢
parent_review: 2026-05-24 /mol:review --axis=perf,arch on cuda-md scope
chain: cuda-md-review-cleanup-01..09
---

# CUDA MD Review Cleanup — 总索引 (reference, not implementable)

## 用途

本文件 **不是** 可执行 spec — `status: reference`,无 `<slug>.acceptance.md`,`/mol:impl` 会拒绝。它是 2026-05-24 `/mol:review --axis=perf,arch` 在 CUDA MD scope 下产出的 **47 项发现的完整清单 + 修复路线图**,作为下面 9 个可执行 sub-spec 的总索引。

每个 sub-spec 单独可 `/mol:impl`,但都按 cluster 编号引用本文件。当 `cuda-md-review-cleanup-09` 的 umbrella `/mol:review` 回归 gate 通过后,本 inventory 可改为 `status: done` 或删除。

## Scope

```
include/atomiverse/cuda/md/          ← 所有 .cuh
src/cuda/md/                          ← 所有 .cu(含 neighborlist/ 子目录)
include/atomiverse/cuda/neighbor_list.cuh + src/cuda/neighbor_list.cu
include/atomiverse/cuda/potentials/   ← 所有 .cuh
src/cuda/potentials/                  ← 所有 .cu
include/atomiverse/cuda/observer/stress.cuh + src/cuda/observer/stress.cu
include/atomiverse/cuda/{cuda_driver,dynamics_driver,component,context,frame,box,async_bridge,launch}.cuh
src/cuda/{driver,async_bridge}.cu
benchmarks/bench_md_step.cu
```

明确 **out-of-scope**(不含):SCF / `cuda::scf::*`、minimizer / `cuda::minimizer::*`、recorder I/O 格式、transfer/materialize、多 GPU halo-exchange(cluster 13,仅 NOTES 跟踪),以及 user MEMORY 已记录的三项 correctness 修复(`[[cuda-driver-step-counter-parity]]` / `[[cuda-nbnxn-smart-prune-missing]]` / `[[atv-cuda-angle-theta0-double-deg2rad]]` 本身,但其文档/测试 follow-up 在 cluster 11 中)。

## 9-link chain 概览

| Link | Slug | Cluster | 发现 | 严重度 | 关键交付 |
|---|---|---|---|---|---|
| 01 | `cuda-md-review-cleanup-01-component-capability` | C1 | 5 | 3 🔴 / 2 🟡 | `ComponentBase` `reads()`+`needs_neighbor_list()` 虚函数;`install_infra_` 类型化重写 |
| 02 | `cuda-md-review-cleanup-02-nlist-skin-gating` | C3+C4 | 8 | 1 🚨 / 3 🔴 / 4 🟡 | Verlet-skin 门控 + 邻居列表过启动修复(预计 ~28× 重建加速) |
| 03 | `cuda-md-review-cleanup-03-launch-purity` | C2 | 8 | 3 🚨 / 3 🔴 / 2 🟡 | hot-path 零 D2H/零 sync(NPT/LJ/Driver/AsyncBridge/ZeroMomentum) |
| 04 | `cuda-md-review-cleanup-04-box-by-value-pbc` | C6 | 5 | 1 🔴 / 4 🟡 | `cuda::Box` 按值传递 + `bool ortho` 模板特化 |
| 05 | `cuda-md-review-cleanup-05-reductions-and-atomics` | C5+C7 | 7 | 4 🔴 / 3 🟡 | warp-shuffle 规约 + cluster-staged 力写回 |
| 06 | `cuda-md-review-cleanup-06-launch-config-streams` | C8 | 5 | 4 🔴 / 1 🟡 | `__launch_bounds__` + 力流池(2-4× force-kernel 提速) |
| 07 | `cuda-md-review-cleanup-07-kernel-header-and-microfixes` | C9+C12 | 8 | 2 🔴 / 5 🟡 / 1 🟢 | 共享 kernel 头签名对齐 + 6 项 micro-fixes |
| 08 | `cuda-md-review-cleanup-08-placement-and-conventions` | C10+C11 | 7 | 4 🟡 / 3 🟢 | 文件重定位 + docstring/alias header 清理 |
| 09 | `cuda-md-review-cleanup-09-bench-and-regression` | C14+umbrella | 1 + ac-015 | 1 🟡 | PE-NPT benchmark + `/mol:review` 回归 gate |

**总计:47 项发现,9 个可执行 sub-spec,1 个 umbrella regression gate。**

## 链节依赖与执行顺序

```
01 (capability) ─┐
                 ├─→ 03 (purity) ──┐
02 (skin) ───────┘                 │
                                   ├─→ 09 (regression)
04 (Box) → 05 (reduce) → 06 (launch) → 07 (header+micro) → 08 (placement) ─┘
```

- **01 必须在 03 前**:03 的零-sync 测试需要 01 引入的类型化 `reads()` 才能让 `install_infra_` 正确注入 VirialStress,否则 NPT 跑不起来。
- **02 必须在 03 前**:不门控邻居列表重建,每步 radix sort 是不可消除的 sync barrier,03 的零-sync 断言无法通过。
- **04→05→06 顺序而非依赖**:都触碰相同的 hot kernel(LJ pair / bonded),串行避免 merge 冲突。
- **07 = 低风险扫尾**:kernel header 签名对齐 + 6 项 micro-fixes,放在 hot kernel 稳定后。
- **08 必须放最后**:重定位 `src/cuda/md/neighborlist/` 和 `cuda/observer/stress.{cuh,cu}` 会和 01-07 的 diff 全部冲突。
- **09 = umbrella gate**:跑 `/mol:review --axis=perf,arch` 并断言 0 🚨 / ≤3 🔴 / ≤8 🟡 残留;PE-NPT benchmark 记录每个 cluster 的 wall-clock delta。

---

## Cluster 详表(47 项发现完整记录)

每项含:`严重度 file:line — 症状 (1 句) | 修复 (1 句)`。引用源:2026-05-24 `/mol:review` aggregated report。

---

### Cluster 1 — `ComponentBase` capability-typed auto-injection(结构根)

**根因**:`cuda::ComponentBase`(`include/atomiverse/cuda/component.cuh`)缺 `reads()` 和 `needs_neighbor_list()` 虚函数;`cuda::md::MD::install_infra_`(`include/atomiverse/cuda/md/md.cuh:243-256`)退化为 `name()=="virial_stress"` / `name()=="neighbor_list"` 字符串代理 + `dynamic_cast<const cuda::Pair*>`。CPU 端契约用 typed `reads()` —— CUDA 端必须对齐。修复后多个下游 🟡 自动消失。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `include/atomiverse/cuda/observer/stress.cuh:22` | `cuda::VirialStress` `role()` 和 `capabilities()` 都默认 None;CPU 同胞类正确声明 `Observer`+`Stress` | 重载 `role()→Observer`、`capabilities()→Stress` |
| 🔴 | `include/atomiverse/cuda/potentials/angle.cuh:19-21` | 声明 `Energy\|Force` 但 `src/cuda/potentials/angle_harmonic.cu:145-200` 写入 `f.scalars->stress[]` | `capabilities()` 加 `\| Stress` |
| 🔴 | `include/atomiverse/cuda/potentials/dihedral.cuh:19-21` | 声明 `Energy\|Force` 但 `src/cuda/potentials/dihedral_multi_harmonic.cu:213-282` 写入 Voigt stress | `capabilities()` 加 `\| Stress` |
| 🟡 | `include/atomiverse/cuda/async_bridge.cuh:39` | 基类无 `role()`;CUDA `Thermo`/`RDF`/`Mulliken`/`Recorder`/`DensityRecorder` 全部继承 `role()==None` | 每个具体子类各自声明 `role()` |
| 🟡 | `include/atomiverse/cuda/neighbor_list.cuh:48` + `include/atomiverse/cuda/md/nose_hoover_npt.cuh:84` + `include/atomiverse/cuda/md/md.cuh:243-256` | 结构根缺 `reads()`/`needs_neighbor_list()`;`NeighborList` 无 role/capabilities;`NoseHooverNPT` 兼任 Integrator+Thermostat+Barostat 却无 `reads()==Stress` | `ComponentBase` 加两个虚函数 → `NeighborList`/`NoseHooverNPT` 声明 → `install_infra_` 用 capability 位匹配,弃用字符串/`dynamic_cast` |

---

### Cluster 2 — `launch()` GPU-purity(零 D2H / 零 sync)

**根因**:`CLAUDE.md ## Architecture` § "CUDA execution model" 与 `component.cuh:122-129` 规定 `launch()` 必须 GPU-pure(零 D2H、零 `cudaStreamSynchronize`)。前序 spec `msu-pe-cuda-components-07-async-correctness` 已落地 `cuda::stream_synchronize` 计数器,但测试循环只覆盖 LJ-pair + VelocityVerlet + `CudaCsvRecorder` —— **未覆盖 NPT / 邻居列表重建 / `ZeroMomentum`**。下列违例必须修复。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🚨 | `src/cuda/potentials/pair_lj_cut.cu:298-301` | `compute()` 每步 `cudaMemcpyAsync` D→H `n_pairs` + 阻塞 `cudaStreamSynchronize` 来定 grid 尺寸 | 用 `grid_size(pair_capacity)` 上界启动;kernel 已有 `if (idx < n_pairs)` 早退 |
| 🚨 | `src/cuda/md/nose_hoover_npt.cu:248,266,505-512,540-545,598-602` | NPT 每步 3-5 次 sync + ~12 次小 D2H(KE 规约/box 对角线读回/extended_hamiltonian/end_of_step) | 单 device 规约 kernel → 8-double pinned buffer + `cudaMemcpyAsync` + event polling;box 对角线用 device 版本计数器跟踪 |
| 🚨 | `src/cuda/driver.cu:71-72` | 每相位双流阻塞 sync(6 phase × 2 stream = 12 syncs/step ≈ 80 µs/step on Hopper),串行化所有相位 | 删除 per-slot sync;只在 `run()` 结束或 CPU-bridge 组件需要时 sync |
| 🔴 | `src/cuda/md/init.cu:233` | `ZeroMomentum::launch()` 每 `freq_` 步 D2H+sync 用于 COM 移除 | 在 device 上算 v_cm 写到 resident scratch;`remove_com` 直接读;或正式把 "integrity component" 加入 NOTES 分类表 |
| 🔴 | `src/cuda/md/nose_hoover_npt.cu:300-334, 336-382` | NHC chain 在 host 跑 `std::exp` + 交错 kernel launch(mtchain=3 → ~30-50 µs/step) | 融合为单 "scale-velocity-N-times" kernel,接 host 预算的 `eta_t_dot[]` 数组 |
| 🔴 | `src/cuda/md/nose_hoover_npt.cu:540-545, 598-602` | `extended_hamiltonian_` 2 次 D2H + `end_of_step` 2 次 H2D | 整体上 device 算 H,或合并到一次 pinned + event sync |
| 🟡 | `src/cuda/md/nose_hoover_npt.cu:399-404` | `update_eps_dot_` 每次 D2H 6 个 stress double(NPT 每步两次) | 单 pinned buffer + event polling |
| 🟡 | `src/cuda/async_bridge.cu:39` | `cudaStreamSynchronize(defer_stream)` 在 `on_data_ready` 之前,串行化 CPU recorder 与下游 device 工作 | 高频 recorder 改用 `cudaLaunchHostFunc` / stream callback |

**验证仪表**:扩展 `tests/cuda/async/test_hot_path_no_sync.cu` 驱动 NPT 循环(NoseHooverNPT + NeighborList + LJ + 全 bonded + VirialStress),断言 `n_stream_syncs_during_md_step_ == 0`。

---

### Cluster 3 — NeighborList Verlet-skin 重建门控(🚨 单项最大 perf 收益)

**根因**:`src/cuda/neighbor_list.cu:177-320` 每步重建整条管线。死代码 `launch_cluster_displacement`(`src/cuda/md/neighborlist/fused_search.cu:116-136`)+ `ref_cluster_{x,y,z}` + `rebuild_flag` buffers(`src/cuda/md/neighborlist/nlist_impl.cuh:54,68`)已分配但 **从未被读写**。接好后约 28× 重建成本加速。关联 user MEMORY `[[cuda-nbnxn-smart-prune-missing]]`。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🚨 | `src/cuda/neighbor_list.cu:177-320` + `src/cuda/md/neighborlist/fused_search.cu:116-136` | 整条管线(Morton + radix sort + AABB + grid + count + fill + compact + csr_to_flat + prune + export)每步重跑;~30-50% 步进时间浪费 | 重建后快照 cluster 中心到 `ref_cluster_{x,y,z}`;接通 `launch_cluster_displacement` → `rebuild_flag`;管线入口 gate `rebuild_flag != 0`(1 步滞后可接受) |
| 🔴 | `src/cuda/md/neighborlist/sort_and_cluster.cu:107-123` | 每步全 CUB radix sort(~30-50 µs/call) | 与 skin gating 共同生效 |
| 🔴 | `src/cuda/neighbor_list.cu:267-269` | `cudaMemsetAsync(interaction_mask, 0xFF, pair_capacity*8)` 每步(5M-pair cap ≈ 40 MB,HBM3 ~10 µs) | 移到 `build()`,或只 memset csr_to_flat 实际写入的范围 |

---

### Cluster 4 — NeighborList 过启动 / 单线程 kernel 浪费

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `src/cuda/neighbor_list.cu:190-194` | `launch_setup_cell_device<<<1,1>>>` 每步运行;NVE/NVT box 固定 | host 缓存 dims/n_stencil;只在 box volume 变化超阈值时重启 |
| 🔴 | `src/cuda/md/neighborlist/fused_search.cu:233-235` + `src/cuda/neighbor_list.cu:312` | `smart_rolling_prune` 和 `export_atom_pairs_kernel` 用 `pair_capacity` 定 grid;~3× 实际 pair 数 → ~67% 线程浪费 | 滞后读 device-side `flat_pair_counter` 给下一轮 grid sizing,或 config-write + cooperative-group |
| 🟡 | `src/cuda/md/neighborlist/setup.cu:178,193` | `setup_cell_from_*<<<1,1>>>` 单线程 kernel,O(stencil³) 三重循环 + O(N²) 周期 dedup | host 在 `build()` 算一次,`cudaMemcpyToSymbolAsync` 上传 `__constant__`;只在 box 变化时重传 |
| 🟡 | `src/cuda/neighbor_list.cu:233, 293` | `cut_skin_sq` / `cutoff_sq` 每步重算 | 缓存在 `State`,`build()` 时算一次 |
| 🟡 | `src/cuda/md/neighborlist/sort_and_cluster.cu:112-120` | CUB `temp.resize(needed)` 首次可能 `cudaMalloc` | 在 `build()` 预分配到最大尺寸 |

---

### Cluster 5 — Stress / energy 规约模式(warp-shuffle 迁移)

**根因**:hot kernel 用 6 次顺序 shared-mem 规约 + 完整 `__syncthreads()` 算 stress tensor;sm_90 有硬件加速的 `__reduce_add_sync` / `cooperative_groups::reduce`,一个 barrier 即可。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `src/cuda/potentials/pair_lj_cut.cu:140-155` + `bond_harmonic.cu:97-110` + `angle_harmonic.cu:145-158` + `dihedral_multi_harmonic.cu:237-250` | 4 个 force kernel 各做 6 次顺序 shared-mem stress 规约 + 完整 barrier(累计 ~4-8 µs/step) | 7 个规约(1 E + 6 stress)合并到单次 `__reduce_add_sync` 或 `cooperative_groups::reduce` |
| 🔴 | `src/cuda/observer/stress.cu:103-105` | `kinetic_stress_k` 独立 launch(NPT-only,~5 µs launch + 12 KiB shared) | 与 `KineticEnergy` 融合 |
| 🟡 | `src/cuda/observer/stress.cu:43-66` | `kinetic_stress_k` 用 12 KiB shared/block(限制 4-5 blocks/SM) | warp-shuffle (`__shfl_xor_sync`)→ ~384 bytes |
| 🟡 | `src/cuda/md/nose_hoover_npt.cu:129-167` | `reduce_ke_axes_k` 用 6 KiB shared/block | warp-shuffle 3-axis 规约 |
| 🟡 | (多处 KE/stress 规约) | `if (s_energy[0] != 0) atomicAdd(...)` 代码味 | 删除条件;`atomicAdd(0)` 在现代 GPU 上免费;或迁移到 `__reduce_add_sync` |

---

### Cluster 6 — `cuda::Box` 按值传递(PBC kernel 改造)

**根因**:integrator kernel 各自 inline 重写 PBC wrap:每线程算 `inv_Lx = 1/Lx`(FP64 div ~20-30 cycle × 3 轴 × N atoms)+ 读 `box[0],box[4],box[8]` 3 个间隔很远的 uncoalesced load,尽管 `cuda::Box`(`include/atomiverse/cuda/box.cuh:31-91`)已经缓存了 `inv_l/half_l/ortho`。bonded kernel 也用 `round(dx/bx)` 而非缓存。triclinic `minimum_image` 在 ortho=true 时仍走完整路径。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `src/cuda/md/velocity_verlet.cu:42-47` + `nve_limit.cu:62-67` + `nose_hoover_npt.cu:101-110` (合并 perf B + perf H + arch §8) | per-thread PBC 算 `1/Lx` FP64 div + uncoalesced `box[0/4/8]` 读 + 跨 kernel 算法重复 | 改为按值传 `cuda::Box`,直接用缓存的 `inv_l[]` |
| 🟡 | `include/atomiverse/cuda/box.cuh:31-91` | `cuda::Box` 有缓存但 kernel 不用;triclinic 始终走完整路径 | 整体迁移 kernel 至按值传 `cuda::Box`(`box.cuh:8-14` 已注 defer) |
| 🟡 | `src/cuda/potentials/bond_harmonic.cu:53-55` + `angle_harmonic.cu:62-67` + `dihedral_multi_harmonic.cu:72-80` | bonded kernel 用 `round(dx/bx)` 而非缓存 `inv_l` | `dx - bx * nearbyint(dx * inv_lx)` 或采用 `cuda::Box::min_image` |
| 🟡 | `src/cuda/md/neighborlist/fused_search.cu:186-203` | `smart_rolling_prune` 内循环对 ortho box 也跑完整 triclinic(4-5× 工作量,每 cluster pair 64×) | 按 `bool ortho` 模板特化;同样应用于 `export_atom_pairs_kernel`、`fill_pairs`、`min_aabb_dist_sq` |
| 🟡 | `src/cuda/md/neighborlist/pair_search.cu:124-214` | `fill_pairs` warp-per-cluster 内 125 stencil 串行 | 按 ortho 模板化,或 persistent CTA + `cooperative_groups` 负载均衡 |

---

### Cluster 7 — 力 atomic 争用(cluster-staged 写回)

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `src/cuda/potentials/pair_lj_cut.cu:102-107` | 每 pair 6 个 atomicAdd(i 3 个 + j 3 个);稠密 PE → ~1000 atomicAdds/atom × 10k atoms ≈ 10-30% step time on Hopper | cluster-staged 写回 — shared-mem 内累加分量力,每 block 每 atom 1 个 atomicAdd(~64× 争用降低) |
| 🔴 | `src/cuda/potentials/bond_harmonic.cu:69-74` + `angle_harmonic.cu:111-119` + `dihedral_multi_harmonic.cu:197-208` | N 个 atomicAdd/topology(bond 9、angle 9、dihedral 12)→ ~2M/step ≈ 60 µs/step | 按 dest atom 分组(LAMMPS per-atom owner list);寄存器累加,最后写一次 |

---

### Cluster 8 — Launch 配置 / 寄存器压力 / 多流

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `src/cuda/potentials/pair_lj_cut.cu:43-156` + `kinetic_energy.cu` + `reduce_ke_and_update` + `bond_harmonic.cu` + `angle_harmonic.cu` + `dihedral_multi_harmonic.cu` + `observer/stress.cu` (`kinetic_stress_k`) | hot kernel 无 `__launch_bounds__`,编译器过度保守地分配寄存器,占用率受限 | 加 `__launch_bounds__(256, 4)`;用 `nvcc --ptxas-options=-v` 测量 |
| 🔴 | `src/cuda/potentials/dihedral_multi_harmonic.cu:30-251` | ~30+ double 局部变量,寄存器压力逼近 255-reg 上限 → 1-2 blocks/SM | `__launch_bounds__(128, 4)` 强制 spill;若仍欠占用则分解 kernel |
| 🔴 | `src/cuda/md/velocity_verlet.cu:67-82` + `nve_limit.cu:94-113` | 两阶 Verlet 每步 2 次 launch(~10 µs);小 N 时 launch 占主 | (可选)小 N(<10k)时用 cooperative-groups grid-sync 一融合 kernel(复杂度高) |
| 🔴 | `src/cuda/driver.cu:54-78` | 所有组件都在 `context_.stream` / `defer_stream`,无 per-component 隔离;bonded + LJ 独立 compute 但串行,SM 空闲 | 力流池 — Pair/Bond/Angle/Dihedral 各一个流;`cudaEventRecord` force-start + `cudaStreamWaitEvent` fan-in;预期 2-4× force-kernel 提速 |
| 🟡 | `src/cuda/potentials/pair_lj_cut.cu:92` | `inv_r = 1/r` FP64 div + pair search 已算 sqrt(r²)≈ 100 µs/step | `export_atom_pairs` 存 r²(不开根),LJ 用 `inv_r2 = 1/r2`(单 rcp);`__drcp_rn` / `__frcp_rn` |

---

### Cluster 9 — 共享 kernel 头 / 定义签名分歧

**根因**:`*_kernel.cuh` 头声称"single source of truth"但签名与 `.cu` 定义不符;torch 后端调用 phantom 符号。

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🔴 | `include/atomiverse/cuda/potentials/lj_kernel.cuh:32-41` | 头 14 参;`src/cuda/potentials/pair_lj_cut.cu:43-55` 18 参(加 stress + special-bonds CSR);torch 后端 `src/torch/potentials/lj.cu:96-105` 调 phantom 14 参 | 定义 14 参 overload 转发到 18 参,或更新头 + torch 至 18 参 |
| 🔴 | `include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh:29-37` | 同样的"共享 kernel"头分歧模式 | 校验 torch 调用点,同样修复 |

---

### Cluster 10 — 目录 / 文件归置

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🟡 | `src/cuda/md/neighborlist/`(整个目录) | 私有 nlist impl 在 MD-domain 文件夹下,但公共组件在 `cuda/` root(跨域) | 整体迁移至 `src/cuda/neighborlist/` |
| 🟡 | `include/atomiverse/cuda/observer/stress.cuh` + `src/cuda/observer/stress.cu` | `cuda::VirialStress` 在 `cuda/observer/` 下,而 CPU 同胞类是 `cpu/md/virial_stress.h`;CLAUDE.md 把 VirialStress 列在 MD 组件 | 重定位至 `cuda/md/virial_stress.{cuh,cu}` |

---

### Cluster 11 — 惯例 / 文档 / 蓝图刷新

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🟡 | `include/atomiverse/cuda/potentials/angle_harmonic.cuh:18` | docstring 说 `[degrees] (内部转弧度)` 但实现(`angle_harmonic.cu:167-174`)说 `theta0_ 已经是弧度,API 契约`(MEMORY `[[atv-cuda-angle-theta0-double-deg2rad]]`) | 更新头 docstring 为 `[rad]` |
| 🟡 | `include/atomiverse/cuda/md/md.cuh:54-55` | 多 GPU halo-exchange 在 CUDA 端 out-of-scope(无 `#ifdef ATV_BUILD_MPI`),与 CPU 契约真实分歧 | 加 NOTES 跟踪条目(实现不在本 spec) |
| 🟢 | `.claude/architecture.md:91` | 蓝图称 `cuda::Frame` 嵌入 `SCFWorkspace`;当前 `frame.cuh:60-137` 不再如此(纯文档分歧) | 重跑 `/mol:map` |
| 🟢 | `include/atomiverse/cuda/recorder/thermo.cuh:25` | `cuda::Thermo` 在 `ATV::cuda::`(非 `ATV::cuda::recorder::`);CPU 同胞类一致 | 确认这是有意模式;无代码改动 |
| 🟢 | `include/atomiverse/cuda/potentials/lj.cuh:15` + `harmonic_bond.cuh:13` + `harmonic_angle.cuh:10` + `multi_harmonic_dihedral.cuh:10` | 2026-05-19 LAMMPS-命名 refactor 后保留的后向兼容 alias headers;CLAUDE.md "no backwards-compat shims" 规则 | 下次 breaking-change PR 中删除 |

---

### Cluster 12 — Mini-fixes(低影响)

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🟡 | `src/cuda/md/fix_langevin.cu:43-49` | `box_muller` 每 Gaussian pair 用 cos/sin/sqrt/log(4 transcendental,~50 µs/step @ 100k atoms) | `sincos(theta, &s, &c)` 单 PTX op;考虑 `sincospi` |
| 🟡 | `src/cuda/md/barostat.cu:40-73` | `berendsen_scale_k` 每个 block 的 thread 0 重算 `cbrt(mu)`(100k atoms = 400 个 block × redundant cbrt) | 在 host 或 `<<<1,1>>>` 预算 μ,作为 kernel arg 传入 |
| 🟡 | `src/cuda/potentials/bond_harmonic.cu:151` + `angle_harmonic.cu:196` + `dihedral_multi_harmonic.cu:277` | 每步 3 个独立 bonded launch(~15 µs);无法真正融合 | (可选)下游 E/stress gather kernel |
| 🟡 | `src/cuda/potentials/dihedral_multi_harmonic.cu:115-119` | `dv_dcos` 每 dihedral 重算 `2*C[2], 3*C[3], 4*C[4]` | 预算系数数组存储 |
| 🟡 | `src/cuda/md/neighborlist/triclinic.cuh:47-62` | `invert_3x3` device 端每步重算(NVE/NVT 常量) | 只在 box 变化时重算 |
| 🟢 | `src/cuda/md/deform.cu:113` | `if (i==0) { box[0] = ...; }` 在 `propagate_x_k` 内分裂 hot path | 独立 `<<<1,1>>>` follow-up 或并入 box-update kernel |

---

### Cluster 13 — 多 GPU halo-exchange 平价(**OUT OF SCOPE**,仅 NOTES 跟踪)

CPU 端在 `set_pair()` 调用时 `#ifdef ATV_BUILD_MPI` 自动注入 halo exchange;CUDA 端不支持。User 已明确本 spec 不实现此项;cluster 11 中将其作为 NOTES 条目跟踪。

---

### Cluster 14 — Benchmark 覆盖

| Sev | File:line | 症状 | 修复 |
|---|---|---|---|
| 🟡 | `benchmarks/bench_md_step.cu` | 只测 N=100/1000/10000,cutoff=3.0,density=0.8;缺稠密 PE-scale + NPT(🚨/🔴 集中之处) | 加 PE-realistic 用例 — N=10000-50000,cutoff 10.5 Å + skin 0.3,完整 bonds+angles+dihedrals 拓扑,`cuda::NoseHooverNPT` |

---

## Umbrella 回归 gate(在 link 09 中实现)

- 跑 `/mol:review --axis=perf,arch` 在本 inventory 同样 scope 上
- 断言残留:**0 🚨 / ≤3 🔴 / ≤8 🟡**
- 残留必须显式枚举在 link 09 的 `acceptance.md` 中作为已知 backlog,否则失败
- 在 GH200 PDC 节点记录 baseline-vs-after wall-clock,每 cluster 至少:
  - cluster 2 ≥ 1.20× 加速
  - cluster 3 ≥ 1.30× 加速
  - cluster 5+6+7 合并 ≥ 1.10× 加速
  - 其他 ≥ 0.95×(不允许回归)

## 已采纳的可重用规则建议(在 link 08 中落地到 CLAUDE.md / NOTES)

1. `ComponentBase` 必须暴露 typed `reads()` / `needs_neighbor_list()`;auto-injection **不得** 用 `name()` 字符串代理。
2. `launch()` 是 GPU-pure(零 D2H、零 `cudaStreamSynchronize`);例外仅限"integrity component",必须在 `.claude/NOTES.md` 显式分类。
3. `NeighborList::launch()` 必须用 Verlet-skin / cluster-displacement gate 全重建;每步重建是 hot-path bug。
4. 共享 kernel 头(`*_kernel.cuh`)声明的"single source of truth"签名 **必须** 与 `.cu` 定义完全一致,否则必须提供转发 overload。
5. 每 kernel 重写 PBC wrap 是反模式;按值传 `cuda::Box` 或传 `inv_Lx/Ly/Lz` 作为 kernel arg。

## 关联资源

- 触发 review:本对话上下文中 2026-05-24 `/mol:review --axis=perf,arch` 报告(已 archive 到本 inventory)
- 前序 spec:`msu-pe-cuda-components-07-async-correctness`(已 code-complete,提供 `cuda::stream_synchronize` + counter 仪表,本 chain 复用)
- User MEMORY:`[[cuda-driver-step-counter-parity]]` / `[[cuda-nbnxn-smart-prune-missing]]` / `[[atv-cuda-angle-theta0-double-deg2rad]]`(correctness,本 spec 仅作为 cluster 11 / cluster 3 的文档/接线 follow-up)
- 架构契约:`CLAUDE.md ## Architecture` + `## Conventions`
- 项目蓝图:`.claude/architecture.md`(由 cluster 11 的 `/mol:map` 重跑刷新)
