---
title: arch-refactor-05-context — Context 分级：按 driver 域嵌套 service 子结构体（cpu/cuda 同构）
status: approved
date: 2026-07-10
chain: arch-refactor (05/09)
depends_on: [arch-refactor-03-neighbor, arch-refactor-04-engines]
---

# Context 分级：按 driver 域嵌套 service 子结构体（cpu/cuda 同构）

## Summary

`cpu::Context`（include/atomiverse/cpu/component/context.h:48-128）当前把核心执行配置（n_threads / mpi_rank / mpi_size / timer / outer_status）与四个互斥 driver 域的 service 指针（MD/DD：dd、pair_list、is_last_step；MC：proposal、rng；Optimize：reeval_forces、line_search_trial；Analysis：executor）平铺在一个结构体里，每加一个 driver 类别就在全体组件可见的公共面上多一批"别的 workflow 才用"的 nullable 字段。本 spec 按用户钦定形态把这些域字段重组为**嵌套 service 子结构体**：`Context{ …核心字段平铺…; DynamicsServices md; SamplingServices mc; OptimizeServices opt; AnalysisServices analysis; }`，组件改为 `ctx.md.pair_list` / `ctx.mc.proposal` / `ctx.opt.reeval_forces` / `ctx.analysis.executor` 访问；`cuda::Context` 按同构分组对齐（core 平铺 + `md` slice）。这是零成本分组：单结构体、纯聚合、栈上构造、常量偏移取址全部不变，不引入任何指针间接层；全部读点 clean-break 一次 commit 迁移，无兼容 shim。"未填 slice 的字段保持 nullptr = 服务不可用、组件用前判空"的既有语义逐字保留。

## Domain basis

本 spec 不引入新物理、不触及任何数值路径——它是纯字段重组，约束是**行为逐位不变**：同种子、同配置下 MD/MC/Minimize 轨迹与能量/力逐位相同。被守护的既有物理契约（作为回归锚，不新增公式）：

- Metropolis 细致平衡与单 RNG 流可复现性（`ctx.mc.rng` / `ctx.mc.proposal`；物理门禁沿用 `.claude/specs/mc-sampling-03-validate.acceptance.md` 的 P(r)∝r²e^{−βV} / NIST SRSW 锚，Metropolis et al. 1953, DOI 10.1063/1.1699114）。
- DD/halo 力完整性（`ctx.md.dd` / `ctx.md.pair_list`；owned/ghost 分解下力与能量与单 rank 一致）。
- Line-search 试探几何能量一致性（`ctx.opt.reeval_forces` 重入求值 + `ctx.opt.line_search_trial` 抑制迁移，LBFGS/CG 收敛行为不变）。

单位约定不变（eV, Å, fs, amu, K；core/units.h）。

## Design

### 新结构（cpu）

`include/atomiverse/cpu/component/context.h` 内定义四个纯聚合子结构体，`Context` 组合之；核心字段**保持平铺**（钦定形态允许"Core core 或平铺核心字段"，取平铺：n_threads/mpi_rank/mpi_size/timer/outer_status 的读点遍布全体组件与 SCF/PW 域，远多于域字段读点，平铺使本 spec 迁移面收敛到域字段本身，且 `outer_status` 契约（facade/iteration_status.h）字面不动）：

```cpp
struct DynamicsServices {   // 经典力管线 service（发布方：Backend / NeighborList / DynamicsDriver）
    md::mpi::DomainDecomposition* dd = nullptr;
    const md::PairList* pair_list = nullptr;
    bool is_last_step = false;
};
struct SamplingServices {   // 发布方：SamplingDriver::run
    mc::Transaction* proposal = nullptr;
    std::mt19937_64* rng = nullptr;
};
struct OptimizeServices {   // 发布方：OptimizeDriver::run / Optimizer 组件
    std::function<void(bool need_grad)> reeval_forces = nullptr;
    bool line_search_trial = false;
};
struct AnalysisServices {   // 发布方：AnalysisExecutor::install
    AnalysisExecutor* executor = nullptr;
};

struct Context {
    int n_threads = 0;  int mpi_rank = 0;  int mpi_size = 1;
    ATV::TimerRegistry* timer = nullptr;
    IterationStatus outer_status{};
    DynamicsServices md;  SamplingServices mc;
    OptimizeServices opt; AnalysisServices analysis;
};
```

头文件内置守卫：`static_assert(std::is_trivially_copyable_v<DynamicsServices/SamplingServices/AnalysisServices>)`、`static_assert(!std::is_polymorphic_v<Context>)`（`OptimizeServices` 因既有 `std::function` 豁免——它在拆分前就在 Context 里，非本 spec 新增间接）。

### CUDA 同构对齐

`cuda::Context`（include/atomiverse/cuda/component/context.cuh）：core 平铺保持 `stream` / `defer_stream` / `outer_status` / `host_box`；新增 `cuda::DynamicsServices md` 收纳 `pair_list`、`flat_pairs_requested`、exclusion CSR（`excl_offsets`/`excl_partners`/`s_excl`/`exclusions_baked`）、`n_force_rows`/`force_rows`。**裁决：exclusion CSR 与 force_rows 归 `md` slice**（而非另设 device slice）——它们由 `PairLJCut::build` 发布、`NeighborList::launch` 与多 GPU 力分派消费，和 `pair_list` 同属经典力管线的发布-消费链，同 slice 保持"按发布域分组"的一致性。**`host_box` 留 core**：它不是 nullable service 指针而是按值传入全部 kernel 族（integrator/bonded/pair/kspace/minimizer）的胞快照，语义属于核心执行状态。CUDA 无 SamplingDriver、无 `reeval_forces`、无 executor（NOTES 2026-06-06：CUDA 异步走 AsyncBridge，明确否决过 `cuda::Context.executor`），故无 mc/opt/analysis slice——同构指"同名 slice 同语义"，不是字段集相等。`torch::Context`（grad_enabled 单字段）无域字段，不动；但 `src/torch/neighbor/neighbor_list.cu` 借用 cuda::Context 发布 pair_list，属 cuda 读点一并迁移。

### 填充责任（语义不变，仅路径改写）

- `Backend::init/set_frame`（src/cpu/driver/backend.cpp:124-152,363）：`context_.md.dd`、`context_.analysis.executor`（timer 留 core）。
- `NeighborList::build/launch`：`ctx.runtime.md.pair_list`。
- `DynamicsDriver::run`（dynamics_driver.h:76,82）：`md.is_last_step`。
- `SamplingDriver::run`（sampling_driver.h:100-101）：`mc.proposal` / `mc.rng`。
- `OptimizeDriver::run`（optimize_driver.h:50）：`opt.reeval_forces`；LBFGS/CG line search（lbfgs.cpp:286,293、cg.cpp:210,217）：`opt.line_search_trial`。
- `AnalysisExecutor::install`（analysis_executor.cpp:71）：`ctx.analysis.executor = this`。
- CUDA：`cuda::PairLJCut::build` 发布 excl CSR、多 GPU comm 设 force_rows，全部走 `md.`。

各 category driver 只填自己的 slice；**未填 slice 的字段保持 nullptr = 服务不可用，组件用前判空**——既有语义不变。跨 slice 读是合法的（如 `comm::Migrate` 读 `opt.line_search_trial`；Minimize 的 force phase 读 `md.pair_list`）：slice 按**发布域**分组，不按独占消费者分组。

### 热路径不变量（链级门禁）

- `LaunchContext` 每 op 栈构造（backend.cpp:256）逐字不变；`BuildContext`/`LaunchContext`/`DestroyContext`（component.h:25-44）保持 `{Context&, Frame&, stage, step}` 形状，component.h 零改动。
- `ctx.runtime.md.pair_list` 编译为与今日相同的常量偏移载入——无 vtable、无堆、无 registry。
- Driver 命名全部保留（DynamicsDriver/SamplingDriver/OptimizeDriver 等）。

### 已否决备选

1. **拆成多个独立 Context 类型按 driver 传不同类型**——破坏 `ComponentBase` 单一 `launch(LaunchContext)` 签名，等于按 workflow fork 组件层次；一个 `PairLJCut` 将被迫为 MD/Minimize/MC 各编一份。否决。
2. **虚 service 接口 / service registry / service-locator**——热路径指针间接层 + 潜在分配，直接违反链级不变量（Context 保持 POD 形态、每 op 栈构造零成本）。否决。

## Files to create or modify

结构本体：
- include/atomiverse/cpu/component/context.h
- include/atomiverse/cuda/component/context.cuh

CPU 填充点：
- src/cpu/driver/backend.cpp
- include/atomiverse/cpu/driver/backend.h（doc 注释 `Context::dd` → `Context::md.dd`）
- include/atomiverse/cpu/driver/dynamics_driver.h
- include/atomiverse/cpu/driver/optimize_driver.h
- include/atomiverse/cpu/driver/sampling_driver.h
- src/cpu/async/analysis_executor.cpp
- include/atomiverse/cpu/async/analysis_executor.h

CPU 读点（组件头）：
- include/atomiverse/cpu/component/owned_range.h
- include/atomiverse/cpu/neighbor/neighbor_list.h
- include/atomiverse/cpu/potentials/resolve_view.h
- include/atomiverse/cpu/potentials/bonded_cache.h
- include/atomiverse/cpu/potentials/coul_long.h
- include/atomiverse/cpu/minimizer/global_reduce.h
- include/atomiverse/cpu/async/async_component_base.h
- include/atomiverse/cpu/mc/translate_move.h

CPU 读点（实现）：
- src/cpu/comm/halo_exchange.cpp
- src/cpu/comm/atom_migration.cpp
- src/cpu/comm/reduce.cpp
- src/cpu/comm/rcb.cpp
- src/cpu/comm/gather_for_io.cpp
- src/cpu/md/nose_hoover.cpp
- src/cpu/md/nose_hoover_npt.cpp
- src/cpu/md/init.cpp
- src/cpu/md/fix_csvr.cpp
- src/cpu/md/barostat.cpp
- src/cpu/minimizer/lbfgs.cpp
- src/cpu/minimizer/cg.cpp
- src/cpu/minimizer/opt_log.cpp
- src/cpu/mc/translate_move.cpp
- src/cpu/potentials/pair_lj_cut.cpp
- src/cpu/potentials/coul_long.cpp
- src/cpu/kspace/ewald.cpp
- src/cpu/kspace/pme.cpp
- src/cpu/async/async_component_base.cpp
- src/cpu/async/snapshot_publisher.cpp

CUDA / torch 读点与填充点：
- include/atomiverse/cuda/neighbor/neighbor_list.cuh
- include/atomiverse/cuda/observer/rdf.cuh
- src/cuda/neighbor/neighbor_list.cu
- src/cuda/potentials/pair_lj_cut.cu
- src/cuda/potentials/coul_long.cu
- src/cuda/comm/multi_gpu.cu
- src/torch/neighbor/neighbor_list.cu

tests / tools / benchmarks：
- tests/cpu/component/test_context_slices.cpp (new)
- tests/CMakeLists.txt
- tests/cpu/md/test_domain_decomposition.cpp
- tests/cpu/mc/test_mc_move.cpp
- tests/cpu/test_mc_scaffold.cpp
- tests/cpu/kspace/test_pme_coul.cpp
- tests/cpu/kspace/test_ewald_coul.cpp
- tests/mpi/test_topology.cpp
- tests/mpi/test_atom_migration.cpp
- tests/mpi/test_rcb.cpp
- tests/mpi/test_gather.cpp
- tests/mpi/test_halo_exchange.cpp
- benchmarks/md/pe/performance/polyethylene_tension.cpp
- tools/debug_amorphous/prof_amorphous.cpp
- tools/debug_amorphous/gather_frame.h

## Tasks

- [ ] Write failing tests for the sliced Context layout (tests/cpu/component/test_context_slices.cpp (new), registered as `test_cpu_component_context_slices` in tests/CMakeLists.txt): assert `ctx.md/.mc/.opt/.analysis` slices exist with the exact field mapping, default-constructed slices are all-null/false, `DynamicsServices`/`SamplingServices`/`AnalysisServices` are trivially copyable, and `Context` is non-polymorphic (RED before the header change)
- [ ] Implement the nested slices in `cpu::Context` (include/atomiverse/cpu/component/context.h): add `DynamicsServices{dd, pair_list, is_last_step}` / `SamplingServices{proposal, rng}` / `OptimizeServices{reeval_forces, line_search_trial}` / `AnalysisServices{executor}` + header static_asserts, delete the flat domain fields, keep n_threads/mpi_rank/mpi_size/timer/outer_status flat (clean break, no shims)
- [ ] Implement `cuda::DynamicsServices md` in `cuda::Context` (include/atomiverse/cuda/component/context.cuh): move pair_list, flat_pairs_requested, excl_offsets/excl_partners/s_excl/exclusions_baked, n_force_rows/force_rows into the slice; stream/defer_stream/outer_status/host_box stay flat
- [ ] Migrate the CPU driver/backend fill points to slice writes: src/cpu/driver/backend.cpp (md.dd, analysis.executor), include/atomiverse/cpu/driver/dynamics_driver.h (md.is_last_step), include/atomiverse/cpu/driver/optimize_driver.h (opt.reeval_forces), include/atomiverse/cpu/driver/sampling_driver.h (mc.proposal/mc.rng), src/cpu/async/analysis_executor.cpp + include/atomiverse/cpu/async/analysis_executor.h, include/atomiverse/cpu/driver/backend.h doc refs
- [ ] Migrate all CPU component read points to slice reads: comm (src/cpu/comm/{halo_exchange,atom_migration,reduce,rcb,gather_for_io}.cpp), md (src/cpu/md/{nose_hoover,nose_hoover_npt,init,fix_csvr,barostat}.cpp), minimizer (src/cpu/minimizer/{lbfgs,cg,opt_log}.cpp + include/atomiverse/cpu/minimizer/global_reduce.h), mc (src/cpu/mc/translate_move.cpp + include/atomiverse/cpu/mc/translate_move.h), potentials (src/cpu/potentials/{pair_lj_cut,coul_long}.cpp + include/atomiverse/cpu/potentials/{resolve_view,bonded_cache,coul_long}.h), kspace (src/cpu/kspace/{ewald,pme}.cpp), neighbor (include/atomiverse/cpu/neighbor/neighbor_list.h), async (src/cpu/async/{async_component_base,snapshot_publisher}.cpp + include/atomiverse/cpu/async/async_component_base.h), component (include/atomiverse/cpu/component/owned_range.h)
- [ ] Migrate the CUDA + torch read/fill points to `md.` slice access: src/cuda/neighbor/neighbor_list.cu, src/cuda/potentials/{pair_lj_cut,coul_long}.cu, src/cuda/comm/multi_gpu.cu, include/atomiverse/cuda/neighbor/neighbor_list.cuh, include/atomiverse/cuda/observer/rdf.cuh, src/torch/neighbor/neighbor_list.cu
- [ ] Migrate test/tool/benchmark read points: tests/cpu/md/test_domain_decomposition.cpp, tests/cpu/mc/test_mc_move.cpp, tests/cpu/test_mc_scaffold.cpp, tests/cpu/kspace/{test_pme_coul,test_ewald_coul}.cpp, tests/mpi/{test_topology,test_atom_migration,test_rcb,test_gather,test_halo_exchange}.cpp, benchmarks/md/pe/performance/polyethylene_tension.cpp, tools/debug_amorphous/{prof_amorphous.cpp,gather_frame.h}
- [ ] Add docstring per doxygen for every new slice struct and moved field in include/atomiverse/cpu/component/context.h and include/atomiverse/cuda/component/context.cuh (@brief with publisher, null-semantics, ownership; carry over the Frame-purity rationale verbatim where it applies)
- [ ] Verify against the regression anchors: MPI suite (ctest -R test_mpi, DD/halo via ctx.md), MC Metropolis (ctest targets from tests/cpu/mc/ + test_mc_scaffold, ctx.mc.rng/proposal), minimizer line-search (ctest -R test_cpu_minimizer, ctx.opt.reeval_forces/line_search_trial), and the PE stage-5 performance gate (benchmarks/md/pe/performance/polyethylene_tension, per-step wall time vs parent-commit baseline)
- [ ] Run full check + test suite (cmake --build build -j + ctest --test-dir build --output-on-failure; CUDA build on the GH200 node when available)

## Testing strategy

- **Happy path：** 全量 ctest（~16 s 门禁）在 CPU 构建上 100% 通过，含新 `test_cpu_component_context_slices`（slice 布局、默认空值、trivially-copyable/non-polymorphic 守卫）。
- **边界情形：** 默认构造 `Context` 的全部 slice 指针为 nullptr（服务不可用路径——非 MD driver 下组件判空跳过，既有各域测试覆盖）；`md.is_last_step` 仅在最后一步为 true（reduce cadence 测试）；`opt.line_search_trial` 在试探区间外恒为 false（Migrate 不被误抑制）；MC 字段在非 Sampling driver 下恒 null。
- **域校验（science.required，回归锚——数值必须逐位不变）：**
  - MPI 回归：`ctest -R test_mpi`（halo/migration/topology/gather/rcb + pw k-point），DD/halo 读 `ctx.md.*`。
  - MC Metropolis 回归：tests/cpu/mc/ 全部目标 + test_mc_scaffold（`ctx.mc.rng`/`ctx.mc.proposal` 单流可复现）。
  - Minimizer line-search 回归：`ctest -R test_cpu_minimizer`（LBFGS/CG 经 `ctx.opt.reeval_forces` 的试探求值与迁移抑制）。
  - CUDA：GH200 节点 `ATV_BUILD_CUDA=ON` 构建 + cuda ctest 目标通过（neighbor/potentials/multi-GPU 读 `ctx.md.*`）。
- **性能门禁：** PE stage-5 基准（polyethylene_tension）每步耗时相对父 commit 基线无回归（同机、Release、-march=native）——LaunchContext 栈构造与字段常量偏移访问不变的实证。

## Out of scope

- 文档表更新（CLAUDE.md 架构节 / NOTES 提炼 / developer docs 的 Context 描述）→ 统一归 arch-refactor-09-harness。
- `torch::Context` 重设计（现仅 grad_enabled，无域字段可分组）。
- 核心字段的进一步分组（`Core core` 嵌套形态）——已考虑并弃用，理由见 Design；如后续 driver 类别增至需要时再议。
- namespace 与目录对齐、`cpu::pw` 脱离 `scf::Phase`（既有登记的独立 follow-up）。
- Frame 字段/residency 变更 → arch-refactor-06-frame。
- 已否决备选（per-driver Context 类型、虚 service 接口/registry）不实现。
- 任何行为、数值、性能语义改动：零（纯重组）。
