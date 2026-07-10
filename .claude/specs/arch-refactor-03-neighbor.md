---
title: arch-refactor-03-neighbor — neighbor 契约从 ComponentBase 拆出为 opt-in 协议
status: code-complete
date: 2026-07-10
chain: arch-refactor (03/09)
depends_on: [arch-refactor-02-halo]
---

# arch-refactor-03-neighbor — neighbor 契约从 ComponentBase 拆出为 opt-in 协议

## Summary

当前 neighbor 契约以三个虚函数摊在 `ComponentBase` 上（`needs_neighbor_list()` component.h:97、`provides_neighbor_list()` :108、`neighbor_cutoff_hint()` :130），CUDA 侧缺 `neighbor_cutoff_hint`，且 `CoulLong` 只 override 了 `needs_neighbor_list` 而未转发 cutoff——halo 尺寸探测对它静默失效。本篇沿用 arch-refactor-02-halo 建立的 opt-in 协议模式（协议接口 + `ComponentBase` 上单个可空 `as_<protocol>()` 类型化访问器，无 `dynamic_cast`），把 neighbor 三件套拆为消费侧 `NeighborConsumer` 与供给侧 `NeighborProvider` 两个独立协议，CPU/CUDA 同形状对齐（CUDA 顺带补齐 cutoff 上报），`scan_pending_infra_()` / `cuda::md::MD::install_infra_()` 的探测随迁为协议查询。clean-break：旧虚函数直接删除，无兼容垫层。用户可见结果不变——`set_pair` 后 NeighborList 缺失照旧在 `build()` 抛错点名、MPI halo 照旧自动注入且尺寸公式不变（1.1× pair cutoff）；唯一有意的行为增量是 `CoulLong` 开始上报真实 cutoff。协议查询只发生在 `build()` 期扫描，step 热路径零新增虚调用，`cell_list.cpp` 重建算法一行不动。

## Domain basis

本篇不引入新物理。需保持不变的领域契约：Verlet 邻表语义对齐 LAMMPS `neighbor`/`neigh_modify`（neighbor_list.h 头注释审计表，对照 lammps/src/neighbor.cpp `Neighbor::decide()`）；MPI halo cutoff = 1.1 × pair cutoff 的 skin 带（md.h:279-312，防 ghost 在 rebuild 周期间漂出 `tag_to_local`，PE NPT 下实测 +1–3 eV/rebuild 的 MTK 扩展哈密顿量跳变依据）；Ewald 实空间分裂要求邻表 cutoff ≥ `CoulLong` 实空间 cutoff（coul_long.h 头注释）。单位一律 Å（`core/units.h` 约定）。

## Design

**协议实体（新增，纯接口、无状态、不改任何所有权）：**

- `ATV::cpu::NeighborConsumer` — `virtual Real neighbor_cutoff() const = 0`：该组件要求邻表覆盖的相互作用 cutoff [Å]。返回 `0` 保留旧 `neighbor_cutoff_hint` 语义 = "不声明 cutoff，halo 注入跳过尺寸推导"。
- `ATV::cpu::NeighborProvider` — 供给侧标记协议（本篇为纯 tag；发布途径仍是 `NeighborList::build()/launch()` 写 `Context::pair_list`，非虚函数，归属不动——见 Out of scope）。
- `ComponentBase` 上各加一个可空类型化访问器（02 的探测机制同款）：
  `virtual const NeighborConsumer* as_neighbor_consumer() const { return nullptr; }`、
  `virtual const NeighborProvider* as_neighbor_provider() const { return nullptr; }`；
  协议实现者 override 返回 `this`。净效果：ComponentBase 上 3 个 neighbor 虚函数 → 2 个访问器，行为虚函数全部移入协议。
- CUDA 侧同形状：`cuda::NeighborConsumer::neighbor_cutoff()`（`double`，与 `cuda::Pair::pair_cutoff()` 一致）+ `cuda::NeighborProvider` + 两个访问器。`cuda::Pair` 借此首次上报 cutoff，补齐 CPU/CUDA 探测面差异。

**迁移映射：**

- `cpu::PairEvaluator`（pair_evaluator.h:36-49，既有 "NeighborConsumer 雏形"）→ 继承 `NeighborConsumer`，`neighbor_cutoff()` 转发既有纯虚 `pair_cutoff()`，`as_neighbor_consumer()` 返回 `this`；删除 `needs_neighbor_list()`/`neighbor_cutoff_hint()` override。
- `cpu::PairMorse`（O(N²) 演示，现以 `needs_neighbor_list() final { return false; }` 退出门禁）→ 改为 `as_neighbor_consumer() const final { return nullptr; }`，opt-out 语义保留且更直白。
- `cpu::CoulLong` → 实现 `NeighborConsumer`，`neighbor_cutoff()` 返回 `cutoff_`。修复现状 hint==0 的静默缺口（deliberate delta，见 Testing）。
- `cpu::md::NeighborList` / `cuda::NeighborList` → 实现 `NeighborProvider`，删除 `provides_neighbor_list()` override。
- `CPUDriver::InfraScan`（cpu_driver.h:333-377）字段随语义改名：`pair_evaluator_name/pair_evaluator_cutoff` → `neighbor_consumer_name/neighbor_consumer_cutoff`；扫描改为 `if (const auto* nc = c.as_neighbor_consumer()) …` / `if (c.as_neighbor_provider()) …`。**保持"取第一个 consumer 的 cutoff"语义不变**（不改成 max——那是行为变更，不属本篇）。`throw_if_missing_neighbor_list_` 报错文案随迁协议措辞，保留"点名 offender + 给出 set_neighbor 指引"契约。四个 CPU 调用点（md.h / minimize.h / evaluate.h / mc.h）随字段名更新，halo 注入逻辑（1.1× 系数、已注册 halo 去重）逐行不变。`cuda::md::MD::install_infra_()`（md.cuh:274-327）同构迁移。

**备选与决策理由（定案：消费/供给分开的双协议）：**

1. *单一 `NeighborParticipant` 双角色协议* — 否。消费者群体（pair styles、CoulLong）与供给者群体（NeighborList）不相交，双角色协议迫使每个实现者背上另一半的无意义默认实现，恰是本链要清除的 fat-base 病；且与 02 的 halo send/recv 分侧先例不对称。
2. *保留 ComponentBase 布尔虚函数、只补 CUDA hint* — 否。违背链目标（ComponentBase 瘦身、契约按协议聚合），且 cutoff 与 needs 两个虚函数天然属于同一消费契约，分离存放已经造成过 CoulLong 的转发遗漏。
3. *`dynamic_cast` 探测协议接口* — 否。仓库既有明文契约"capability-typed、无 dynamic_cast"（cpu_driver.h:350-355、2026-07-02 note），访问器方案零 RTTI 且与 02 探测机制同构。

**热路径论证：** 协议访问器仅在 `build()` 期的一次 pending 扫描中被调用；`NeighborList::launch()`（rebuild 判据 + CSR 复用）与 `Backend::step()` 虚调度结构均无改动；`cell_list.cpp` 不在触碰范围（其优化是独立工作流 `.claude/specs/cell-list-omp-finalize.md`）。Neigh 占 24-44% step 的既有评估不受影响，以 PE stage-5 基准门禁背书。本篇假定 02 已落地（halo 探测已走 `HaloParticipant` 协议）；若实际落地顺序反转，本篇不得代做 02 的 halo 迁移。

## Files to create or modify

- include/atomiverse/cpu/component/neighbor_protocol.h (new) — `NeighborConsumer` + `NeighborProvider`（component 支柱层，与 pair_evaluator.h 同层，镜像 02 的 halo 协议存放）
- include/atomiverse/cpu/component/component.h — 删三虚函数，加两个 `as_*` 访问器
- include/atomiverse/cpu/component/pair_evaluator.h — PairEvaluator 实现 NeighborConsumer
- include/atomiverse/cpu/potentials/pair_morse.h — opt-out 改为访问器返回 nullptr
- include/atomiverse/cpu/potentials/coul_long.h — 实现 NeighborConsumer（补 cutoff 上报）
- include/atomiverse/cpu/neighbor/neighbor_list.h — 实现 NeighborProvider
- include/atomiverse/cpu/driver/cpu_driver.h — InfraScan 字段改名 + 扫描/报错随迁协议查询
- include/atomiverse/cpu/md/md.h — install_defaults 调用点 + 头注释随迁
- include/atomiverse/cpu/minimizer/minimize.h — 调用点随迁
- include/atomiverse/cpu/driver/evaluate.h — 调用点随迁
- include/atomiverse/cpu/mc/mc.h — 调用点随迁
- include/atomiverse/cuda/component/neighbor_protocol.cuh (new) — CUDA 镜像协议（CUDA-free 头）
- include/atomiverse/cuda/component/component.cuh — 删两虚函数，加两个访问器
- include/atomiverse/cuda/potentials/pair.cuh — Pair 实现 NeighborConsumer（cutoff 上报补齐）
- include/atomiverse/cuda/neighbor/neighbor_list.cuh — 实现 NeighborProvider
- include/atomiverse/cuda/md/md.cuh — install_infra_ 随迁 + 头注释更新
- tests/cpu/test_install_infra_hook.cpp — 默认值断言迁到访问器
- tests/cpu/md/test_md_set_api.cpp — LJ/NPT/邻表断言迁到协议
- tests/cuda/md/test_install_infra_capability.cu — CUDA 侧同构迁移 + cutoff 断言
- .claude/architecture.md — 刷新引用已删虚函数的段落

## Tasks

- [x] Write failing tests for the neighbor opt-in protocol: update tests/cpu/test_install_infra_hook.cpp（`as_neighbor_consumer()`/`as_neighbor_provider()` 默认 nullptr）、tests/cpu/md/test_md_set_api.cpp（PairLJCut → consumer 且 `neighbor_cutoff()==` 构造 cutoff；NoseHooverNPT → nullptr；NeighborList → provider；PairMorse → nullptr opt-out；CoulLong → consumer 上报 cutoff）、tests/cuda/md/test_install_infra_capability.cu（CUDA 同形状 + `cuda::Pair` cutoff 上报）
- [x] Implement CPU protocols in include/atomiverse/cpu/component/neighbor_protocol.h (new) and thin ComponentBase in include/atomiverse/cpu/component/component.h（删 needs_neighbor_list/provides_neighbor_list/neighbor_cutoff_hint，加双访问器；doxygen `@brief`/`@return` 含 Å 单位）
- [x] Migrate CPU consumers and supplier: include/atomiverse/cpu/component/pair_evaluator.h、include/atomiverse/cpu/potentials/pair_morse.h、include/atomiverse/cpu/potentials/coul_long.h、include/atomiverse/cpu/neighbor/neighbor_list.h
- [x] Migrate the CPU driver scan: include/atomiverse/cpu/driver/cpu_driver.h（InfraScan 字段改名 + 协议查询 + throw 文案）及调用点 include/atomiverse/cpu/md/md.h、include/atomiverse/cpu/minimizer/minimize.h、include/atomiverse/cpu/driver/evaluate.h、include/atomiverse/cpu/mc/mc.h（保持 first-consumer cutoff 语义与 1.1× halo 公式逐行不变）
- [x] Implement the CUDA mirror: include/atomiverse/cuda/component/neighbor_protocol.cuh (new)、include/atomiverse/cuda/component/component.cuh、include/atomiverse/cuda/potentials/pair.cuh、include/atomiverse/cuda/neighbor/neighbor_list.cuh、include/atomiverse/cuda/md/md.cuh
- [x] Verify auto-injection regressions: set_pair 无 set_neighbor 时 MD/Minimize/Evaluate/MC `build()` 照旧抛错点名 offender；set_neighbor 后构建成功；MPI halo 注入（尺寸/去重）与迁移前一致，CoulLong-only 的 halo 尺寸增量以断言固化
- [x] Verify hot-path non-regression [cell_list.cpp/plan.h/backend.cpp git-diff 零改动已验；PE stage-5 计时对照移交 ac-007 owed /mol:bench]: `git diff` 确认 src/cpu/neighbor/cell_list.cpp 与 NeighborList::launch() rebuild 判据零改动；运行 benchmarks/md/pe/performance/polyethylene_tension（CPU，--scale 缩短）对比迁移前基线的 step 计时与 Neigh 份额
- [x] Update docs: .claude/architecture.md 及各头文件注释中对已删虚函数的引用全部随迁（rg 三个旧符号零命中）
- [x] Run full check + test suite [CPU 109/111（2 个已档案化的 macOS 既有失败）；CUDA/GPU 节点移交 ac-008]（登录节点 CPU 全量 ctest；GH200 批作业 CUDA build + ctest）

## Testing strategy

- **Happy path：** 全量 `ctest --test-dir build --output-on-failure`（~16 s 门禁，Release + `-DATV_BUILD_RHF=ON`，MPI module 环境）；CUDA 侧在 GH200 批作业内 `-DATV_BUILD_CUDA=ON` 构建 + ctest。
- **边界：** PairMorse（consumer opt-out）无邻表构建成功；四类 driver 缺邻表均抛错且报错点名 offender、附 `set_neighbor` 指引；显式注册 halo 时不重复注入（02 协议路径，行为沿用）；`set_neighbor` 幂等/异名拒绝不受影响。
- **链级不变量（热路径门禁）：** step 虚调度结构不变——协议查询仅在 build 期扫描，`launch()` 零新增虚调用；`cell_list.cpp` 重建算法零改动（git diff 为空）；PE stage-5 基准回归门禁（polyethylene_tension，见验收 ac-007）。
- **领域校验：** 无新物理；以 MPI 能量 bruteforce 校验 halo 尺寸链路未破——tests/mpi/test_lj_energy_bruteforce.cpp、tests/mpi/test_coul_energy_bruteforce.cpp、tests/mpi/test_halo_exchange.cpp 全绿。**有意行为增量：** CoulLong 现在上报真实 cutoff，CoulLong-only 多 rank 场景从"halo 尺寸推导被跳过"变为按 1.1×cutoff 注入——这是修复而非漂移，在 test_coul_energy_bruteforce 语境下以断言固化。
- **命名保留：** Driver 命名不动（`scan_pending_infra_` / `install_infra_` / `install_defaults` 钩子名保持）。

## Out of scope

- `cell_list.cpp` 重建算法优化（OMP finalize、候选环 hoisting）→ 独立工作流 `.claude/specs/cell-list-omp-finalize.md`，本篇只动接口声明层。
- `Context::pair_list` 的归属与发布途径（非虚函数发布、Frame-purity 取舍）→ arch-refactor-05-context。
- 引擎层（`cuda::md::MultiGpu`）的邻表/行分割处理 → arch-refactor-04-engines。
- "first-consumer cutoff" → "max over consumers" 的语义升级（行为变更，需独立评估）。
- torch 后端（自有 LJ 层级，不含本组虚函数，零触碰）。
