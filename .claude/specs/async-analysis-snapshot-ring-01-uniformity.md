---
title: 异步分析快照环统一化（链 1/3：单一路径 + 冻结-派发不变量）
status: code-complete
created: 2026-06-06
---

# 异步分析快照环统一化（链 1/3：单一路径 + 冻结-派发不变量）

## Summary

本链将 CPU 端异步分析框架（`AnalysisExecutor` 快照环 + `Publish` 发布器 + `AsyncComponentBase` 派发器）收敛为唯一一条统一路径：所有组件走同一冻结-派发流程，不再存在任何按组件粒度的特化。背压策略被强制为单一的 `drop_oldest`（移除 `block` 与 `grow`，以及组件级 `backpressure()` 覆盖），环深度 K 成为唯一一个全局配置字段；捕获掩码（capture mask）在 `build()` 时按所有订阅者求并集一次性确定，删除每步重算路径。同时以一个先红（RED-first）的回归测试锁定冻结-派发正确性不变量：任何被投递的快照槽，其全部原子字段都来自同一个仿真步 S(N)，且 `slot.step == N`。用户可见结果是更简单、更可预测的异步分析行为，并为链 02（CUDA 设备环）确定其依赖的执行器形状。

## Domain basis

`science.required = true`，但本链**不引入任何新的物理方程**。其领域基础是**快照一致性不变量**（snapshot-consistency invariant），可被证伪：

- 设解析轨迹 `x_i(step) = f(tag_i, step)`（例如 `f(tag, s) = tag + s`），即每个原子的坐标由其全局 tag 与步号唯一确定。
- 不变量：对任意被投递到 worker 的快照槽 `slot`，其所有原子 `i` 满足 `slot.x[i] == f(slot.tag[i], slot.step)`。换言之，一个槽内不允许出现「撕裂」（部分原子来自步 N、部分来自步 N+1）或步号错配（`slot.step` 与原子内容不符）。

CPU 已天然满足该不变量的原因（须在 Design 中重申）：`Publish` 在热路径线程上于 `md::Phase::end_of_step` 同步完成对帧的内存拷贝（冻结），该冻结按程序顺序先于下一步对 `frame` 的任何修改，因此投递给 worker 的快照是某一步的一致快照——这正是「冻结-然后-派发」（freeze-then-dispatch）。

单位：坐标 [Å]，速度 [Å/fs]，力 [eV/Å]（沿用 `frame_snapshot.h` 既有约定，本链不改）。

## Design

涉及实体与改动：

1. **`core/types.h` — `BackpressurePolicy` 枚举处置（推荐：整体删除）。**
   核查结果：该枚举的全部引用仅来自 CPU 分析子系统——`analysis_executor.{h,cpp}`、`async_component_base.h`、`backend.cpp:132`、`tests/cpu/test_analysis_executor.cpp`；Grep 命中的 CUDA 文件中的 `::block`/`::grow` 为无关 token，与本枚举无关。既然背压收敛为唯一策略，单值枚举属于死重量。**推荐方案：整体删除 `enum class BackpressurePolicy`，并移除 `AnalysisExecutor::Config::policy` 字段**，drop_oldest 成为隐含且唯一的行为。本 spec 后续 Files/Tasks 按此推荐执行。

2. **`AnalysisExecutor` — 折叠为纯 drop_oldest。**
   `acquire_slot_locked()` 退化为单分支：槽空闲则占用并环形推进；槽忙则计入 `n_dropped_` 后覆盖最旧槽。删除 `grow_ring_locked`、`wait_for_free_slot_locked`、`cv_slot_free_`、`ring_grew_`、`ring_size_grew()`、`policy_` 成员及 `grow()` 分支引用的 `ring_capacity_` 动态扩张语义（容量在 `init()` 后恒定 = K）。`Config::ring_size` 重命名语义为「全局环深度 K」并补 Doxygen：K 越大丢帧越少、内存越高。`shutdown()` 中对 `cv_slot_free_.notify_all()` 的调用一并移除。

3. **`AsyncComponentBase` — 删除组件级背压覆盖。**
   删除 `virtual BackpressurePolicy backpressure()`（约第 69 行）。`launch()` 仍为 `final` 的 ~30ns 纯派发器（仅 submit，无物理）——保持不变。`capture_mask()` 纯虚保留。

4. **`Publish` — 仅用 build 期静态掩码。**
   `launch()` 使用在 `build()` 期由 `mask_union_all` 计算并缓存的静态并集掩码；删除 `mask_union_firing`（约第 45 行）与每步掩码重算。注意：`Publish::build()` 当前因发布器先于订阅者 build 而拿不到订阅者列表——本链将掩码并集的求值时机移到首次 `launch()` 的惰性一次性计算（或在订阅者全部 build 完成后的统一时点），缓存到 `Publish` 成员，后续 `launch()` 复用，绝不每步重算。`subscriber_fires` 仍用于决定本步是否有订阅者触发（是否走快照），但触发与否不再改变所拷贝的字段集合。

5. **重申冻结-派发不变量在 CPU 成立**（见 Domain basis），并以 RED 测试固化为可执行回归。

生命周期/所有权：`AnalysisExecutor` 仍由 `cpu::Backend` 独占持有，`init()`/`shutdown()` 括起，环容量恒定，无 `launch()` 内分配。

## Files to create or modify

- `include/atomiverse/core/types.h` — 删除 `BackpressurePolicy` 枚举。
- `include/atomiverse/cpu/analysis_executor.h` — 删除 `policy_`、`ring_grew_`、`cv_slot_free_`、`ring_size_grew()`、`grow_ring_locked()`、`wait_for_free_slot_locked()` 声明；`Config` 删 `policy` 字段，`ring_size` 补 K 语义 Doxygen；`acquire_slot()`/`acquire_slot_locked()` Doxygen 改为单策略。
- `src/cpu/analysis_executor.cpp` — `acquire_slot_locked()` 折叠为纯 drop_oldest；删除 `grow_ring_locked`、`wait_for_free_slot_locked` 定义及 `init()`/`shutdown()` 中相关引用。
- `include/atomiverse/cpu/async_component_base.h` — 删除 `backpressure()` 虚函数。
- `src/cpu/snapshot_publisher.cpp` — 删除 `mask_union_firing`；`Publish` 改用一次性缓存的 `mask_union_all`。
- `include/atomiverse/cpu/snapshot_publisher.h` — 为 `Publish` 增加缓存掩码成员 / 一次性求值标志（如需）。
- `src/cpu/backend.cpp` — 删除第 132 行 `cfg.policy = ...` 赋值。
- `tests/cpu/test_analysis_executor.cpp` — 删除 `BlockBackpressure`、`GrowBackpressure` 用例及所有 `cfg.policy` 赋值；保留并强化 drop_oldest / round-robin / K 全局深度断言。
- `tests/cpu/test_snapshot_publisher.cpp` — 新增 RED 回归：每步捕获 + 人为放慢 worker，在解析轨迹上断言每个被投递槽的原子全部属于 `slot.step`；新增静态掩码断言（掩码于 build 期定一次、无每步路径）。

## Tasks

- [x] Write failing snapshot-consistency regression for freeze-then-dispatch (tests/cpu/test_snapshot_publisher.cpp): every-step capture + slowed on_snapshot over analytic trajectory x_i(step)=f(tag_i,step), assert each delivered slot's atoms all belong to slot.step
- [x] Write failing static-mask + drop_oldest-only tests (tests/cpu/test_snapshot_publisher.cpp, tests/cpu/test_analysis_executor.cpp): assert mask computed once at build, remove Block/Grow cases so referencing BackpressurePolicy::block|grow fails to compile, assert ring depth K is global
- [x] Remove BackpressurePolicy enum from include/atomiverse/core/types.h and drop Config::policy field in include/atomiverse/cpu/analysis_executor.h, delete cfg.policy assignment in src/cpu/backend.cpp
- [x] Collapse acquire_slot_locked() to pure drop_oldest in src/cpu/analysis_executor.cpp and delete grow_ring_locked/wait_for_free_slot_locked/cv_slot_free_/ring_grew_/ring_size_grew()/policy_ from analysis_executor.{h,cpp}
- [x] Remove backpressure() virtual from include/atomiverse/cpu/async_component_base.h
- [x] Replace per-step mask with build-time cached mask_union_all in src/cpu/snapshot_publisher.cpp and delete mask_union_firing, adding cache member to include/atomiverse/cpu/snapshot_publisher.h
- [x] Add Doxygen for Config::ring_size (global ring depth K: larger K = fewer drops, more memory) and updated acquire_slot semantics with units
- [x] Run full check + test suite (CPU, ATV_BUILD_CUDA=OFF, login node)

## Testing strategy

- **Happy path:** 既有 `AcquireSlotRoundRobin`、`DropOldestBackpressure`、`InitShutdownClean`、`ShutdownWithInFlightJobs` 在折叠后仍通过（后者去掉 `cfg.policy=block`，改用默认 drop_oldest）；`Publish` 仍正确将 positions/box/step 写入快照。
- **边界 / 回归（RED-first）：**
  - 冻结-派发一致性：每步捕获、worker 在 `on_snapshot` 中 sleep 人为落后，解析轨迹 `x_i(step)=f(tag_i,step)`；断言每个被投递槽满足 `slot.x[i] == f(slot.tag[i], slot.step)`，无撕裂、无步号错配。
  - drop_oldest-only：移除 `block`/`grow` 用例后，源码中任何 `BackpressurePolicy::block|grow` 引用编译失败（枚举已删）；环深度 K 为全局单值且容量恒定（`ring_capacity()` 在加载后不随饱和变化）。
  - 静态掩码：掩码仅在 build 期求值一次；不存在每步重算路径（断言并集字段集合在不同步保持稳定，且无 `mask_union_firing` 符号）。
- **领域验证：** 上述一致性断言即为 Domain basis 不变量的可证伪检验，构成本链的科学验证（无新方程，故仅此一项）。

## Out of scope

- CUDA 设备端快照环（属链 02 `async-analysis-snapshot-ring-02-*`，依赖本链确定的执行器形状）。
- 任何按组件粒度的背压/掩码特化——已被设计裁决永久排除，不在任何链内恢复。
- completeness（may-drop-frames）的按系统配置轴——本链仅锁定 accuracy 不变量与单一路径；completeness 轴的暴露留待后续链（已确认非本链目标）。
- 性能/吞吐量调优（K 的自动选择、worker 数自适应）——本链不改既有默认（K=3、1 worker）。
