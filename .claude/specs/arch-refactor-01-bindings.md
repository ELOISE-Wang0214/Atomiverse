---
title: arch-refactor-01-bindings — DefaultBinding 增加 PhaseKind 判别，注册期拦截跨 workflow 误绑
status: code-complete
date: 2026-07-10
chain: arch-refactor (01/09)
depends_on: []
---

# arch-refactor-01-bindings：DefaultBinding 增加 PhaseKind 判别，注册期拦截跨 workflow 误绑

## Summary

`DefaultBinding` 目前只有两个裸 `int`（phase、order），把五个数值全部从 0 起步、彼此重叠的
phase 枚举擦除成不可判别的整数。后果是：一个声明了 `evaluate::Phase::evaluate(=1)` 默认绑定的
Evaluator 组件被 `add()` 进 MD driver 时，静默落在 `md::Phase::post_integrate(=1)`；进
Minimize 时落在 `minimizer::Phase::line_search(=1)` —— 组件被注册、被执行，但执行在错误的
生命周期槽位上，既不报错也不产出正确物理。本 spec 做两件事：其一，把两个违约组件
（`cpu::semiempirical::AM1Evaluator`、`cuda::scf::SCFEvaluator`）的硬编码 `default_bindings()`
删除，对齐"Evaluator 一律显式 `.on()` 绑定"的既有约定；其二，给 `DefaultBinding` 增加
`PhaseKind` 判别字段，`Driver<PhaseEnum>::add()` 在注册期校验判别与本 driver 的 workflow
匹配，不匹配立即 throw 且错误信息点名双方 workflow。全树约 50 处 `default_bindings()`
override 机械加 tag，一次 commit 迁完所有 call site，无兼容层、无 `[[deprecated]]`。
本篇是 arch-refactor 链（9 篇）的链首，不依赖任何后续 sub-spec。

## Domain basis

本 spec 纯架构改造、不涉及物理方程；依据为下列已核实的代码事实（file:line）：

- `include/atomiverse/core/binding.h:8-11` — `DefaultBinding{int phase; int order}`，类型擦除源头。
- `include/atomiverse/facade/driver.h:137-160` — `Driver::add()` 两个所有权分支均只查
  `bindings.empty()`（:141-147 / :150-154），不校验 phase 归属。
- `include/atomiverse/facade/plan.h:78-109` — `PlanBuilder::build()` 仅做 `[0, phase_count)`
  范围检查（:91-94）；重叠数值全部落在合法范围内，永远查不出来。
- `include/atomiverse/facade/stage.h` — 5 个 phase 枚举全部 0 起步数值重叠：`md::Phase`
  10 相（:22-34）、`minimizer::Phase` 5 相（:59-66）、`evaluate::Phase` 3 相（:89-94）、
  `vib::Phase` 7 相（:129-138）、`sampling::Phase` 3 相（:167-174）。
- 违约组件①：`include/atomiverse/cpu/semiempirical/am1_evaluator.h:31-33` 硬编码
  `evaluate::Phase::evaluate(=1)`。
- 违约组件②：`include/atomiverse/cuda/scf/scf_evaluator.cuh:105-107` 同样硬编码。
- 既有约定：CPU `SCFEvaluator`/`PWEvaluator`/`GPWEvaluator` 不声明 `default_bindings`；
  `include/atomiverse/cpu/async/async_component_base.h:10-16` 的注释明言该隐患并规定
  "本基类不声明默认绑定，跨 driver 必须显式 `.on()`"。
- 门禁安全性依据：跨 driver 复用的 comm 组件（Migrate/HaloSend/HaloRecv/Reduce）进入
  Minimize 走的是显式 `.on(minimizer::Phase::force).add(...)` 路径
  （`include/atomiverse/cpu/minimizer/minimize.h:120-139`）；MD 自身的自动注入
  （`include/atomiverse/cpu/md/md.h:96,273-311`）添加的均为 MD-kind 组件；categorical
  setter 走 `this->on(phase).add(...)`（`include/atomiverse/cpu/driver/cpu_driver.h:324`）。
  故判别门禁不会误伤任何现存自动注入路径。

## Design

**新增判别类型（放 `core/binding.h`，理由见下）：**

```cpp
/// core/binding.h
enum class PhaseKind { MD, Minimizer, Evaluate, Vib, Sampling };

constexpr const char* phase_kind_name(PhaseKind k);  // "MD"/"Minimizer"/"Evaluate"/"Vib"/"Sampling"

struct DefaultBinding {
    PhaseKind kind; ///< 该 phase 所属 workflow（判别字段，置于首位、无默认值）
    int phase;      ///< Target phase index.
    int order;      ///< Ordering within phase (lower = earlier).
};
```

位置选 `core/binding.h` 而非 `facade/stage.h` 的理由：`DefaultBinding` 定义在
`core/binding.h`，且被 CPU/CUDA 两个 `ComponentBase` 与 `facade/driver.h` 共同包含；层次方向
是 core → facade → backend，若 `PhaseKind` 定义在 `facade/stage.h`，`core/binding.h` 就必须
反向包含 facade 头，层次倒置。判别枚举与承载它的结构体同头即可保持 `core/binding.h`
零依赖、CUDA-free。

`kind` 置于首位且**无默认成员初始化**是 clean-break 的机械保障：所有未迁移的旧式两元
`{phase, order}` 聚合初始化会因 `int` 无法隐式转换为 scoped enum 而编译失败，编译器替我们
枚举出全部 ~50 个待迁移站点，不存在静默漏迁。

**每个 PhaseEnum 关联其 PhaseKind（放 `facade/stage.h`）：** 枚举定义在哪里，映射就放在哪里；
且 `stage.h` 已有 `phase_names(Phase)` 的 ADL tag-dispatch 惯例，照抄该形状：

```cpp
/// facade/stage.h（每个 phase namespace 内各一）
namespace md        { constexpr PhaseKind kind_of(Phase) { return PhaseKind::MD; } }
namespace minimizer { constexpr PhaseKind kind_of(Phase) { return PhaseKind::Minimizer; } }
namespace evaluate  { constexpr PhaseKind kind_of(Phase) { return PhaseKind::Evaluate; } }
namespace vib       { constexpr PhaseKind kind_of(Phase) { return PhaseKind::Vib; } }
namespace sampling  { constexpr PhaseKind kind_of(Phase) { return PhaseKind::Sampling; } }

/// facade/stage.h（ATV 命名空间层）——组件声明默认绑定的唯一入口
template <typename PhaseEnum>
constexpr DefaultBinding make_binding(PhaseEnum p, int order = 0) {
    return {kind_of(p), static_cast<int>(p), order};  // kind 由枚举类型推出，杜绝手写抄错
}
```

`stage.h` 新增 `#include "atomiverse/core/binding.h"`（facade→core，方向正确）。

**注册期校验（`facade/driver.h` `Driver::add()`）：** 两个所有权分支（lvalue/rvalue）共用一个
私有校验助手；对每条默认绑定检查 `db.kind == kind_of(static_cast<PhaseEnum>(0))`（ADL 按
`PhaseEnum` 解析），不匹配即 throw，错误信息点名绑定所属 workflow 与本 driver 类别，例如：

```
Component 'am1_evaluator' declares default bindings for the Evaluate workflow,
but this driver executes MD phases. Bind explicitly with .on(phase).add().
```

`.on(phase).add()` 路径（`PhaseBuilder`）的 phase 参数是 `PhaseEnum` 强类型，编译期即安全，
不需要也不做运行期判别。`detail::PendingBinding`、`PlanBuilder::Entry`、`Plan`/`StageSlot`/
`Operation` 全部保持 `int phase` 不变——判别在 `add()` 处消费完毕即擦除，Plan 编译与执行
路径零改动。

**违约组件对齐约定：** `AM1Evaluator` 与 `cuda::scf::SCFEvaluator` 删除 `default_bindings()`
override，回落到基类空实现，与 CPU `SCFEvaluator`/`PWEvaluator`/`GPWEvaluator` 一致；此后
`add()` 走既有的 "no default bindings" throw，裸 `add()` 调用点（3 个测试/示例文件）迁到
`.on(evaluate::Phase::evaluate).add(...)` 或 `single_point`。删除后这两个组件不再可能触发
kind 门禁——门禁的存在保护的是其余 ~50 个声明了默认绑定的组件及未来新增组件。

**已否决备选：**

1. **`Binding<PhaseEnum>` 模板化（编译期强类型）**——不可行：`default_bindings()` 是共享
   `ComponentBase` 上的 virtual，模板化意味着按 workflow fork 整个组件层次，毁掉跨 workflow
   组件复用（Thermo/recorder/comm 组件靠同一基类 + `.on()` 跨 MD/Minimize 使用）。因此验收
   目标定为"注册时明确报错"，而非编译失败。
2. **phase 枚举数值分段（MD 0-99、Minimizer 100-199…）**——仍是裸 int + magic offset，破坏
   `Plan` 的 `[0, COUNT)` 槽位索引（`plan.h:80-82` 按 phase 直接下标），且旧值静默合法。
3. **`PhaseKind` 定义在 `facade/stage.h`**——`core/binding.h` 需反向包含 facade，层次倒置
   （见上）。
4. **在 `PlanBuilder::build()` 校验**——报错点太晚（`build()` 时才炸，离误绑的 `add()` 调用
   行很远），且 `PlanBuilder` 只见 `int`，判别信息在 `add()` 之后已被擦除。

## Files to create or modify

**(A) 判别机制**

- include/atomiverse/core/binding.h
- include/atomiverse/facade/stage.h
- include/atomiverse/facade/driver.h

**(B) 违约组件与其裸 `add()` 调用点**

- include/atomiverse/cpu/semiempirical/am1_evaluator.h
- include/atomiverse/cuda/scf/scf_evaluator.cuh
- tests/cpu/semiempirical/test_am1_scaffold.cpp
- tests/cuda/scf/test_scf.cu
- tests/cuda/scf/test_rhf.cu

**(C) 机械加 tag（`default_bindings()` 改用 `make_binding`）**

CPU 头（34）：
- include/atomiverse/cpu/md/init.h
- include/atomiverse/cpu/md/velocity_verlet.h
- include/atomiverse/cpu/md/nve_limit.h
- include/atomiverse/cpu/md/fix_langevin.h
- include/atomiverse/cpu/md/fix_csvr.h
- include/atomiverse/cpu/md/nose_hoover.h
- include/atomiverse/cpu/md/nose_hoover_npt.h
- include/atomiverse/cpu/md/barostat.h
- include/atomiverse/cpu/md/deform.h
- include/atomiverse/cpu/comm/reduce.h
- include/atomiverse/cpu/comm/atom_migration.h
- include/atomiverse/cpu/comm/halo_exchange.h
- include/atomiverse/cpu/comm/rcb.h
- include/atomiverse/cpu/kspace/ewald.h
- include/atomiverse/cpu/kspace/pme.h
- include/atomiverse/cpu/observer/virial_stress.h
- include/atomiverse/cpu/async/async_observer.h
- include/atomiverse/cpu/async/snapshot_publisher.h
- include/atomiverse/cpu/minimizer/lbfgs.h
- include/atomiverse/cpu/minimizer/cg.h
- include/atomiverse/cpu/minimizer/opt_log.h
- include/atomiverse/cpu/minimizer/trajectory_writer.h
- include/atomiverse/cpu/minimizer/final_geometry.h
- include/atomiverse/cpu/minimizer/async_observer.h
- include/atomiverse/cpu/vib/reference_geometry.h
- include/atomiverse/cpu/vib/mass_weighting.h
- include/atomiverse/cpu/vib/displace_step.h
- include/atomiverse/cpu/vib/hessian_accumulator.h
- include/atomiverse/cpu/vib/eigen_solver.h
- include/atomiverse/cpu/vib/symmetrize.h
- include/atomiverse/cpu/vib/thermochem.h
- include/atomiverse/cpu/vib/normal_mode_writer.h
- include/atomiverse/cpu/vib/frequency_report.h
- include/atomiverse/cpu/mc/translate_move.h

CUDA 头（12）+ 类外定义（3）：
- include/atomiverse/cuda/md/velocity_verlet.cuh
- include/atomiverse/cuda/md/nve_limit.cuh
- include/atomiverse/cuda/md/fix_langevin.cuh
- include/atomiverse/cuda/md/nose_hoover.cuh
- include/atomiverse/cuda/md/nose_hoover_npt.cuh
- include/atomiverse/cuda/md/barostat.cuh
- include/atomiverse/cuda/md/deform.cuh
- include/atomiverse/cuda/md/init.cuh
- include/atomiverse/cuda/md/kinetic_energy.cuh
- include/atomiverse/cuda/minimizer/sd.cuh
- include/atomiverse/cuda/minimizer/cg.cuh
- include/atomiverse/cuda/minimizer/lbfgs.cuh
- src/cuda/md/barostat.cu
- src/cuda/md/deform.cu
- src/cuda/md/nose_hoover.cu

测试/工具/基准内的本地组件（8）：
- tests/cpu/driver_test_helpers.h
- tests/cpu/test_snapshot_publisher.cpp
- tests/cuda/async/test_async_bridge_tearing.cu
- tests/cuda/async/test_hot_path_no_sync.cu
- tools/debug_amorphous/run_atomiverse.cpp
- tools/debug_amorphous/prof_amorphous.cpp
- tools/debug_amorphous/gather_frame.h
- benchmarks/md/pe/performance/polyethylene_tension.cpp

**(D) 文档 / doc-comment 对齐（滞后示例 `eval.add(SCFEvaluator)` → 显式绑定）**

- include/atomiverse/cpu/driver/evaluate.h
- include/atomiverse/cpu/scf/scf_solver.h
- include/atomiverse/cuda/driver/evaluate.cuh
- include/atomiverse/cuda/scf/scf_solver.cuh
- include/atomiverse/cuda/observer/mulliken.cuh
- include/atomiverse/cuda/recorder/density_recorder.cuh
- docs/user-guide/backend-switching.md
- docs/developer/extending-cpu-subsystems.md

**(E) 新测试**

- tests/cpu/test_driver_binding_kind.cpp (new)
- tests/cuda/test_driver_binding_kind.cu (new)
- tests/CMakeLists.txt

## Tasks

- [x] Test: 编写失败测试 tests/cpu/test_driver_binding_kind.cpp (new) 与
      tests/cuda/test_driver_binding_kind.cu (new)，注册进 tests/CMakeLists.txt（Files E）——
      覆盖：`cpu::md::MD::add(AM1Evaluator{})` throw；`cuda::md::MD::add(cuda::scf::SCFEvaluator{…})`
      throw；`cpu::Evaluate::add(cpu::md::VelocityVerlet{…})` / `cpu::minimizer::Minimize::add(
      cpu::md::VelocityVerlet{…})` / MD + Evaluate-kind 探针组件 → throw 且信息同时含双方
      workflow 名（"MD"/"Evaluate"/"Minimizer"）；同 kind `add()` 行为不变；Evaluate/MD/Minimize/
      Hessian/MC 五类 driver 最小构建后经 `Driver::plan()` 断言 phase 顺序与自动注入组件集不变。
      测试直接使用最终 API（`make_binding`），RED 阶段以编译失败或断言失败呈现。
- [x] Add: 在 include/atomiverse/core/binding.h（Files A）新增 `enum class PhaseKind {MD,
      Minimizer, Evaluate, Vib, Sampling}`、`phase_kind_name(PhaseKind)`，并将 `DefaultBinding`
      改为 `{PhaseKind kind; int phase; int order}`（kind 首位、无默认值，旧式两元初始化编译失败）。
- [x] Add: 在 include/atomiverse/facade/stage.h（Files A）为 5 个 phase namespace 各加
      `constexpr PhaseKind kind_of(Phase)`（沿用 `phase_names` 的 ADL tag-dispatch 惯例），并新增
      `ATV::make_binding(PhaseEnum, int order = 0)` 助手；新增 include `core/binding.h`。
- [x] Refactor: 在 include/atomiverse/facade/driver.h（Files A）的 `Driver::add()` 两个所有权
      分支加入判别校验（私有助手函数）：`db.kind != kind_of(PhaseEnum)` 即 throw
      `std::runtime_error`，信息点名组件名、绑定所属 workflow 与本 driver 类别，并提示
      `.on(phase).add()`；`PendingBinding`/`PhaseBuilder`/`PlanBuilder` 保持不变。
- [x] Refactor: 删除 include/atomiverse/cpu/semiempirical/am1_evaluator.h 与
      include/atomiverse/cuda/scf/scf_evaluator.cuh 的 `default_bindings()` override，对齐
      "Evaluator 显式绑定"约定；将 tests/cpu/semiempirical/test_am1_scaffold.cpp、
      tests/cuda/scf/test_scf.cu、tests/cuda/scf/test_rhf.cu 中的裸 `add()` 调用点迁到
      `.on(evaluate::Phase::evaluate).add(...)`（Files B，一次 commit 迁全）。
- [x] Refactor: 机械迁移 Files C 全部 57 个文件的 `default_bindings()` 返回体到
      `make_binding(<Phase 枚举值>, <order>)`；以全量编译作为完备性门禁（未迁站点必编译失败）；
      测试/工具本地组件若确属跨 driver 复用，同步改为显式 `.on()`。
- [x] Docs: 修正 Files D 全部 8 个文件中滞后的 `eval.add(SCFEvaluator…)` 示例（改
      `.on(evaluate::Phase::evaluate).add()` 或 `single_point`），更新
      docs/developer/extending-cpu-subsystems.md 的 `default_bindings` 示例到 `make_binding`；
      为 `PhaseKind`/`phase_kind_name`/`kind_of`/`make_binding` 补 doxygen `@brief`/`@param`/`@return`。
- [x] Verify: 确认冷路径改造零热路径影响——`git diff --name-only` 不含
      include/atomiverse/facade/plan.h 与 src/cpu/driver/backend.cpp（Plan 的 build 期
      scheduling/frequency 缓存、backend.cpp:264 零开销 timing 路径、step 循环结构均未触碰）；
      现有 test_step_timing / test_timing_report / test_md_plan_equivalence 全绿。
- [x] Run: full check + test suite（`cmake --build build -j$(nproc)` + `ctest --test-dir build
      --output-on-failure`；CUDA 测试按 CLAUDE.md 于 GH200 GPU 节点跑）。

## Testing strategy

**Happy path**

- 同 kind 注册不变：`md.add(VelocityVerlet{…})`、`minimize.add(LBFGS{…})`、
  `hessian.add(DisplaceStep{…})`、`mc.add(TranslateMove{…})` 均成功，`build()` 后
  `Driver::plan()` 中 phase 槽位与 order 与改造前一致。
- 五类 driver（Evaluate / MD / Minimize / Hessian / MC）最小回归：各自最小组件组合构建后，
  断言 phase 顺序与自动注入组件（Migrate/Reduce/Halo/ZeroForces/ZeroEnergy/Publish/Thermo）
  不变；既有 test_md_plan_equivalence、test_md_set_api、test_evaluate_driver、
  test_cpu_driver_role_validation、vib/mc 套件全绿构成回归网。

**误绑报错（本 spec 的核心验收）**

- `cpu::md::MD::add(AM1Evaluator{})` → throw `std::runtime_error`（删除默认绑定后走
  "no default bindings" 拒绝路径；回归点：改造前静默落 `post_integrate`）。
- `cuda::md::MD::add(cuda::scf::SCFEvaluator{2, 0.1, 1e-6, 100})` → 同理 throw（GPU 节点）。
- kind 门禁路径（信息含双方 workflow 名）：`cpu::Evaluate::add(cpu::md::VelocityVerlet{…})` →
  信息同时含 "MD" 与 "Evaluate"；`cpu::minimizer::Minimize::add(cpu::md::VelocityVerlet{…})` →
  含 "MD" 与 "Minimizer"；测试本地 Evaluate-kind 探针组件 `add()` 进 MD → 含 "Evaluate" 与 "MD"。
- lvalue 与 rvalue 两个 `add()` 分支均被上述用例覆盖（各至少一例）。

**Edge cases**

- 空 `default_bindings`（`AsyncComponentBase` 等）仍走既有 "no default bindings" throw，
  test_async_component_base 保持绿。
- 跨 driver 复用组件（comm::Migrate/HaloSend/HaloRecv/mpi::Reduce 进 Minimize）走
  `.on()` 显式路径（minimize.h:120-139），不触发门禁——由现有 MPI/minimizer 套件回归。
- `.on(phase).add()` 路径不做运行期判别（编译期强类型），行为与改造前完全一致。

**链级不变量（arch-refactor 全链共享）**

- clean-break：无兼容层、无 `[[deprecated]]`，一次 commit 迁全部 call site；`kind` 首位无默认值
  使漏迁 = 编译失败。
- 热路径门禁：`facade/plan.h`（build 期缓存 scheduling/frequency）与 `src/cpu/driver/backend.cpp`
  （:264 零开销 timing 路径）零改动，step 循环结构不变；`Driver` 命名保留；方法 = Solver +
  Evaluator、物性 = Capability 的分层不受本 spec 影响（未触碰 role/capability 校验）。
- **为何不需要 PE stage-5 性能门禁**：本 spec 的全部运行期改动位于 `Driver::add()` 与既有
  `PlanBuilder::build()` 之前的注册期——每组件执行一次、发生在 `run()` 之前的冷路径；step
  循环内消费的 `Plan`/`Operation` 数据结构逐字节不变。现有 test_step_timing /
  test_timing_report 即足够的守护。

**Domain validation**：不适用——纯架构 spec，无新物理；物理正确性由未改动的既有单测
（能量守恒、SCF 收敛、PE xcheck 等）回归保证。

## Out of scope

- `ComponentBase` 瘦身与 `HaloParticipant` 抽出（链 02-halo）。
- neighbor 结构、engine 层、Context、Frame、molrs 桥、Python 绑定改造（链 03–08）。
- CLAUDE.md 重写与知识固化（链 09-harness，链尾统一处理；本篇不改 CLAUDE.md）。
- `Binding<PhaseEnum>` 编译期强类型化（已否决，见 Design）。
- phase 枚举本身的增删、改值或数值分段；`PlanBuilder`/`Plan`/`StageSlot` 数据结构变更。
- role/capability 校验逻辑（`EvaluateDriver`/`DynamicsDriver` 的既有验证规则不动）。
- namespace 与目录对齐（CLAUDE.md 已声明为 deliberate follow-up，非本篇顺手改）。
