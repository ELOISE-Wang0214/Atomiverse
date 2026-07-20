---
title: SCF 驱动化重构 · 04b 外层接线（HF 纯 tag，删 build/launch，MeanFieldOutputs 收窄，run/launch 债务终结）
slug: scf-driver-refactor-04b-wire
status: code-complete
chain: scf-driver-refactor
link: 04b/08
created: 2026-07-20
---

# SCF 驱动化重构 · 04b 外层接线

## Summary
本环把高斯 SCF 方法翻到 04a 的原生 Driver 节点嵌套上，并**终结 03a 起显式记账的 run/launch 共存债务**：`HF` 变纯多态 tag（虚析构，不再派生 `Component`）；`MeanFieldMethod<D> : CPUDriver<scf::Phase>, HF`（不再是 Component），覆写 `default_budget()→max_iter`；**删除 `build()`/`launch()`**（单一入口 `run` 留存），ASPC 预测/校正重宿进 `drive_scf_`、只凭持久成员 `use_predictor_`（去掉 03a 的组件路径门）；**删除 `scf_evaluator.{h,cpp}`**（`using SCF = RKS` 留在 `methods.h`）。外层 `validate_components`（`evaluate.h/md.h/vib/minimize`）改走 04a 的节点感知遍历，把 driver 节点跨投射到 kind tag（`HF` 等）；绑定形态 `outer.on(<force/evaluate 相>).add(method[, budget=max_iter])` **源码稳定**（靠基类推导自动从组件 `add` 重路由到 driver `add`）。`SCFResultHandle` 收窄为 `MeanFieldOutputs{result, shells, ready}`（`IntegralPack` 彻底内部化），迁移 `orbital_info_observer` + 密度/轨道 cube writer（mulliken/density observer 经核实不消费 handle，不动）。迁移示例 05/06/11/12/16 + `single_point` 文档 + 相关测试。用户可见能量/行为不变。

## Domain basis
本环纯架构接线，无物理；领域锚点是**行为保持**：HF-tag 翻转 + 原生绑定 + handle 收窄不改任何数值——AIMD 轨迹（`06_aimd`）、单点能量、cube 密度逐位不变。

据已固化的目标形态：**方法是 Driver、不是 Component；能量源是 kind（经跨投射识别）；Frame 是唯一货币；`IntegralPack` 是私有工作区、绝不上公共面**。本环让这三条对高斯 SCF 全部落地——`run` 与 `launch` 从此不再共存于一类（该债务自 03a 显式记账、指向本环），`SCFResultHandle` 里最后一个内部类型 `IntegralPack` 被收窄成消费方真正需要的 `shells`。

沿用整链决策日志（01–04a 固化）；新用户裁定（Handler→Criterion）在 06 落地，本环不涉。

## Design

**1. `HF` 纯多态 tag（`cpu/component/hf.h`）：** `class HF { public: virtual ~HF() = default; };`——无 `Component` 基，含虚析构故 `dynamic_cast<const HF*>(driver_ptr)` 良定义。`MeanFieldMethod<D> : public CPUDriver<scf::Phase>, public HF`（不再 Component），覆写 `cpu::Driver::default_budget()` 返回 `max_iter`；每方法 `name()`（rhf/uhf/rks/uks）从 Component override 重宿为 driver 级 `name()`（供 plan/timing 标签）。

**2. 删 `build()`/`launch()` + ASPC 重宿（`mean_field_method.{h,cpp}`）：** `MeanFieldMethod::build(BuildContext)`/`launch(LaunchContext)`（03a 暂留的组件路径）**删除**——`run(frame)`（→私有成员 `drive_scf_`）成为唯一入口。ASPC 逻辑（旧 `methods.cpp:259-310` 在 launch 里）重宿进 `drive_scf_`：前置按 `use_predictor_ && 有历史` 预测→喂初始猜测，后置校正+记录，**只凭持久成员**（`aspc_alpha_`/`aspc_beta_`/`last_P_*`），**去掉 03a 的"组件路径"限定**。方法作为嵌套 driver 节点在 MD 里跑时，每个外层力相迭代经 `run_nested→run→drive_scf_`，方法对象跨外层步持久，故 ASPC 历史照旧跨步累积——与旧 launch 路径逐位等价。

**3. `scf_evaluator.{h,cpp}` 删除：** `scf_evaluator.h` 仅是兼容 include（`SCF=RKS` 别名实在 `methods.h:267`），删之无损别名；`scf_evaluator.cpp` 从 root CMake 移除。测试：`test_scf_evaluator`（经组件 build/launch 用 `SCF`）**迁移**到 driver 路径（`single_point(frame, rks)` / `Evaluate + on().add(rks)` / `rks.run(frame)`）；`test_scf_evaluator_no_alloc`（旧断言 launch 无每步分配）**重宿**为 `test_method_driver_no_alloc`——断言方法作为嵌套节点时 `build_on` 一次性分配、`run_nested`（每外层步）不新分配固定缓冲（SCF 循环复用方法的成员状态矩阵，02 各 stage launch 无分配；DIIS deque 增长为既有文档化例外）。这保住了 AIMD 关键不变量（每 MD 步 SCF 力评估不重建整套机器）。

**4. `validate_components` 节点感知 + 跨投射（`evaluate.h:105`、`md.h:359`、`vib/hessian`、`minimize`）：** 这些 validate 从组件专用 `for_each_pending` 切到 04a 的**节点感知遍历**（暴露 driver 节点）：组件节点照旧 `dynamic_cast<const HF*/Potential*/…>(&component)`，**driver 节点跨投射其多态基到 kind tag**（`dynamic_cast<const HF*>(driver_ptr)`），故方法-driver 作为 driver 节点在 force/evaluate 相绑定时被识别为能量源。`evaluate.h:108-116` 的 kind 列表照旧含 `HF`。

**5. 绑定形态源码稳定（载重简化）：** 外层 `outer.on(<phase>).add(method[, budget])` **无需改调用点**。当 `method`（RHF）本是 Component（至 04a）时，`add` 走组件重载；`method` 变 driver（本环去掉 Component 基）后，组件 `add` 因 SFINAE（约束 `Comp:Component`）排除、driver `add`（约束 `Comp:cpu::Driver`，04a 引入）经**派生→基类模板推导**（`RHF` 唯一派生 `CPUDriver<scf::Phase>` ⇒ `InnerPhase=scf::Phase`）选中，budget 缺省取 `default_budget()=max_iter`。故 `06_aimd:111` 的 `sim.on(md::Phase::force).add(cpu::scf::RHF{}.predictor(3))` 一字不改地从组件绑定重路由为 driver 节点绑定。示例迁移主要是注释（如 `06_aimd:57` "SCF::launch reassembles" 已不成立）。

**6. `SCFResultHandle`→`MeanFieldOutputs`（`mean_field_method.h`）：** 收窄为消费方真正所需：
```
struct MeanFieldOutputs {                    // 原 SCFResultHandle 重命名
    std::shared_ptr<const SCFResult> result;
    std::shared_ptr<const std::vector<Shell>> shells;  // 基组，取代 IntegralPack
    std::shared_ptr<const bool> ready;
};
```
去掉 `IntegralPack integrals` 指针 ⇒ `IntegralPack` 彻底不上任何公共面（留纯内部）。`result_handle()` 改填 `shells`（方法后置从 `integrals_sp_->shells` 拷贝共享）。消费方（**核实仅**）：`orbital_info_observer.{h,cpp}`、`recorder/density_cube_writer.{h,cpp}`（`:44,46` `integrals.shells`→`*handle.shells`）、`recorder/orbital_cube_writer.{h,cpp}`（`:46`）。`observer/mulliken`、`observer/density` **经 grep 核实不消费 handle**，不动。

**7. `single_point`（`single_point.h`）去留：** **保留、几乎不改**。它是模板自由函数（`Evaluate + on(phase).add_ref(evaluator)`），对方法-driver 经 04a 的 `add_ref(CPUDriver<InnerPhase>&, order, budget=-1)` 重载工作（budget=`default_budget()=max_iter`），对组件求值器（`PWEvaluator`/`AM1`）经组件 `add_ref` 工作——同一模板两条重载各归其位，无歧义（无类型同时是 Component 与 Driver）。仅更新其 Doxygen：方法-driver 作为嵌套节点绑定、预算 `max_iter`。不设 `SinglePoint` 类（禁），不改签名。

**Reuse decision：**
- `reuse` 04a 的原生嵌套（`on().add`/`add_ref` driver 重载、`cpu::Driver`、节点感知遍历）—— 方法-driver 直接绑定，不新增机制。
- `reuse` `drive_scf_`（03a）—— ASPC 只是并入其前后置，不新写循环。
- `reuse` `single_point` 模板 —— 靠 04a 重载自动适配，几乎零改。
- `generalize` `validate_components`（evaluate/md/vib/minimize）—— 从组件专用改节点感知 + driver 跨投射。
- `generalize`/rename `SCFResultHandle`→`MeanFieldOutputs` —— 字段从 `IntegralPack` 收窄到 `shells`。
- `delete`：`MeanFieldMethod::build`/`launch`、`scf_evaluator.{h,cpp}`、`SCFResultHandle` 旧形。
- `new` —— 无新自由函数、无 `make_*`、无新 `.sh` 门；仅类型/绑定改写 + 一个 no-alloc 重宿测试。

## Files to create or modify
- `include/atomiverse/cpu/component/hf.h` — `HF` 纯多态 tag（删 `: public Component`，加虚析构）。
- `include/atomiverse/cpu/scf/mean_field_method.h` — `MeanFieldMethod<D> : CPUDriver<scf::Phase>, HF`；删 `build`/`launch` 声明；覆写 `default_budget()`/`name()`；`SCFResultHandle`→`MeanFieldOutputs{result, shells, ready}`。
- `src/cpu/scf/mean_field_method.cpp` — 删 `build`/`launch` 定义；ASPC 并入 `drive_scf_`（去组件路径门）；`result_handle()` 填 `shells`。
- `include/atomiverse/cpu/driver/evaluate.h` — `validate_components` 节点感知 + driver 跨投射 HF（`:105-119`）。
- `include/atomiverse/cpu/driver/md.h` — 同（`:359`）。
- `include/atomiverse/cpu/vib/hessian.h` — validate 节点感知 + force 相绑定方法-driver。
- `include/atomiverse/cpu/minimizer/minimize.h` — validate 节点感知（几何优化力源为 driver 节点）。
- `include/atomiverse/cpu/scf/orbital_info_observer.{h,cpp}` — `MeanFieldOutputs`。
- `include/atomiverse/cpu/recorder/density_cube_writer.{h,cpp}` + `orbital_cube_writer.{h,cpp}` — `*handle.shells` 取代 `integrals.shells`。
- `include/atomiverse/cpu/driver/single_point.h` — Doxygen 更新（方法-driver 嵌套绑定，budget=max_iter）。
- `include/atomiverse/cpu/scf/scf_evaluator.h` (delete) + `src/cpu/scf/scf_evaluator.cpp` (delete)。
- `CMakeLists.txt` (root) — 去掉 `scf_evaluator.cpp`。
- 示例：`examples/{05_geom_opt,06_aimd,16_density_cube,11_koopmans,12_delta_scf}.cpp` — 注释/验证 + 16 的 handle 重命名。
- 测试：`tests/cpu/scf/test_scf_evaluator.cpp`（迁移→driver 路径）、`test_scf_evaluator_no_alloc.cpp`（重宿→`test_method_driver_no_alloc`）、`tests/cpu/recorder/test_density_cube.cpp`、`tests/cpu/scf/test_orbital_info_observer.cpp`（`MeanFieldOutputs`）、`tests/cpu/test_evaluate_driver.cpp`、`tests/cpu/test_minimize_driver_swap_optimizer.cpp`、minimizer 测试、vib hessian 测试（driver 节点绑定）。
- `tests/cpu/scf/test_method_driver_wire.cpp` (new) — 失败优先：driver 节点作 HF 被 validate 接受 + 原生绑定 + no-alloc + `MeanFieldOutputs`。
- `tests/CMakeLists.txt` — 注册更新。

## Tasks
- [x] Write failing tests (tests/cpu/scf/test_method_driver_wire.cpp) that bind a method-driver via on(force/evaluate).add + single_point and assert validate accepts it as an HF node, MeanFieldOutputs exposes shells, and the nested run allocates nothing per outer step; update tests/CMakeLists
- [x] Make HF a pure polymorphic tag (cpu/component/hf.h, virtual dtor, no Component) and MeanFieldMethod<D> derive CPUDriver<scf::Phase> + HF (not Component), overriding default_budget()→max_iter and name()
- [x] Delete MeanFieldMethod build()/launch() and re-host ASPC predict/correct/record into drive_scf_ keyed only on use_predictor_ (drop the 03a component-path gate)
- [x] Delete scf_evaluator.{h,cpp} (keep using SCF = RKS in methods.h; update root CMakeLists) and migrate test_scf_evaluator to the driver path + re-home test_scf_evaluator_no_alloc onto the nested method-driver run
- [x] Switch validate_components (evaluate.h, md.h, vib/hessian, minimize) to the 04a node-aware walk and cross-cast driver-nodes to the kind tags (HF etc.)
- [x] Narrow SCFResultHandle → MeanFieldOutputs{result, shells, ready} and migrate orbital_info_observer + density/orbital cube writers + their tests + example 16
- [x] Migrate examples 05_geom_opt/06_aimd/11_koopmans/12_delta_scf + the single_point doc (binding shape is source-stable; update stale SCF::launch comments)
- [x] Migrate tests test_evaluate_driver, test_minimize_driver_swap_optimizer, minimizer tests, and vib hessian tests to the driver-node binding
- [x] Run full check + test suite

## Testing strategy
- **原生绑定为 HF（单元，秒级）**：`test_method_driver_wire` 断言 `Evaluate`/`MD`/`Minimize` 的 `validate_components` 接受一个 driver 节点（方法-driver）作 HF 能量源（跨投射）；`eval.on(evaluate::Phase::evaluate).add(RHF{})` 与 `single_point(frame, rhf)` 均绑定为嵌套节点、跑出与旧组件路径一致的能量。
- **两入口/一致**：`rhf.run(frame)`、`single_point(frame, rhf)`、`Evaluate + on().add(rhf)` 三条路径能量一致（共用 `drive_scf_`）。
- **no-alloc 重宿**：`test_method_driver_no_alloc` 断言方法作嵌套节点时 `build_on` 一次性分配、后续每外层步 `run_nested` 不新分配固定缓冲。
- **ASPC 不变**：`.predictor(3)` 的 AIMD（示例 06 规模）轨迹与旧 launch 路径逐位等价。
- **handle 收窄**：`test_density_cube`/`test_orbital_info_observer`/示例 16 编译+通过，经 `*handle.shells` 得到与旧 `integrals.shells` 相同的基组求值；`IntegralPack` 不在任何公共 handle。
- **债务终结（grep）**：`MeanFieldMethod` 无 `build`/`launch` 覆写；`HF` 不派生任何东西；无 `scf_evaluator` token；无 `SCFResultHandle` 旧形。
- **Edge**：方法-driver 拥有（rvalue `add`）vs 借用（`add_ref`）；`with_gradient(false)` 单点 vs `true` 力；`Minimize` 力源为 driver 节点；`vib::Hessian` 力相绑定方法-driver。
- **Domain validation**：AIMD/单点/几何优化/cube 数值不变即领域验证；文献值 e2e 在 03b。
- **Regression example**：本环**不**含 `regressions/` 案例——纯接线；本链回归在 03b。runtime 关口是完整 ctest 快套件。

## Out of scope
- `ATV_BUILD_RHF→ATV_BUILD_HF` 改名 —— 05。
- Handler→Criterion 全树迁移 + 删 `LoopHandler`/`add_handler` —— 06（本环 SCF 已用 02/03a 落的 `ScfConverged` Criterion，不涉外层 handler 迁移）。
- CUDA driver 节点/嵌套 —— CUDA 保持无嵌套；CUDA SCF 保持 inner-solver。
- PW/GPW/NDDO 驱动化 —— 各域自治，本链不动（其求值器仍是 Component）。
