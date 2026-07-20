---
title: SCF 驱动化重构 · 04a 原生 Driver 节点嵌套（cpu::Driver 擦除基 + PendingBinding 变体，删 SubDriver）
slug: scf-driver-refactor-04a-nest
status: done
chain: scf-driver-refactor
link: 04a/08
created: 2026-07-20
---

# SCF 驱动化重构 · 04a 原生 Driver 节点嵌套

## Summary
把 `SubDriver` 作为"Driver→Component 适配器"的三件事（build 时 `compile()`+`build_plan`、launch 时 `run_nested`、destroy 时 `destroy_plan`）从一个**公开包装类**提升进 Plan/后端的**节点表示**本身：新增非模板多态 `cpu::Driver` 擦除基（由 `CPUDriver<P>` 实现），把 `detail::PendingBinding` 与编译后 `Plan` 的节点由 `Component*` 改为 `variant<Component*, DriverBinding{cpu::Driver*, budget}>`，让 CPU 后端 `step`/`build_components`/`destroy_components` 按节点类分派；新增 `on(phase).add(driver&&/&, order, budget)` 重载（rvalue 经 `owned_drivers_` 拥有、lvalue 借用），budget 默认取子 driver 自己的 `default_budget()`。据此**删除 `cpu/driver/sub_driver.h`**，其嵌套测试迁到原生机制并保留"内层收敛不能停外层"断言。子 verdict 隔离（`run_nested` 里 `ctx.loop` 重指，`cpu_driver.h:432-437`）与跨工作流误绑门（`validate_binding_kinds_`，`driver.h:270-284`）原样保留。**本环不动 SCF**——`MeanFieldMethod` 本环仍经 `HF` 是 Component；把方法-driver 翻到原生节点属 04b。CUDA 不受影响（变体加法式、driver 重载 SFINAE 关闭）。

## Domain basis
本环纯基础设施，无物理；领域锚点是**行为保持**——原生节点嵌套与旧 `SubDriver` **逐位等价**：同一子 driver 在外层每次迭代跑同一预算、写同样的 Frame，内层 verdict 仍只落在内层 `status_`。

设计根据取自本库已固化的架构原则（`sub_driver.h:23-29`）：**类型擦除只发生在节点边界，永不在 stage 边界**——父 loop 走父 Phase 枚举，子 loop 走子 Phase 枚举，二者从不相见，跨界的只是一个"节点"。旧设计让这个节点必须是 `Component`（靠 `SubDriver` 伪装）；本环把节点边界原生化——一个 Driver **就是**一个节点，无需包装。且"内层不能停外层"这一 bug 曾两次上线（LBFGS 无条件写 converged；NEB 首个收敛像停整条带），现由 `run_nested` 的 `ctx.loop` 重指保证结构上不可写，本环必须保住该不变量。

沿用整链决策日志（01–03b 固化）。

## Design

**验证内部（实测行号）：** `pending_` 为 `vector<PendingBinding<Component>>`（`driver.h:247`），`PendingBinding{phase, order, Component* component}`（`:43-47`）；`Driver::build()` 用 `PlanBuilder<Component>` 把每个 `pb.component` 编进 `plan_`（`:219-225`）；组件向探针 `for_each_pending`（`cpu_driver.h:276` `fn(*pb.component)`）、`bind_halo_participants_`（`:303` 直取 `pb.component`）、约束扫描（`:191`）、`validate_long_range_contract_`（`:335`）、`validate_components`（`evaluate.h:107`）皆解引用 `pb.component`；后端有 `build_plan`/`step_plan`/`destroy_plan`（`backend.h:103-112`）；`CPUDriver::compile()`=`build()`（`:634`）、`run_nested`（`:660`）、`setup_pass` 的 `ctx.loop` 重指（`:432-437`）。

**1. `cpu::Driver` 非模板基类（新 `cpu/driver/driver.h`）——用户裁定"不要 DriverNode 这种新概念，就用 Driver 基类"：** `CPUDriver<PhaseEnum>` 是模板，异型容器（Plan 槽）无法直接持有，故必须有一个非模板基——这个基类就命名为 **`cpu::Driver`**，不引入任何新概念名；类型擦除的三件事（`sub_driver.h:94-115`）作为它的纯虚接口：
- `void build_on(cpu::Backend& rt)` = `this->compile(); rt.build_plan(this->plan());`
- `int run_nested_erased(cpu::Context& ctx, cpu::Frame& f, int budget)` = `this->run_nested(ctx, f, budget);`
- `void destroy_on(cpu::Backend& rt)` = `rt.destroy_plan(this->plan());`
- `int default_budget() const` = 子 driver 自配预算（基类虚，preset 覆写；方法-driver 在 04b 覆写返回 `max_iter`；`Minimize` 覆写返回其 `max_iter`）。
`CPUDriver<PhaseEnum>` 额外派生 `cpu::Driver` 并实现之（它已有 `compile`/`run_nested`/`plan`）。基类含虚析构，故 `dynamic_cast<const HF*>(driver_ptr)`（04b 用）良定义。

**2. Plan 节点变体（`facade/driver.h` + `facade/plan.h`）：** `PendingBinding` 与编译后 `Plan` 的节点由裸 `Component*` 改为 `variant<Component*, DriverBinding>`，`DriverBinding{cpu::Driver* driver; int budget}`（driver 节点类经节点特征按设备取：CPU 为 `cpu::Driver`，无嵌套的后端为不支持哨兵）。`Driver::build()`（`:219-225`）遍历 `pending_` 时把两类节点都编进 `PlanBuilder`/`plan_`。所有解引用 `pb.component` 的向组件探针改为节点感知：组件节点照旧 `fn(*component)`，**driver 节点被组件向探针跳过**（`for_each_pending`、halo 绑定、约束扫描、long-range、`validate_components`）——因为容器本身无 Role（旧 `SubDriver` role=0，`sub_driver.h:84-86`），有 Role 的是它所**包含**的东西。另提供一个节点感知遍历（暴露 driver 节点）供 04b 做 HF 跨投射。

**3. `add`/`add_ref`/`on().add` driver 重载（`facade/driver.h` + PhaseBuilder）：** driver 无 `default_bindings()`，故只经 `on(phase)` 绑定（与"无默认绑定的组件必须 `.on(phase).add()`" `:177` 同规）：
- `on(phase).add(CPUDriver<InnerPhase>&& child, int order = 0, int budget = -1)` —— **拥有**：child 移入 `owned_drivers_`（基类新成员 `vector<unique_ptr<cpu::Driver>>`，与 `owned_components_` 并列）；`budget<0` 取 `child.default_budget()`。
- `on(phase).add_ref(CPUDriver<InnerPhase>& child, int order = 0, int budget = -1)` —— **借用**：非拥有，run 后可查 `child.converged()`/`stop_reason()`（与组件 `add_ref` 同借用契约）。
- 裸 `add(driver&&/&)` 删除，报 "drivers have no default bindings — use .on(phase).add()"（镜像组件路径）。
- 这些重载在**无嵌套的后端 SFINAE 关闭**，故 CUDA facade 原样编译、绑不出 driver 节点。

**4. CPU 后端节点分派（`backend.h` + `src/cpu/driver/backend.cpp`）：** 编译后 Plan 的每个节点带类标签；`build_components()` 对 driver 节点调 `driver->build_on(*this)`、组件节点照旧 `component->build`；`step(idx)` 对 driver 节点调 `driver->run_nested_erased(context(), frame, budget)`、组件节点照旧 `component->launch`；`destroy_components()` 对 driver 节点调 `driver->destroy_on(*this)`。`run_nested` 内 `ctx.loop` 重指（`cpu_driver.h:432-437`）与恢复不变，子 verdict 仍只落子 `status_`。

**5. 跨工作流误绑门交互（`validate_binding_kinds_`，`driver.h:270-284`）：** 该门在**注册时**只校验一个**组件**的 `default_bindings()` 各项 `PhaseKind` 是否等于父 driver 的 `kind_of(PhaseEnum{})`——对组件绑定**原样保留、行为不变**。driver 节点**无 `default_bindings()`**，只经 `on(phase)` 用**父的**phase 值绑定，故其槽的 `PhaseKind` 平凡地就是父的、不可能跨工作流错配；**子 driver 自己的 phase 枚举（`scf::Phase`/`minimizer::Phase`/…）完全活在 `run_nested`/子的独立编译 Plan 内，永不触达父的 `validate_binding_kinds_`**——只有槽登记（父 phase）重要，子内部计划不参与父的 kind 检查。

**删除 + 迁移：** 删 `cpu/driver/sub_driver.h`（整类提升进 cpu::Driver/节点分派后再无用）；其嵌套测试迁到原生 API。

**Reuse decision：**
- `reuse` `run_nested` + `ctx.loop` 重指（`cpu_driver.h:432-437,660`）、`build_plan`/`step_plan`/`destroy_plan`（`backend.h:103-112`）、`compile`/`plan`（`:634`）—— 原样复用，不改语义。
- `generalize` `PendingBinding`/`Plan`/`PlanBuilder`（`driver.h:43-47,219-225`，`plan.h`）—— 节点由 `Component*` 泛化为 `variant<Component*, DriverBinding>`；对应"Generalize Plan 节点表示"性质的任务。
- `generalize` 后端 `step`/`build_components`/`destroy_components`（`backend.cpp`）—— 加 driver 节点分派分支。
- `reuse` `SubDriver` 的三件事逻辑（`sub_driver.h:94-115`）—— **提升**进 `cpu::Driver`，而非另写。
- `new` —— `cpu::Driver` 基 + 节点特征 + `owned_drivers_` + PhaseBuilder driver 重载。**无 `make_*`、无新 `.sh` 门。**
- `delete`：`cpu/driver/sub_driver.h`（公开包装类退场）。

## Files to create or modify
- `include/atomiverse/cpu/driver/driver.h` (new) — `cpu::Driver` 擦除基 + 节点特征。
- `include/atomiverse/cpu/driver/cpu_driver.h` — `CPUDriver<P>` 派生并实现 `cpu::Driver`（`build_on`/`run_nested_erased`/`destroy_on`/`default_budget`）；`for_each_pending`（`:276`）与直取 `pb.component` 处（`:191,303,335`）改节点感知（driver 节点跳过）；加节点感知遍历。
- `include/atomiverse/facade/driver.h` — `PendingBinding` 变体；`owned_drivers_`；PhaseBuilder driver `add`/`add_ref` 重载（拥有/借用/budget 默认）；删裸 `add(driver)`；`build()`（`:219-225`）按节点类编 Plan；SFINAE 守卫。
- `include/atomiverse/facade/plan.h` — `Plan`/`PlanBuilder` 节点表示改变体 + 节点类访问器。
- `include/atomiverse/cpu/driver/backend.h` + `src/cpu/driver/backend.cpp` — `step`/`build_components`/`destroy_components` 按节点类分派。
- `include/atomiverse/cpu/driver/sub_driver.h` (delete)。
- `tests/cpu/driver/test_sub_driver.cpp` → `tests/cpu/driver/test_driver_nesting.cpp` (rename+migrate) — 用原生 `on(phase).add(inner, budget)`；保留"内层收敛不能停外层"+ 借用式 run 后 `inner.converged()`/`stop_reason()` 断言。
- `tests/CMakeLists.txt` — 更新该测试注册名。

## Tasks
- [x] Migrate tests/cpu/driver/test_sub_driver.cpp to a failing native-API test tests/cpu/driver/test_driver_nesting.cpp (bind an inner driver via on(phase).add / add_ref, assert it runs its budget each outer iteration, inner convergence cannot stop outer, and borrowed-form post-run inner.converged()/stop_reason()) and update tests/CMakeLists
- [x] Add cpu::Driver erased base + node trait (cpu/driver/driver.h) and make CPUDriver<P> implement build_on/run_nested_erased/destroy_on/default_budget, promoting the sub_driver.h:94-115 ops, with Doxygen briefs
- [x] Make facade PendingBinding + Plan/PlanBuilder nodes a variant<Component*, DriverBinding> and thread both kinds through Driver::build() (facade/driver.h:219-225, plan.h)
- [x] Add on(phase).add(driver&&, order, budget)/add_ref(driver&, order, budget) overloads (own via owned_drivers_, borrow), budget defaulting to child.default_budget(); delete the bare add(driver) with a use-on(phase) message; SFINAE-guard so backends without nesting compile unchanged
- [x] Dispatch node kind in the CPU backend (backend.h/backend.cpp: build_components→build_on, step→run_nested_erased, destroy_components→destroy_on) and make for_each_pending + the direct pending_ probes (cpu_driver.h:191,303,335) node-aware so component probes skip driver-nodes
- [x] Delete cpu/driver/sub_driver.h and the old test_sub_driver.cpp; sweep any SubDriver references
- [x] Run full check + test suite

## Testing strategy
- **原生嵌套（单元，秒级）**：`test_driver_nesting` 用 `cpu::minimizer::Minimize` 作内层、`cpu::md::MD` 作外层，`outer.on(md::Phase::end_of_step).add_ref(inner, /*budget=*/50)`（及 rvalue `add(...)` 拥有式），断言：内层每个外层迭代跑到其预算/收敛；外层跑满自己的步数不被内层收敛打断（**内层收敛不能停外层**，复用 `ctx.loop` 重指）；借用式 run 后 `inner.converged()`/`stop_reason()` 可读。
- **节点分派**：driver 节点经 `build_on`（一次）/`run_nested_erased`（每外层迭代）/`destroy_on` 走通；组件节点行为不变（既有 `test_evaluate_driver`/`test_md`/`test_minimize_*` 全绿）。
- **误绑门保持**：组件跨工作流误绑仍在 `add` 时抛（既有 `test_driver_binding_kind` 绿）；driver 节点只经父 phase 绑定，不入该门。
- **Edge**：同一内层被拥有（rvalue）vs 借用（lvalue）；budget 显式 vs 默认 `default_budget()`；内层 phase 枚举（minimizer::Phase）与外层（md::Phase）数值重叠但互不可见。
- **Domain validation**：与旧 `SubDriver` 逐位等价即本环领域验证。
- **CUDA**：facade 变体为加法式、driver 重载 SFINAE 关闭，CUDA facade 原样编译、绑不出 driver 节点（CUDA 构建在 GPU 节点，本链 CPU-only，验证见 Out of scope）。
- **SCF 不动**：`MeanFieldMethod` 本环仍 Component，既有 SCF/AIMD 绑定路径不变。
- **Regression example**：本环**不**含 `regressions/` 案例——纯基础设施，无端到端新功能；本链回归在 03b。runtime 关口是完整 ctest 快套件。
- **RED/GREEN**：唯一新增可测符号是原生嵌套，由 `test_driver_nesting`（RED，先于基础设施）覆盖。

## Out of scope
- HF 收窄纯 tag、删 `MeanFieldMethod::build`/`launch`、ASPC 重宿、删 `scf_evaluator`、`validate_components` 对 driver 节点跨投射 HF、`SCFResultHandle`→`MeanFieldOutputs` 收窄、示例/测试迁移 —— 属 04b。
- CUDA driver 节点分派/嵌套 —— CUDA 保持无嵌套（重载 SFINAE 关闭）；CUDA 构建验证在 GPU 节点，本链 CPU-only。
- `ATV_BUILD_RHF→ATV_BUILD_HF` 改名 —— 05；Handler→Criterion 全树迁移 + 删 `LoopHandler` —— 06。
