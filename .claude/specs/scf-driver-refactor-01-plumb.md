---
title: SCF 驱动化重构 · 01 门面管道铺设（scf::Phase / PhaseKind::SCF / PhaseTraits / Criterion）
slug: scf-driver-refactor-01-plumb
status: code-complete
created: 2026-07-20
revised: 2026-07-20
chain: scf-driver-refactor
link: 01/05
---

# SCF 驱动化重构 · 01 门面管道铺设

## Summary
本子规格只做一件事：为后续把高斯均场 SCF 自洽循环改造成 `CPUDriver<scf::Phase>` 驱动，把所有**纯增量**的门面/核心/基础设施类型一次性铺好，并删除唯一会阻挡它落地的硬性 grep 门（`tests/cuda/scf_phase_scrub.sh`）。落地后代码树新增了 `ATV::scf` 相位枚举、`PhaseKind::SCF` 判别值、`PhaseTraits<scf::Phase>` 特化、以及 `Criterion` 组件类别 + driver 拉取式轮询，但**没有任何消费方**——现有 `SCFSolver`/`ScfLoop` 仍是唯一驱动 SCF 的路径，四个方法 RHF/UHF/RKS/UKS 的能量逐位不变。用户可见行为零变化；这一环的价值是让 02（阶段组件）、03（方法驱动）能在一个已经就绪的类型底座上开工。

## Domain basis
本环无新物理、无方程。领域锚点是**行为保持**：所有改动都是加类型或删门，不触碰任何数值代码路径，因此 SCF 能量在定义上逐位不变（真正的"迁移前后 RHF/UHF/RKS/UKS 能量 ≤1e-12 Ha 一致"验收与文献值回归案例落在 03，那才是第一次真正走新驱动跑出能量）。

本环遵循并固化贯穿整条链的架构决策日志（供 02–05 引用）：
- ① 一切驱动循环皆 Driver、可组件化嵌套 —— 故 SCF 需要自己的相位枚举 `scf::Phase`，`stage.h` 里"故意没有 scf::Phase"的旧结论就此反转。
- ② Solver = 封闭数值内核的命名约定（EigenSolver/Davidson/FFT 类），判据是"中间步骤是否为可组合的科学阶段"；SCF 有可组合阶段 ⇒ 必须是 Driver，不是 Solver。
- ③（用户修订）Driver 内部状态通道 = **构造器直连引用**：SCF 迭代态（F/P/C/S/X…）是方法 driver 的**私有成员**（真名，如 `F_alpha_`/`P_alpha_`），方法构造器把引用直接接进各 stage 组件——与"参数 ctor 直传"同一条纪律推广到状态。**不设 `ScfWorkspace` 中间袋、不设 `ctx.scf` 切片**；跨 Driver 边界仍 Frame-only，重资产观测仍走具名 handle。
- ④ 无共享 SCF 驱动基类；`MeanFieldMethod<D>` CRTP 保留为参数家 + 唯一 `run(frame)` 入口（改造在 03），本环不引入任何 `SCFDriver` 之名。
- ⑧ 方法裸名沿用树内先例（MD/MC/Minimize/Evaluate）；相位枚举命名沿用 `evaluate::Phase` 先例。
- ⑨ 终止判据也是 Component：新增 `Criterion` 组件类别（`verdict()` 查询、driver 拉取合成），独立的 Handler 顶层概念退场——本环仅增量引入类别，既有 `LoopHandler` 全树迁移与删除属 06-criterion-unify。

## Design
新增五处**纯增量**类型/字段，外加一处删除，彼此无运行期耦合，合起来构成 SCF-as-Driver 的类型底座：

1. **`ATV::scf::Phase`（+ `kind_of`）** —— 在 `facade/stage.h` 里，用一个新的 `namespace ATV::scf { enum class Phase { fock, mix, diagonalize, density, check, COUNT }; }` 替换第 132–137 行"deliberately no scf::Phase"的 NOTE 段。相位序反映一次自洽迭代内的组合：构建 Fock/哈密顿 → DIIS 混合（**Fock 空间**外推，必须先于对角化——现实现即如此：`scf_solver.cpp:241→252→284` 建 F → 对 F 做 DIIS 外推 + level shift → 对角化；Pulay, *Chem. Phys. Lett.* **73**, 393 (1980)；*J. Comp. Chem.* **3**, 556 (1982)）→ 对角化 → 构建密度 → 收敛判定。`COUNT` 作哨兵供 `slot_count`。同文件按 `evaluate::Phase` 的 `kind_of`（stage.h:122–128）逐字仿写一个 `constexpr PhaseKind kind_of(Phase) { return PhaseKind::SCF; }`（ADL 标签派发）。这样 `ATV::slot_of` / `slot_count` / `make_binding` 三个既有入口对 `scf::Phase` 立即可用，无需任何 `static_cast<int>`。

2. **`PhaseKind::SCF`** —— 在 `core/binding.h:14` 的 `enum class PhaseKind { MD, Minimizer, Evaluate, Vib, Sampling, Path }` 追加 `SCF`，并在 `phase_kind_name`（:19 switch）补一个返回 `"SCF"` 的 case。这让 03 里 `make_binding(scf::Phase::…, order)` 携带正确判别值，04 里外层 `validate_binding_kinds_` 能识别 SCF 工作流。

3. **（已撤销——用户裁定"是什么就用什么命名"）** 原计划给 `facade/loop.h` 的 `LoopMetrics` 增 `delta_e`/`diis_error`：**撤销**。SCF 不使用 LoopMetrics；`ScfConverged`（02）自持 `prev_energy_`/`delta_energy_`/`diis_error_` 真名成员，从构造器接线的引用自己测量、自己判定。`facade/loop.h` 本环**不动**（既有 LoopMetrics 的诚实更名 `ConvergenceMetrics` 属 06）。

4. **（已撤销——用户裁定"ScfWorkspace 无存在必要"）** 原计划的 `ScfServices` 切片 + `ctx.scf` + `ScfWorkspace` 前置声明：**撤销**。SCF 迭代态是方法 driver 的私有成员，经**构造器引用直连**流入各 stage 组件（见决策③修订）。`cpu/component/context.h` 本环**不动**。

5. **`PhaseTraits<scf::Phase>` 特化** —— 在 `cpu/driver/cpu_driver.h`（紧邻 :105–155 的既有特化，遵从 :145–149 的警告注释）加一个特化：`eval_phase = fock`（能量/哈密顿产出槽，类比 force 槽）、`neighbor_phase = fock`（SCF 无邻居表，折叠进 fock 以便分类 setter 实例化）、`end_phase = check`。没有它，`CPUDriver<scf::Phase>` 在 03 无法实例化——这正是它作为独立步骤存在的原因。

6. **删除 scrub 门（硬阻塞）** —— `tests/cuda/scf_phase_scrub.sh` 会在 `include/`/`src/`/`examples/`/`tests/` 里 grep 到任何 `scf::Phase` token 即 fail build（仅豁免 `facade/stage.h`）。第 1 条一旦引入枚举，该门必 fail，故必须在**同一子规格**内删除该脚本及其 `test_scf_phase_scrub` ctest 注册（`tests/CMakeLists.txt:592–599 含注释块`），不再设任何替代 `.sh` 门。顺带清扫两处对已删脚本的悬挂注释引用（`tests/cuda/scf_solver_isolation.sh:9`、`include/atomiverse/cuda/pw/pw_scf.cuh:14`）。

7. **`Criterion` 组件类别（终止判据即组件，决策⑨）** —— 新增 `cpu/component/criterion.h`：`class Criterion : public Component`，带纯虚 `ATV::Verdict verdict(const cpu::LoopState&)` 查询接口。`CPUDriver` 增 `add_criterion(rvalue)`/`add_criterion_ref(lvalue)` 注册（复数、可组合，遵循既有 add/add_ref 所有权规则），并在每迭代结束后**拉取**各已注册 Criterion 的 `verdict()`，按严重度 Stop > Converged > Continue 合成、与既有 handler 轮询共用同一次 MPI 归约。裁决通道保持"driver 拉取"——组件依旧写不了循环状态（`ctx.loop` 借用语义不变）。本环纯增量：既有 `LoopHandler`/`add_handler` 原样保留，全树迁移与删除属 06-criterion-unify。**注册契约（spec-audit 澄清）：** `add_criterion(rvalue)` 拥有该对象**一次**（存入 `owned` 存储），并做两件事——把指针放进 `criteria_` 轮询列表；若该 Criterion 还声明了 phase 绑定（`default_bindings()` 非空，如 `ScfConverged` 绑 `check` 相做测量），同一指针**同时**进入 Plan 的相位槽。一次所有权、两处登记；verdict-only 的 Criterion（`default_bindings()` 为空）只进轮询列表。

**Reuse decision**（逐一了结 librarian_report 的复用候选）：
- `reuse` `CPUDriver<PhaseEnum>` 全套设施（cpu_driver.h 的 run_nested/poll_handlers_/setup·step·teardown/metrics/compile）—— 本环不复制任何循环代码，只补 `PhaseTraits<scf::Phase>` 特化，让 03 能原样继承这套机器。
- `generalize` `PhaseKind` + `phase_kind_name`（core/binding.h）—— 追加 `SCF` 值 + name case，扩展既有枚举/switch 覆盖 SCF 工作流；对应一条 "Generalize `PhaseKind` …" 任务。
- `reuse` 阶段类型化入口（`slot_of`/`slot_count`/`make_binding`，stage.h:300/308/331）—— `scf::Phase` + `kind_of` 逐字仿 `evaluate::Phase`（stage.h:122–128），不新增任何转换路径。
- `new — ATV::scf::Phase`：确无既有 SCF 相位枚举（历史那个已随 pw 迁移删除），按 `evaluate::Phase` 形状新造。
- `new — PhaseTraits<scf::Phase> 特化`：既有模板的新特化，逐字仿 md/minimizer 特化。
- 五个阶段组件、`ScfConverged`（librarian 归到 cpu/scf/ 域目录）—— 本环**不**创建；定义属 02（无 Workspace、无切片——状态经构造器引用直连）。

无任何新符号复现 `reuse`/`generalize` 候选。

## Files to create or modify
- `include/atomiverse/facade/stage.h` — 新增 `namespace ATV::scf` 的 `Phase` 枚举 + `kind_of`，替换 132–137 行的 NOTE 为描述新枚举的 Doxygen 段。
- `include/atomiverse/core/binding.h` — `PhaseKind` 追加 `SCF`；`phase_kind_name` 追加 `"SCF"` case。
- `include/atomiverse/cpu/driver/cpu_driver.h` — 追加 `PhaseTraits<ATV::scf::Phase>` 特化；增 `add_criterion`/`add_criterion_ref` 注册与每迭代 Criterion 轮询（严重度合成 + 单次 MPI 归约，与既有 handler 轮询并存）。
- `include/atomiverse/cpu/component/criterion.h` (new) — `Criterion` 组件类别（`Verdict verdict(const LoopState&)` 查询接口）。
- `tests/cpu/driver/test_scf_phase_plumbing.cpp` (new) — 覆盖以上全部管道的失败优先单测。
- `tests/CMakeLists.txt` — 删除 `test_scf_phase_scrub` 的 `add_test` + 其注释块（592–599）；注册 `test_scf_phase_plumbing`。
- `tests/cuda/scf_phase_scrub.sh` (delete) — 硬阻塞 grep 门，连同 ctest 注册一并删除，不设替代。
- `tests/cuda/scf_solver_isolation.sh` — 删去第 9 行 "Mirrors tests/cuda/scf_phase_scrub.sh" 悬挂注释引用。
- `include/atomiverse/cuda/pw/pw_scf.cuh` — 删去第 14 行对已删 scrub 的注释引用（纯注释，不入 CPU 编译）。

## Tasks
- [x] Write failing tests for scf::Phase plumbing (tests/cpu/driver/test_scf_phase_plumbing.cpp) and register it in tests/CMakeLists.txt
- [x] Add ATV::scf::Phase enum + kind_of in include/atomiverse/facade/stage.h, replacing the 132-137 NOTE
- [x] Generalize PhaseKind with SCF value + phase_kind_name case in include/atomiverse/core/binding.h
- [x] Specialize PhaseTraits<ATV::scf::Phase> (fock/fock/check) in include/atomiverse/cpu/driver/cpu_driver.h
- [x] Add the Criterion component kind (Verdict verdict(const LoopState&)) in include/atomiverse/cpu/component/criterion.h plus add_criterion/add_criterion_ref registration and per-iteration severity-composed, MPI-reduced polling in cpu_driver.h (additive; existing LoopHandler untouched)
- [x] Delete tests/cuda/scf_phase_scrub.sh + its test_scf_phase_scrub ctest registration and sweep the two stale comment references (scf_solver_isolation.sh, pw_scf.cuh)
- [x] Add Doxygen @brief per doc.style for scf::Phase, kind_of, Criterion, PhaseTraits<scf::Phase>
- [x] Run full check + test suite

## Testing strategy
- **Happy path**（test_scf_phase_plumbing）：断言 `slot_of(scf::Phase::fock)==0 … slot_of(check)==4`、`slot_count<scf::Phase>()==5`；`kind_of(scf::Phase::fock)==PhaseKind::SCF` 且 `phase_kind_name(PhaseKind::SCF)=="SCF"`；`make_binding(scf::Phase::fock, 0).kind==PhaseKind::SCF`；实例化 `CPUDriver<ATV::scf::Phase>`（必要时用最小本地子类）证明 `PhaseTraits` 特化就位；注册一个玩具 Criterion，断言 driver 每迭代拉取 `verdict()` 并按 Stop > Converged > Continue 合成（组件自身写不了循环状态）。
- **Edge / 编译期门**：枚举 `COUNT` 与 `slot_count` 一致；`Criterion` 为纯虚接口（不可实例化）。
- **既有门保持**：改动全部走 `slot_of`/`kind_of`，不含任何 `static_cast<int>` over Phase，故既有 `stage_typing_scrub` 保持绿；`kind_of` 是 `constexpr` 标签派发返回枚举值、非返回构造对象的自由函数，故 `no_factory_scrub` 保持绿。
- **Domain validation**：无。本环不触碰任何数值代码路径，SCF 能量在定义上不变；逐位一致的迁移验收与文献值回归属 03。
- **Regression example**：本环**不**含 `regressions/` 案例。01 不暴露任何端到端用户可见功能（SCF 驱动首次真正跑出能量在 03），故按整条链的切分，唯一的 `regressions/scf/<case>` 端到端案例及其 `type: runtime` 验收由 03 交付。01 的 runtime 关口是完整 ctest 套件（秒级）。

## Out of scope
- 五个阶段组件（含 `ScfConverged` Criterion）、`DualConvergenceCheck` 改造、状态成员与构造器引用接线 —— 均属 02/03a。
- 既有 `LoopHandler`（`ForceConverged`/`StopOnNaN`/`AllReplicasConverged`）+ `add_handler` API 的全树迁移与删除 —— 属 06-criterion-unify。
- 把 RHF/UHF/RKS/UKS 接到 `CPUDriver<scf::Phase>`；本环 `MeanFieldMethod` 仍拥有 `ScfLoop`/`SCFSolver`，方法能量路径不变 —— 属 03。
- 删除 `ScfLoop`/`SCFSolver`/`SCFConfig`/`make_config`/`scf_evaluator.h` —— 属 03/04。
- HF 收窄为纯 kind tag、原生 driver 嵌套、删 `sub_driver.h`、`SCFResultHandle` 收窄、外层消费方迁移 —— 属 04。
- `ATV_BUILD_RHF→ATV_BUILD_HF` 宏改名及 CMake/CI/slurm 同步，以及 CLAUDE.md/AGENTS.md 选项表与"scf::Phase stays banned"相位命名说明的文档同步 —— 属 05（本环仅改 `stage.h` 内那条被 scrub 豁免的代码注释）。
- CUDA `scf::Phase`：CUDA SCF 保持 inner-solver 形态，`scf::Phase` 本链为 CPU-only。
- 不新增任何 `.sh` 门；结构性断言由 evaluator 运行，不落地为提交脚本。
