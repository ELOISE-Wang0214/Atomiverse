---
title: SCF 驱动化重构 · 02 阶段组件化（五阶段组件，状态经构造器引用直连 + ScfConverged）
slug: scf-driver-refactor-02-stages
status: code-complete
created: 2026-07-20
chain: scf-driver-refactor
link: 02/05
---

# SCF 驱动化重构 · 02 阶段组件化

## Summary
本环把驱动 SCF 自洽循环的五个内部 helper（`HamiltonianBuilder`/`EigenSolver`/`DensityBuilder`/`DIISMixer`/`DualConvergenceCheck`）包装成五个绑定 `scf::Phase` 相位的 stage 组件（继承 `cpu::Component`，构造器直接吃离散参数**与状态引用**，**无任何 config 结构、无 Workspace 中间袋、无 ctx.scf 切片**——迭代态（F/P/C/S/X…）将是方法 driver 的私有成员，本环隔离测试中由测试局部矩阵扮演），同时把 `DualConvergenceCheck` 改造为一个 `ScfConverged` **Criterion 组件**（绑 `check` 相位：`launch()` 自算 ΔE/DIIS 残差存入自身真名成员，`verdict()` 判定容差——终止判据也是组件（决策⑨），无独立 Handler 概念，**不写 LoopMetrics**）。所有改动**纯增量**：现有 `ScfLoop`/`SCFSolver` 仍是唯一真正驱动四个方法的路径，新组件只在单测里被隔离验证（unit tier，秒级，无 e2e），四个方法的能量逐位不变。落地后，03 就能在一套已就绪、已验证 bit-identical 的 stage 组件之上，把 `MeanFieldMethod<D>` 接成 `CPUDriver<scf::Phase>`。

## Domain basis
本环无新物理；领域锚点是**行为保持**——每个 stage 组件对同一输入产出与其所包装 helper **逐位一致**的结果（单测直接对拍），因此把它们组合起来的 03 才可能复现四个方法逐位一致的能量。

一处关键领域事实决定了相位次序：本代码库的 DIIS 作用在 **Fock 矩阵**上、且**先于对角化**（Pulay 直接迭代子空间反演；`scf_solver.cpp:241→252→284` 建 F → 从未外推的 F 求 DIIS 误差并外推 → level shift → 对角化）。DIIS 参考：
- Pulay, *Chem. Phys. Lett.* **73**, 393 (1980)。
- Pulay, *J. Comp. Chem.* **3**, 556 (1982)。

因此 `mix`（DIIS）相位在数据依赖上**必须夹在 `fock` 与 `diagonalize` 之间**；（修订后的）01 已据此把 `scf::Phase` 定为 **`fock → mix → diagonalize → density → check`**，本环沿用，不再触碰 `stage.h`。

沿用整条链的决策日志（01 已固化 ①②③④⑧⑨）：② Solver = 封闭数值内核命名约定 —— `EigenSolver` 保持封闭内核、被 `diagonalize` 组件私有持有，**不**继承 `Component`；③（用户修订）状态通道 = **构造器直连引用** —— 迭代态是方法 driver 的私有成员（03a 定义），stage 组件构造时拿到各自读/写的具名引用；无 Workspace 中间袋、无 ctx.scf 切片。

## Design
五个 stage 组件全部继承 `cpu::Component`（它们是 SCF 循环的管线胶水，不是能量 KIND；能量 KIND 是方法本身，归 04 的 HF 纯 tag）。生命周期 `build()`（一次性分配）→ `launch()`（每迭代计算，不做固定缓冲分配）→ `destroy()`；每个组件 `default_bindings()` 用 `ATV::make_binding(scf::Phase::X, order)` 声明其相位，`launch()` 内用 `ctx.at(scf::Phase::X)` 自守，**经构造器接线的具名引用读写状态**（ctx 只用于相位自守，不承载 SCF 状态）。

**状态归属与接线契约（用户裁定：无 ScfWorkspace）** —— 相位间流动的基本物理量是**方法 driver 的私有成员**（03a 定义，真名：`S_`,`Hcore_`,`X_`,`F_alpha_`,`F_beta_`,`err_alpha_`,`err_beta_`,`C_alpha_`,`C_beta_`,`eig_alpha_`,`eig_beta_`,`P_alpha_`,`P_beta_`,`J_`,能量标量 `e_total_`/`e_one_`/`e_coulomb_`/`e_xc_grid_`/`e_exact_exchange_`，计数 `n_basis_`/`n_alpha_`/`n_beta_`/`n_dropped_`/`iter_`，模式 `restricted_`/`is_hf_`/`diis_active_`，`v_nuc_`）。**接线规则**：每个 stage 组件的构造器接收它**读**的 `const RealMatrix&`/`const double&` 与它**写**的 `RealMatrix&`/`double&`——xtensor 对象身份稳定，per-run resize 不失效引用。本环的隔离单测中，这些成员由测试局部矩阵扮演（组件不知道也不关心引用来自谁——这正是可测性来源）。
- **DIIS 历史留在 `mix` 组件内部**（`DIISMixer` 已自持 `std::deque` 历史）。理由：DIIS 子空间历史是该组件的私有加速器状态（派生量），非相位间共享的基本物理量。历史在每 run 起始由方法 driver 调 `mix` 组件的 `reset()` 清空；`launch()` 中 deque 逐迭代增长是加速器固有状态（文档化例外，非固定缓冲分配）。同理 `ScfConverged` 的 `prev_energy_` 留其内部，每 run reset。

**五个 stage 组件**：
1. `FockStage`（bind `fock`）—— 包装 `HamiltonianBuilder`。构造器吃离散方法参数（`restricted,is_hf,XCType,GridLevel,n_radial,n_angular,hybrid_alpha_override,const MolecularGrid* fixed_grid`），**不吃 SCFConfig**。为使新路径全程无 config，给 `HamiltonianBuilder` 增一个**离散参数入口**（prepare/step 存离散成员或加离散重载，与既有 SCFConfig 重载共享同一私有核心 → 逐位一致）；旧 `SCFSolver` 仍走其 SCFConfig 路径（SCFConfig 删除属 03）。`build()` 调 `prepare`；`launch()` 经构造器引用读 P、写 F/J/E_xc/E_exact_exchange/int_n。
2. `MixStage`（bind `mix`）—— 包装 `internal::diis_error()` + `DIISMixer` + `apply_level_shift`。构造器吃 `use_diis,diis_start_iter,diis_max_vectors,level_shift`。`launch()` 经构造器引用读 F/P/S/X：算 err_alpha/err_beta 写入其输出引用；满足 `iter>=diis_start_iter` 时 DIIS 外推 F；`|level_shift|` 超阈时施加 level shift。restricted 单通道、unrestricted 双通道分支与旧循环 `:245-282` 一致。
3. `DiagonalizeStage`（bind `diagonalize`）—— **私有持有** `internal::EigenSolver`（封闭内核，不继承 Component）。构造器无方法参数（路径选择 `n_dropped>0 ? step_orthogonal(F,X) : step(F,S)` 由引用的 `n_dropped` 决定）。`launch()` 经构造器引用读 F（+S 或 X）、写 C/eig；restricted 复用 α 结果。
4. `DensityStage`（bind `density`）—— 包装 `DensityBuilder` + 密度阻尼 + 能量求值。构造器吃 `restricted,is_hf,damping`。`launch()` 经构造器引用读 C/旧 P，建新 P（restricted factor 2）、施阻尼（DIIS 未激活时，"保留旧密度比例"约定，与 `:303-316` 一致）、按四分支（restricted×HF/KS、unrestricted×HF/KS）算 `e_total/e_one/e_coulomb` 写入其输出引用。注意 HF 能量用**旧 P**、KS 用**新 P**（与 `:318-348` 逐位一致）。
5. `ScfConverged`（bind `check`，**Criterion 类别**，01 落地的 `cpu/component/criterion.h`）—— `DualConvergenceCheck` 的整体化身，测量与判定一体：**测量**（`launch()` 经构造器引用读 `e_total` 与 err 矩阵，算 `dE = e_total - prev`（iter 0 为 +∞）与 `err_norm`（合并 α/β Frobenius），存入**自身真名成员** `delta_energy_`/`diis_error_`（供 trace 追加；**不写 LoopMetrics**——用户裁定"是什么就用什么命名"）+ **判定**（`verdict()` 只读自身成员：`|delta_energy_|<tol_energy_ && diis_error_<tol_grad_ ? Converged : Continue`，不看迭代计数——iter 0 因 delta_energy_=+∞ 自然 Continue）。构造器吃 `tol_energy,tol_grad`；`prev_energy` 留组件内部，每 run `reset()`。driver 每迭代**拉取** verdict 并按严重度合成——组件写不了循环状态，裁决通道不变。

**放置**：`ScfConverged` 位于 flag-gated 的 `cpu/scf/stages/scf_converged.h`——`Criterion` 类别基（`cpu/component/criterion.h`）无条件编译，而 SCF 具体判据受 `ATV_BUILD_RHF` 门控，故留在域目录。

**Reuse decision**（逐一了结 librarian 复用候选）：
- `reuse` 五个内部 helper —— `EigenSolver`/`DensityBuilder`/`DIISMixer`/`DualConvergenceCheck`/`diis_error()`/`apply_level_shift` 原样调用，不重实现任何数值；stage 组件只做"读输入引用 → 调 helper → 写输出引用"的胶水。
- `generalize` `HamiltonianBuilder` —— 增离散参数入口以脱离 SCFConfig（既有 SCFConfig 重载保留给旧循环），对应一条"Generalize `HamiltonianBuilder`"性质的任务，而非另造一个 Fock 构建器。
- `reuse` `cpu::Component` 生命周期与 `make_binding`/`ctx.at`/`default_bindings` 惯用法（component.h）—— stage 组件逐字仿既有组件形状，不新增基类。
- `reuse` `Criterion` 类别（01 落地）—— `ScfConverged` 是 Criterion 组件（测量+判定一体）；`ForceConverged`（loop_handler.h:91）仅作"只判 metrics"的结构参照。
- `new — 五个 stage 组件（含 ScfConverged）`：确无既有等价物；它们不重实现任何 `reuse`/`generalize` 候选，只组合之。无 Workspace 类型、无任何新自由函数。
- `EigenSolver` 明确**不**升格为组件（封闭数值内核，决策②），由 `DiagonalizeStage` 私有持有。

## Files to create or modify
- `include/atomiverse/cpu/scf/stages/fock_stage.h` (new) + `src/cpu/scf/stages/fock_stage.cpp` (new) — `FockStage`。
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h` + `src/cpu/scf/internal/hamiltonian_builder.cpp` — 增离散参数入口（与 SCFConfig 重载共享私有核心）。
- `include/atomiverse/cpu/scf/stages/mix_stage.h` (new) + `src/cpu/scf/stages/mix_stage.cpp` (new) — `MixStage`（DIIS + level shift，历史内置）。
- `include/atomiverse/cpu/scf/stages/diagonalize_stage.h` (new) — `DiagonalizeStage`（私有持 `EigenSolver`）。
- `include/atomiverse/cpu/scf/stages/density_stage.h` (new) + `src/cpu/scf/stages/density_stage.cpp` (new) — `DensityStage`（密度 + 阻尼 + 能量）。
- `include/atomiverse/cpu/scf/stages/scf_converged.h` (new) — `ScfConverged`（Criterion 类别：`launch()` 测量存自身真名成员，`verdict()` 判定容差）。
- `tests/cpu/scf/stages/test_scf_stages.cpp` (new) — 五组件的隔离失败优先单测（测试局部矩阵经引用接线）。
- `CMakeLists.txt` (root) — 在 `ATV_BUILD_RHF` 源列表（约 382–409 行）登记新增 `src/cpu/scf/stages/*.cpp`。
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_stages`。

## Tasks
- [x] Write failing tests for the SCF stage components + ScfConverged (tests/cpu/scf/stages/test_scf_stages.cpp), wiring test-local matrices into the stage constructors by reference, and register them in tests/CMakeLists.txt
- [x] Implement FockStage (bind fock) wrapping HamiltonianBuilder and generalize HamiltonianBuilder with a config-free discrete-parameter entry, in cpu/scf/stages/fock_stage.{h,cpp} + internal/hamiltonian_builder.{h,cpp}
- [x] Implement MixStage (bind mix) wrapping diis_error + DIISMixer + apply_level_shift with component-owned DIIS history, in cpu/scf/stages/mix_stage.{h,cpp}
- [x] Implement DiagonalizeStage (bind diagonalize) privately holding internal::EigenSolver (closed kernel, not a Component) in cpu/scf/stages/diagonalize_stage.h
- [x] Implement DensityStage (bind density) wrapping DensityBuilder + damping + the four-branch energy evaluation in cpu/scf/stages/density_stage.{h,cpp}
- [x] Implement ScfConverged (Criterion kind, bind check) measuring delta_energy_/diis_error_ into its own plainly-named members in launch() and judging tolerances in verdict(), in cpu/scf/stages/scf_converged.h
- [x] Register the new stage sources in the root CMakeLists.txt ATV_BUILD_RHF list
- [x] Run full check + test suite

## Testing strategy
- **Happy path（隔离对拍，unit tier）**：每个 stage 组件对同一输入产出与其所包装 helper **逐位一致**的输出——`FockStage.launch(ws)` 的 F/J/E_xc == `HamiltonianBuilder`；`MixStage` 的 err 与外推 F == `diis_error()`+`DIISMixer`+`apply_level_shift`；`DiagonalizeStage` 的 C/eig == `EigenSolver`；`DensityStage` 的 P/E == `DensityBuilder`+阻尼+能量式；`ScfConverged.launch()` 算出的 `delta_energy_`/`diis_error_` 数值正确，`verdict()` 在给定量下返回 `Converged`/`Continue`。测试用小尺寸合成矩阵（如 2×2 或 H₂/STO-3G 规模的单步输入），**绝不跑整轮 SCF 自洽**（那是 e2e，禁入 tests/）。
- **Edge cases**：restricted vs unrestricted、HF vs KS 四分支能量；DIIS 未激活（走阻尼）vs 激活；`iter==0`（`delta_e=+∞` → `ScfConverged` 必 Continue）；`n_dropped>0` 走 `step_orthogonal` 路径 vs 健康基走 dsygvd。
- **生命周期 / 契约**：`build()` 一次性分配、`launch()` 不做固定缓冲分配（DIIS deque 增长为文档化例外）；组件经 `ctx.at(scf::Phase::X)` 自守、经**构造器接线的具名引用**访问状态（ctx 不承载 SCF 状态）；`EigenSolver` 不继承 `Component`（编译期/结构断言）；新 stage 构造器**无 SCFConfig token、无内部 config 结构、无 Workspace 参数**（签名审查）。
- **Domain validation**：上述逐位对拍即本环的领域验证（每 stage 复现参考 helper）。整分子逐位一致能量属 03。
- **Additive 保持**：旧 `ScfLoop`/`SCFSolver` 仍驱动四方法，既有 `tests/cpu/scf/*`（含 `test_scf_solver`/`test_uhf_solver`/`test_uks_solver` 等）保持全绿、能量不变。
- **Regression example**：本环**不**含 `regressions/` 案例。02 不接任何端到端路径（唯一能跑出能量的仍是旧循环）；驱动组合后的 SCF 首次 e2e 跑文献值在 03，该 `regressions/scf/<case>` 及其 `type: runtime` 验收由 03 交付。本环 runtime 关口是完整 ctest 套件（秒级）。

## Out of scope
- 把五个 stage 组件接入 `CPUDriver<scf::Phase>`、`MeanFieldMethod<D>` 改造、per-method 文件拆分、run 级 setup（X/初始猜测）与 post-SCF（去移位轨道能级/自旋污染）编排 —— 属 03。
- 删除 `ScfLoop`/`SCFSolver`/`SCFConfig`/`make_config`、以及 `HamiltonianBuilder` 的 SCFConfig 重载 —— 属 03（本环仅**增**离散入口，不删 config）。
- `HF` 收窄为纯 tag、原生 driver 嵌套、删 `sub_driver.h`、`SCFResultHandle` 收窄、外层消费方迁移 —— 属 04。
- `ATV_BUILD_RHF→ATV_BUILD_HF` 改名及 CMake/CI/slurm/CLAUDE.md 同步 —— 属 05。
- CUDA：CUDA SCF 保持 inner-solver，`scf::Phase` 本链 CPU-only。
- ROHF：仅作 03/04 的可扩展性验收范例，不实现。
