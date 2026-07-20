---
title: SCF 驱动化重构 · 03a 方法驱动化（MeanFieldMethod→CPUDriver<scf::Phase>，逐位对拍，参考实现留存）
slug: scf-driver-refactor-03a-convert
status: approved
chain: scf-driver-refactor
link: 03a/06
created: 2026-07-20
---

# SCF 驱动化重构 · 03a 方法驱动化

## Summary
本环把 `MeanFieldMethod<D>` CRTP（`RHF`/`UHF`/`RKS`/`UKS` 的共享参数家 + 唯一 `run(frame)` 入口）改造成同时派生 `CPUDriver<scf::Phase>`，用 02 的五个 stage 组件（Fock/Mix/Diagonalize/Density/ConvergenceMeasure）驱动自洽循环：构造器把方法成员**直接**接进各 stage 组件构造器（全程无 config bag）、经 `add_criterion` 注册 `ScfConverged` 默认收敛判据（Criterion 组件，测量+判定一体）、以 `max_iter` 为预算。所有对外入口（两个 `run` 重载、`forces`、组件 `build`/`launch`+ASPC、`mp2`/`fci`）统一汇入一个**私有成员** `drive_scf_(frame, guess, opts)` = 前置装配 → `CPUDriver<scf::Phase>` 循环 → 后置收尾。方法按文件拆成 `cpu/scf/{rhf,uhf,rks,uks}.{h,cpp}`（各文件只留构造器默认 + `resolve_occupation` + `name()`），`methods.h` 转为伞头。**关键：`ScfLoop`/`SCFSolver`/`SCFConfig`/`make_config` 全部留在树里**——它们不再驱动方法，只作**逐位对拍的活参考**（仍跑 `test_scf_solver*`），`03b` 在四方法能量逐位一致被证明后才删除。四方法能量对旧 `SCFSolver` 参考 ≤1e-12 Ha，既有 `tests/cpu/scf/*` 全绿。

## Domain basis
本环无新物理；领域锚点是**行为保持**——迁移前后 `RHF`/`UHF`/`RKS`/`UKS` 的收敛能量逐位一致（≤1e-12 Ha），因为新驱动路径调用的是 02 已逐位对拍过的同一批 stage 组件、其内核又是同一批数值 helper。

**"参考实现先于替换物存活"原则（记入本环设计根据）** —— 逐位对拍的唯一硬关口要求有一个可信参考。若在单一 03 里同时 (i) 把方法改成驱动、(ii) 删除 `SCFSolver`、(iii) 迁移求解器单测，一旦对拍失手，你已经删掉了用来对比的黄金实现。故本链把"转换"（03a）与"退役"（03b）分开：**03a 保留 `SCFSolver`/`ScfLoop`/`SCFConfig` 作为在树的活参考（其单测继续绿）、并据此证明新路径逐位一致；03b 才删除这套已被证冗余的旧机器。** 这是"新路径证毕之前不动黄金"的纪律，不是任务数腾挪。

**`drive_scf_` 编排（单一私有成员漏斗；`drive_scf_` 是私有成员函数，不是自由函数，全环不新增任何 `make_*` 命名）** —— 四个对外入口都是薄壳，只设置 `DriveOpts{write_energy, forces, aspc, unwrap}`，共用同一个私有 `drive_scf_`：
- **前置装配（per-run，循环之前，方法自持，非伪组件）**：解包快照（仅可变路径，经私有 RAII 守卫类——用户裁定无自由 helper 函数）→ 重建私有积分 → 正交化矩阵 X（Löwdin vs 规范正交 + `n_dropped`）→ 初始猜测（MINAO/Hcore/外部）→ 装填方法的成员状态矩阵（`S_`/`Hcore_`/`X_`/`P_*` 等，stage 组件已在构造时接好引用）→ `MixStage`/`ScfConverged` 复位。
- **循环（唯一的驱动部分）**：`CPUDriver<scf::Phase>` 每迭代走 `fock → mix → diagonalize → density → check`，poll `ScfConverged`，预算 `max_iter`。
- **后置收尾（per-run，循环之后，方法自持）**：去移位轨道能级（收敛且启用移位时）→ UHF 自旋污染 → 填 `SCFResult` → 能量 eV 写 `frame.meta` → `on_evaluate_` 回调 → 恢复解包。

沿用整链决策日志（01/02 已固化 ①②③④⑤⑧）。为何前置/后置留在方法而非做成伪组件：积分重建 + 初始猜测 + X 依赖具体方法的 basis/charge/occupation 知识、每 run 只跑一次、且不是可组合的每迭代科学阶段——按"单次一处的代码不抽函数""per-run 装配是方法的职责"，它们是方法 `run()` 的前后言，不是 stage。详见 Design 的前置/后置表。

## Design
`MeanFieldMethod<D>` 本环同时派生 `CPUDriver<scf::Phase>`（循环节点）**与** `HF`（仍是 Component，本环过渡）。构造器仿 `vib::Hessian` 预设惯用法：按固定有序链 `fock→mix→diagonalize→density→check` 用 `on(scf::Phase::X).add(Stage{成员...})` 组装五个 02 stage 组件、`add_criterion(ScfConverged{tol, conv_tol_grad>0 ? conv_tol_grad : std::sqrt(tol)})`（Criterion 组件，01 增设的注册通道）、`CPUDriver::build()` 编译 Plan。方法成员（`tol`/`use_diis`/`diis_max_vectors`/`diis_start_iter`/`level_shift`/`damping_factor`/`grid_level`/`n_radial`/`n_angular`/`hybrid_alpha_override`/`fixed_grid`/`spin_`/`method_`/`xc_type_`/…）**直接**流入各 stage 构造器——**stage 接线全程无 `SCFConfig`、内部亦无**。

**驱动路径删除**（方法不再用旧循环）：`sync_loop_config_`/`ensure_loop_`/`run_on_pack_` 及 `loop_`/`corrector_loop_`（`ScfLoop`）从方法的**驱动路径**移除，由 `CPUDriver<scf::Phase>` + `drive_scf_` 取代。`drive_scf_` 经 librarian 点名的 `CPUDriver` 设施（`setup`/`step`/`teardown`、`poll_handlers_`、`metrics`、`compile`、`run_nested`）驱动 Plan，不复制循环。

**前置装配（`drive_scf_` 循环前，均为私有成员步骤）：**

| 步骤 | 现位置 | 归属 + 根据 |
|---|---|---|
| 位置解包 + 快照（原自由函数 `unwrap_for_scf`/`restore_wrapped` 改为**方法私有 RAII 守卫类**，ctor 解包、dtor 恢复——用户裁定无自由 helper 函数） | `methods.cpp:35,62` | 仅可变路径（`run(Frame&)`/`launch`）；`run(const Frame&)` 保持几何原样——保行为 |
| 重建积分（`IntegralPack{frame,basis,charge,compute_eri}` + `resolve_occupation`） | `methods.cpp:153` | 复用现有 `integrals_sp_`（单一积分库，同时供 handle）；积分是方法私有工作区（Frame 纯度），依赖成员，每 run 重建 |
| X（Löwdin vs 规范正交 + `n_dropped`） | `scf_solver.cpp:122-139` | 前置 → `ws.X`/`ws.n_dropped`；每 run 一次特征分解，几何相关，非可组合每迭代阶段 |
| 初始猜测（MINAO / Hcore→对角→密度 / 外部） | `scf_solver.cpp:154-213` | 前置 → `ws.P_alpha/P_beta`；方法据 `init_guess` 成员决定，复用对角化+密度内核一次 |
| 装填成员状态矩阵（`n_basis_`/计数/`restricted_`/`is_hf_`/`v_nuc_`/`S_`/`Hcore_`/`X_`/`P_*`） | 散落 | 前置——状态是方法私有成员，stage 组件构造时已接好引用，装填即就位 |
| `MixStage.reset()` + `ScfConverged.reset()` | `scf_solver.cpp:151-152` | 前置——每 run 清 DIIS 历史 + `prev_energy` |

**后置收尾（`drive_scf_` 循环后，均为私有成员步骤）：**

| 步骤 | 现位置 | 归属 + 根据 |
|---|---|---|
| 去移位轨道能级（重建未移位 Fock，对角化） | `scf_solver.cpp:402-421` | 收敛且移位有效时，复用 Fock+对角化内核一次——收敛后一次性清理，非循环阶段 |
| UHF 自旋污染 | `scf_solver.cpp:424-427` | 仅非限制；收敛 C 的派生性质 |
| 填 `SCFResult`（能量分解、eig/C/P、trace、iterations、converged） | `scf_solver.cpp:364-386` | 后置读成员状态终态；每迭代 trace 由 `ScfConverged.launch()` 追加进方法 `SCFResult.trace` |
| 能量 eV → `frame.meta[kEnergy]`（`write_energy`/组件路径） | `methods.cpp:185,323,325` | 后置（可变路径）——Frame 是外层货币 |
| `on_evaluate_` 回调 | `methods.cpp:328` | 后置（组件路径）——用户钩子 |
| ASPC 预测/校正/记录 | `methods.cpp:259-310` | 包裹前置（预测→喂猜测）+ 后置（校正+记录），组件路径；与初始猜测、收敛密度联动 |
| `restore_wrapped` | `methods.cpp:333` | 后置（可变路径）——撤销解包 |

**`mp2()`/`fci()` 仍为成员**：`compute_eri=true` 下 `drive_scf_` 收敛并在 `integrals_sp_` 物化稠密 ERI，二者随后照旧调 `mp2_energy(*integrals_sp_, last_result())` / `fci_energy(...)`；仅改管线（改调 `drive_scf_`，不再 `run_on_pack_`）。

**分阶段 run/launch 债务（显式，04 结清）**：本环 `MeanFieldMethod<D>` 同时是 `CPUDriver<scf::Phase>`（`run`）与 `HF`/`Component`（`build`/`launch`），`run` 与 `launch` **共存**——刻意、临时地违反"永不共存"规则，以让把它当 HF 组件绑定并调 `launch` 的外层 `MD`/`Evaluate` 消费方继续编译；两入口都汇入同一 `drive_scf_`。代码处留注释指向 04（04 删 `build`/`launch`、`HF` 收窄为纯 tag、只留 `run`）。

**留存的过渡债务（不进新驱动路径）**：`SCFSolver`/`ScfLoop`（参考 + `test_scf_solver*`）、`SCFConfig`/`make_config`/`publish_config_`（仅填旧 `SCFResultHandle.config` 快照给 observer）原样保留；`HamiltonianBuilder` 保留其 SCFConfig 重载（旧参考用）+ 02 的离散入口（新路径用）。`SCFResultHandle` 本环**不动**（其 `config` 字段删除 + 消费方迁移由 03b 随 `SCFConfig` 删除一并处理）。

**Reuse decision**（逐一了结复用候选）：
- `reuse` `CPUDriver` 设施（`setup`/`step`/`teardown`/`poll_handlers_`/`metrics`/`compile`/`run_nested`）—— `drive_scf_` 据此驱动 Plan，不复制循环。
- `reuse` 02 的五个 stage 组件（含 `ScfConverged`）—— 方法构造器把成员状态引用 + 离散参数接进各 stage；`drive_scf_` 装填成员并复位有状态 stage。
- `reuse` `vib::Hessian` 预设惯用法（构造器固定有序链 + `validate_components`）—— `MeanFieldMethod` 构造器仿之。
- `reuse` 现有方法积木：`IntegralPack`、`CanonicalOrth`/Löwdin（linalg）、`minao`/`scf_guess`、`scf_gradient`、`compute_spin_contamination`、`mp2_energy`/`fci_energy`、`ASPCPredictor` —— 前置/后置原样调用，不重实现；`unwrap_for_scf`/`restore_wrapped` 两个自由函数**吸收进私有 RAII 守卫类**（无自由 helper 函数）。
- `reuse` `MeanFieldMethod<D>` CRTP 本体（成员、CRTP 链、`with_gradient`/`on_evaluate`/`predictor` 链式）—— 保留，仅换内部驱动。
- `new — drive_scf_ + 前置/后置私有成员` + 每方法 `rhf/uhf/rks/uks` 文件。**不新增任何自由函数、任何 `make_*`。**
- 明确**不复用进新路径**（仅留作参考/handle 债务）：`SCFSolver`/`ScfLoop`/`SCFConfig`/`make_config`。

## Files to create or modify
- `include/atomiverse/cpu/scf/methods.h` — 转为伞头：包含 `mean_field_method.h` + `rhf/uhf/rks/uks.h`；保留 `using SCF = RKS;` 别名。
- `include/atomiverse/cpu/scf/mean_field_method.h` (new) — 从 methods.h 迁入并改造 `MeanFieldMethod<D>`（派生 `CPUDriver<scf::Phase>` + `HF`）、`SCFResultHandle`；声明 `drive_scf_` + 前置/后置私有成员、五个 stage 成员；保留 `build`/`launch`（债务）。
- `src/cpu/scf/mean_field_method.cpp` (new) — `MeanFieldMethod<D>` 模板实现：构造器 stage 接线、`drive_scf_`、前置、后置、两 `run` 重载、`forces`、`mp2`、`fci`、`build`、`launch`、显式实例化。
- `include/atomiverse/cpu/scf/rhf.h` (new) + `src/cpu/scf/rhf.cpp` (new) — `RHF`（构造器默认 + `resolve_occupation` + `name()`）。
- `include/atomiverse/cpu/scf/uhf.h` (new) + `src/cpu/scf/uhf.cpp` (new) — `UHF`。
- `include/atomiverse/cpu/scf/rks.h` (new) + `src/cpu/scf/rks.cpp` (new) — `RKS`。
- `include/atomiverse/cpu/scf/uks.h` (new) + `src/cpu/scf/uks.cpp` (new) — `UKS`。
- `src/cpu/scf/methods.cpp` (delete) — 内容拆入 `mean_field_method.cpp` + 四个 per-method `.cpp`。
- `CMakeLists.txt` (root) — 在 `ATV_BUILD_RHF` 源列表（约 382–409）去掉 `methods.cpp`，加入 `mean_field_method.cpp` + `rhf.cpp` + `uhf.cpp` + `rks.cpp` + `uks.cpp`。
- `tests/cpu/scf/test_method_driver_parity.cpp` (new) — 四方法逐位对拍 + stage 接线单测。
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_method_driver_parity`。

## Tasks
- [ ] Write failing per-method parity + stage-wiring unit tests (tests/cpu/scf/test_method_driver_parity.cpp) and register them in tests/CMakeLists.txt
- [ ] Convert MeanFieldMethod<D> (in methods.h/methods.cpp) to derive CPUDriver<scf::Phase>, wire the five 02 stages config-free in the ctor (member→stage-ctor), register the ScfConverged default criterion via add_criterion, budget=max_iter, and remove the ScfLoop drive path (SCFSolver/SCFConfig/make_config retained only as parity reference + legacy handle snapshot) with Doxygen briefs
- [ ] Implement the private member drive_scf_ prologue (position-unwrap via a private RAII guard class on mutating paths — no free helper functions, rebuild integrals_sp_, X Löwdin-vs-canonical + n_dropped, initial guess, member-state fill, MixStage/ScfConverged reset)
- [ ] Implement the drive_scf_ epilogue (de-shift orbital energies, UHF spin analysis, SCFResult population + per-iteration trace, energy eV write, on_evaluate_, restore_wrapped)
- [ ] Route both run overloads, forces(), and the Component build()/launch()+ASPC path through drive_scf_ (staged run/launch coexistence debt; in-code pointer to 04)
- [ ] Adjust mp2()/fci() to consume integrals_sp_ + converged last_result() after drive_scf_(compute_eri=true), staying members
- [ ] Split the methods into cpu/scf/{rhf,uhf,rks,uks}.{h,cpp} (ctor defaults + resolve_occupation + name() only), extract MeanFieldMethod<D> into cpu/scf/mean_field_method.{h,cpp}, make methods.h an umbrella, delete methods.cpp, and update the root CMakeLists ATV_BUILD_RHF source list
- [ ] Run full check + test suite

## Testing strategy
- **Parity（逐位对拍，unit tier，秒级）**：`test_method_driver_parity` 对 `RHF`/`UHF`/`RKS`/`UKS` 各构造新驱动方法与在树 `SCFSolver` 参考（同一 `IntegralPack`/`SCFConfig`），在既有 `tests/cpu/scf/*` 已用的**同尺度小分子**（H₂/H₂O STO-3G 级）上跑到收敛，断言 `last_result().energy` 与参考差 ≤1e-12 Ha，并抽查 eig/P 一致。**为何仍属 unit tier 而非 e2e**：与已在快套件里的 `test_scf_solver`/`test_uhf_solver`/`test_uks_solver` 同尺度、秒级完成；无 out-of-library 消费方、无文献值、无计时、无第三方包读取。out-of-library 文献值 e2e 是 03b 的 `regressions/scf/<case>`。
- **Stage 接线**：断言 `MeanFieldMethod` 派生 `CPUDriver<scf::Phase>`，构造器把五个 stage 按 `fock→mix→diagonalize→density→check` 接入 Plan、经 `add_criterion` 注册 `ScfConverged`；`drive_scf_` 路径无 `SCFConfig`/`SCFSolver`/`ScfLoop` token（签名/grep 审查）；`drive_scf_` 为私有成员；全环无新 `make_*`。
- **两入口一致**：同一分子经 `run(frame)` 与经组件 `launch`（绑进 `Evaluate`）得到一致能量——证明二者共用 `drive_scf_`。
- **Edge**：`restricted`/`unrestricted`×HF/KS；DIIS 未激活（阻尼）vs 激活；`level_shift≠0` 的去移位后处理；ASPC 预测命中/回退；`mp2`/`fci` 结果不变。
- **Domain validation**：四方法逐位对拍即本环领域验证。整分子 out-of-library 文献值验证属 03b。
- **Additive/参考保持**：`SCFSolver`/`SCFConfig` 及 `test_scf_solver*` 全绿不变（活参考）。
- **Regression example**：本环**不**含 `regressions/` 案例——驱动化后的 e2e 文献值案例连同其 `type: runtime` 验收由 03b 交付。本环 runtime 关口是完整 ctest 快套件（秒级）。

## Out of scope
- 删除 `ScfLoop`/`SCFSolver`/`SCFConfig`/`make_config`/`HamiltonianBuilder` 的 SCFConfig 重载、退役 `test_scf_solver*`/`scf_solver_isolation.sh`、迁移 `delta_scf`/`koopmans`/`level_shift`/`accelerators` 测试、`SCFResultHandle.config` 字段删除 + 消费方迁移、`regressions/scf/<case>` e2e —— 属 03b。
- 删 `build`/`launch`（结束 run/launch 共存）、`HF` 收窄为纯 tag、原生 driver 嵌套、删 `sub_driver.h`、`SCFResultHandle` 剩余 `integrals` 字段收窄 + 外层消费方迁移 —— 属 04。
- `ATV_BUILD_RHF→ATV_BUILD_HF` 改名 —— 属 05。
- CUDA；ROHF 实现（仅作 ac 设计评审范例）。
