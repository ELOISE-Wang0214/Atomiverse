---
title: RHF/UHF virtual-orbital level shift with PySCF-parity DIIS ordering
status: code-complete
created: 2026-07-02
---

# RHF/UHF virtual-orbital level shift with PySCF-parity DIIS ordering

## Summary

为 `cpu::scf::SCFSolver` 增加虚轨道 level shift（能级移动）收敛稳定器：`SCFConfig` 新增 `level_shift`（单位 Hartree，默认 0.0，即完全不改变现有行为），在 SCF 循环中按 PySCF 顺序（damping → DIIS 外推 → level shift → 对角化）对 Fock 矩阵施加 F′ = F + λ(S − S·D_occ·S)，收敛后额外做一次 λ=0 的对角化以保证输出的轨道能不被移位污染。同时把现有 damping/DIIS 次序在代码注释与 Doxygen 中正式化（现有实现是对角化后的**密度**阻尼，与 PySCF 的 Fock 阻尼不同——保留并记录差异，不改动）。用户可见成果还包括 `docs/tutorial/ch02-rhf.md` 新增一节 "SCF convergence is a numerical problem"（教程内容为英文），讲解 DIIS、damping、level shift 与 Jensen 式的收敛观，教程中所有数值由新代码生成并被回归测试锁定。

## Domain basis

- **Level shift（Saunders & Hillier）**：F′ = F + λ(S − S·D_occ·S)，其中 D_occ = C_occ C_occᵀ 关于 S 幂等（RHF 时 D_occ = P/2；PySCF 在 `scf/hf.py:1133` 传 `dm*0.5`）。谱效应：虚轨道本征值整体上移 λ（Hartree），占据轨道不变——拉开 HOMO–LUMO 间隙、抑制占据/虚空间的振荡。参考文献：V. R. Saunders & I. H. Hillier, *Int. J. Quantum Chem.* **7**, 699 (1973), DOI [10.1002/qua.560070407](https://doi.org/10.1002/qua.560070407)。
- **PySCF 2.13.1 `get_fock` 次序（`hf.py:1128-1134`）**：damping → DIIS → level shift。DIIS 外推的是**未移位**的 Fock；误差向量 FPS−SPF 也由未移位 Fock 构造；shift 在外推**之后**施加。|λ| ≤ 1e-4 时跳过。PySCF 的 damping `f = f(1−α) + f_prev·α` 仅在 cycle < `diis_start_cycle`−1 时生效。默认值：`level_shift=0`，`damp=0`，`diis_start_cycle=1`。
- **收敛后处理（`hf.py:214-221`）**：收敛后以 λ=0 额外对角化一次，使报告的轨道能不含 +λ 污染——**必须实现**，否则虚轨道 ε 全部偏移 λ。
- **策略**：默认 λ=0；迭代过程中 λ 恒定（PySCF 不做衰减）；衰减调度属非 parity 扩展 → 见 Out of scope。UHF：PySCF 支持分自旋 λ；本实现保持标量 λ 同时作用于两个自旋通道（记录为限制）。
- **验证目标（本轮 PySCF 生成）**：H2O/6-31G RHF，`level_shift=0.5` 对比 `0`：收敛能量一致到 1.1e-13 Ha；max|Δε| = 2.5e-7。测试判据：|E(λ=0.5) − E(λ=0)| < 1e-9 Ha，且最终 ε 相对 λ=0 无移位（≤ 1e-6 Ha）。几何（Å）：O (0, 0, 0.1173)；H (0, 0.7572, −0.4692)；H (0, −0.7572, −0.4692)。PySCF 锚点（cart=True, conv_tol=1e-12）：E(H2O/6-31G RHF) = −75.9839744727 Ha。

## Design

**触及实体**：`SCFConfig`（`scf_solve.h`）、`SCFSolver::run` 的迭代循环（`scf_solver.cpp:170-314`）与循环后收尾段（`:316-335`）、`internal::` 求解器辅助层。

**新符号**：

- `SCFConfig::level_shift`（`double`，默认 `0.0`，Hartree）。带 Doxygen（含单位），并在注释中声明 PySCF parity 语义与 |λ| ≤ 1e-4 跳过规则。
- `ATV::cpu::scf::internal::apply_level_shift(F, S, D_occ, lambda)` — 纯函数辅助器，返回（或原位改写）F + λ(S − S·D_occ·S)。放入新的 internal 头/源对（与 `internal::diis_error` 自由函数先例一致）：不继承 `ComponentBase`，不注册进 Plan，facade 公有头不包含它（测试可以包含）。所有累加用 `double`（AccReal 约定）。

**循环内接线**（`SCFSolver::run`，现有步骤编号）：在第 2 步 DIIS 外推**之后**、第 3 步对角化**之前**插入：若 `|config_.level_shift| > 1e-4`，对 F_alpha（restricted 时 D_occ = 0.5·P_alpha，因为 restricted 密度以占据数 2 构建）与 F_beta（unrestricted 时 D_occ 为各自旋的 P_sigma，本身逐自旋幂等）施加移位，S 取 `S_cm`。DIIS 误差矩阵（`:179-182`）保持在移位前构造——现有代码顺序天然满足，只需插入点正确并加注释固化。

**收敛后 λ=0 对角化**：循环退出后（`:314` 之后、spin-contamination 之前），若 `result_.converged` 且 |λ| > 1e-4：用收敛密度重建一次 Fock（`hamiltonian_->step`，不经 DIIS、不移位），对角化，覆写 `result_.eigenvalues_alpha/beta` 与 `result_.C_alpha/C_beta`（restricted 时 beta 镜像 alpha）。能量不重算（能量本就由未移位 Fock/密度求得，且已验证与 λ=0 一致到 1e-13）。

**次序正式化（文档而非改码）**：`scf_solver.cpp` 循环加一段结构注释，写明流水线 build F → DIIS error（unshifted）→ DIIS extrapolate → level shift → diagonalize → density → density damping（仅 pre-DIIS），并显式记录：本实现的 damping 是对角化后的**密度**线性混合（`:217-230`，"fraction of old" 约定），不同于 PySCF 的 Fock damping；两者均只在 DIIS 未激活时生效。**不**改为 Fock damping。

**生命周期/所有权**：level shift 完全在 `SCFSolver::run` 栈内生效，无跨迭代状态、无 Frame 触碰、无新 workspace；`SCFResult` 结构不变（`trace` 中迭代期的 ε 允许含移位，最终报告值保证无移位）。

**教程**：`docs/tutorial/ch02-rhf.md`（SCF 循环章）末尾新增英文小节 "SCF convergence is a numerical problem"：(1) 为什么裸迭代会振荡（占据/虚子空间翻转的图像）；(2) damping（含本引擎密度阻尼约定）；(3) DIIS（引用 `internal::DIISMixer` 的误差向量 Xᵀ(FPS−SPF)X）；(4) level shift（Saunders–Hillier 方程 + 谱效应 + "收敛后去移位"规则）；(5) Jensen 式的实践建议（先 shift 稳住、后 DIIS 收尾；病态体系的组合拳）。硬案例：**拉伸 H2O**（两条 O–H 拉至 ~1.7 Å，6-31G RHF，`init_guess=Hcore`、`use_diis=false`、`damping_factor=0`）——λ=0 在 max_iter 内振荡不收敛，λ=1.0 收敛且能量与 DIIS 参考一致；教程给出两组迭代表，全部数字由实现期实际运行产出并被 `test_cpu_scf_level_shift` 锁定（实现期若发现该设定不够判别性，可在测试与教程中同步微调拉伸倍数/λ，但"λ>0 可见地救活或显著加速收敛"这一判据不变）。

## Files to create or modify

- `include/atomiverse/cpu/scf/scf_solve.h` — `SCFConfig` 增加 `level_shift` 字段 + Doxygen
- `include/atomiverse/cpu/scf/internal/level_shift.h` (new) — `internal::apply_level_shift` 声明
- `src/cpu/scf/internal/level_shift.cpp` (new) — 实现
- `src/cpu/scf/scf_solver.cpp` — 循环内接线、收敛后 λ=0 对角化、次序结构注释
- `CMakeLists.txt` — 在 `src/cpu/scf/internal/` 源列表（约 :303）加入 `level_shift.cpp`
- `tests/cpu/scf/test_accelerators.cpp` — helper 级单元测试（谱效应、跳过阈值）
- `tests/cpu/scf/test_level_shift.cpp` (new) — 求解器级 parity / 去移位 / 硬案例 + 教程回归
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_level_shift`（命名约定 `test_{backend}_{domain}_{component}`）
- `docs/tutorial/ch02-rhf.md` — 新增 "SCF convergence is a numerical problem" 小节

## Tasks

- [x] Write failing unit tests for `internal::apply_level_shift` in `tests/cpu/scf/test_accelerators.cpp`: synthetic F/S/D_occ 下虚本征值精确上移 λ、占据本征值不变（generalized eigensolve 验证）；|λ| ≤ 1e-4 时逐元素恒等
- [x] Write failing solver-level tests in `tests/cpu/scf/test_level_shift.cpp` (new) and register `test_cpu_scf_level_shift` in `tests/CMakeLists.txt`: (a) H2O/6-31G RHF 锚点 −75.9839744727 Ha（λ=0，容差 1e-6 Ha）；(b) |E(λ=0.5) − E(λ=0)| < 1e-9 Ha；(c) 最终 max|Δε| ≤ 1e-6 Ha（去移位验证）；(d) 拉伸 H2O 硬案例：λ=0 不收敛（converged==false）、λ=1.0 收敛且 |E − E_ref(DIIS)| < 1e-8 Ha 且迭代数被锁定
- [x] Implement `level_shift` field with Doxygen (units: Hartree, default 0.0, skip rule |λ| ≤ 1e-4) in `include/atomiverse/cpu/scf/scf_solve.h`
- [x] Implement `internal::apply_level_shift` in `include/atomiverse/cpu/scf/internal/level_shift.h` (new) + `src/cpu/scf/internal/level_shift.cpp` (new); add the source to root `CMakeLists.txt`
- [x] Wire level shift into `SCFSolver::run` in `src/cpu/scf/scf_solver.cpp`: after DIIS extrapolation, before diagonalization; restricted 用 0.5·P_alpha，unrestricted 标量 λ 施加于两自旋通道；DIIS 误差保持由未移位 Fock 构造
- [x] Implement post-convergence λ=0 extra diagonalization in `src/cpu/scf/scf_solver.cpp`（重建未移位 Fock → 对角化 → 覆写 `result_.eigenvalues_*` / `result_.C_*`）
- [x] Add pipeline-ordering structural comment in `src/cpu/scf/scf_solver.cpp` documenting damping → DIIS → level shift → diagonalize and the density-damping-vs-PySCF-Fock-damping difference
- [x] Write tutorial section "SCF convergence is a numerical problem" in `docs/tutorial/ch02-rhf.md` (English), covering DIIS / damping / level shift / Jensen-style practice, with all cited numbers generated by `test_cpu_scf_level_shift` runs
- [x] Verify default-ctest wall time impact: `test_cpu_scf_level_shift` 若 >1 s 则移入 `tests/CMakeLists.txt` 尾部的 `ATV_RUN_E2E_TESTS` 禁用列表（helper 单元测试留在默认门内）
- [x] Run full check + test suite

## Testing strategy

- **Happy path**：λ=0 默认路径——现有全部 SCF 测试（`test_cpu_scf_solver`、`test_cpu_scf_accelerators`、UHF/UKS 套件）不改一字通过，参考能量逐位不变；H2O/6-31G RHF 锚点 −75.9839744727 Ha（1e-6 Ha）。
- **Edge cases**：|λ| ≤ 1e-4 跳过（helper 恒等 + solver 路径不额外做收敛后对角化）；restricted 与 unrestricted 两分支（标量 λ 双通道）；`use_diis=false` 下 shift 单独工作；DIIS 激活后 shift 仍施加于外推结果之后。
- **Domain validation（科学）**：(1) 谱定理级单元测试：虚 ε 上移恰为 λ、占据 ε 不变；(2) H2O/6-31G λ=0.5 vs 0 能量一致 <1e-9 Ha、最终 ε 去移位 ≤1e-6 Ha（PySCF 2.13.1 交叉验证目标 1.1e-13 / 2.5e-7，测试判据放宽为规格值）；(3) 拉伸 H2O 硬案例可复现"level shift 救活发散"，收敛能量与 DIIS 参考一致 <1e-8 Ha。
- **教程回归**：ch02 新节引用的每个数字（迭代数、能量、Δε）均由 `test_cpu_scf_level_shift` 断言锁定；教程与测试用同一几何/配置常量。
- **ctest 门**：默认 `ctest` 保持 ~16 s 量级；新增默认门测试合计 <2 s，否则移入 E2E 列表。

## Out of scope

- **Fock damping**：现有实现是密度阻尼（`scf_solver.cpp:217-230`），保留并文档化差异；改为 PySCF 式 Fock damping 是行为变更，另立 spec。
- **分自旋 level shift**（PySCF 的 `(λα, λβ)` 元组）——本实现标量 λ 双通道，作为限制记录在 Doxygen。
- **λ 衰减调度**（动态 shift）——PySCF 不做，属非 parity 扩展。
- **CUDA `cuda::scf::SCFSolver` 的 level shift**——CPU 先行；CUDA 镜像另立 spec。
- **KS/UKS 专项验证**——接线对 KS 分支同样生效（代码路径共享），但科学验证仅覆盖 HF；KS 硬案例另议。
- **GPW/PW 求解器**的收敛稳定器。
