---
title: Canonical orthogonalization with overlap-eigenvalue thresholding for SCF
status: code-complete
created: 2026-07-02
---

# Canonical orthogonalization with overlap-eigenvalue thresholding for SCF

## Summary

当基组接近线性相关时（如弥散基组、近重合原子中心），重叠矩阵 S 的最小本征值趋于零，现有 SCF 广义本征求解路径（LAPACK dsygvd，基于 Cholesky 分解）会数值失稳甚至直接失败。本工作在 `cpu/linalg` 中新增带本征值阈值筛除的正则正交化（canonical orthogonalization）：对 S 本征分解后丢弃本征值低于阈值 τ 的方向，得到 n×m 的矩形变换矩阵 X（m ≤ n），并让 SCF 内部本征求解器在检测到列被丢弃时切换到两步 XᵀFX 路径，健康基组则保持 dsygvd 快速路径不变（所有既有校验数值零漂移）。同时在 `docs/tutorial/ch02-rhf.md` 新增 Jensen 风格教程小节 "Why HF becomes a matrix equation / when the basis goes bad"（Roothaan–Hall、S^{−1/2}、线性相关），教程中的数值由新代码产出并纳入回归测试。

## Domain basis

（scientist 输出，PySCF 2.13.1 已验证；单位：能量 Hartree，长度 Å，S 无量纲）

- 本征分解 S = U s Uᵀ。Löwdin 对称正交化 X = U s^{−1/2} Uᵀ（n×n，当 s_i→0 时病态）；正则正交化 X = Ũ s̃^{−1/2}，仅保留 s_i ≥ τ 的列 → X 为 n×m 矩形矩阵，满足 XᵀSX = I_m。
- Roothaan–Hall：FC = SCε → F′ = XᵀFX（m×m），F′C′ = C′ε，C = XC′（n×m）。推论：n_mo = m < n_ao；C 为矩形；aufbau 占据与密度 P = 2 C_occ C_occᵀ（仍 n×n）、DIIS 误差 FPS − SPF（仍 n×n）在形式上不变；后续 MP2 依赖 n_vir = m − n_occ。
- PySCF 阈值约定：旧式 `remove_linear_dep_` 默认 1e-8（当 cond(S) > 1e10 时触发）；`canonical_orth_(S, thr=1e-7)` 先归一化 S（Snorm = D S D，D = diag(S_ii^{−1/2})）再对归一化本征值取阈值；PySCF ≥ 2.13 将其内建于默认本征求解器，对 S 本征值取阈值 1e-6（scf/hf.py:48-49, 1359-1360）。本工程采用：配置项暴露 τ，默认 1e-8，作用于**原始 S 本征值**；文档中注明 PySCF ≥ 2.13 使用 1e-6。
- 验证靶（PySCF 2.13.1，cart=True，本次生成）：H2/6-31G + 位于 (0, 0, 1e-5 Å) 的重复 6-31G 壳层中心 → n_ao = 6；S 最小本征值 3.53e-12、1.09e-10，其后 0.126；正则正交化恰好丢弃 2 列（X 形状 6×4）；带筛除的 RHF E = −1.1267347374 Ha，对照干净 H2/6-31G E = −1.1267339671 Ha（相差 8e-7）。H2 几何：H (0,0,0)、H (0,0,0.7414) Å。**注意**：该参考值是用零核电荷 ghost 中心生成的；本 spec 选定的 C++ 测试体系是 H2 + H⁺（真实质子，见 Design），核吸引势不同，参考能量必须按同一体系重新生成（见验收 ac-003）。丢弃列数 m = 4 仅取决于基函数重叠结构，与核电荷无关，可直接沿用。

## Design

> **Design correction (2026-07-02, user directive — "no factory functions;
> instantiate and use structs directly")**: 自由函数
> `canonical_orthogonalization(S, tau)` 取消。`CanonicalOrth` 改为带计算
> 构造函数的 struct：`CanonicalOrth(const RealMatrix& S, double tau)` 在构造
> 时完成本征分解与阈值筛除；成员 `X / n_dropped / s_min` 不变；全列被筛除时
> 构造函数抛 `std::runtime_error`。下文原文按此修正解读。

**测试体系选型（显式决定）**：工程的基组表按 Z 索引（`basis_data/`），不支持零核电荷 ghost 原子。为避免引入 ghost 机制，线性相关测试体系选定为 **H2 + H⁺**：三个真实 H 中心 —— H (0,0,0)、H (0,0,0.7414 Å)、H (0,0,1e-5 Å)（与第一个中心近重合），`build_integrals(frame, BasisSet::Pople631G, /*charge=*/+1)` 使电子数保持 2。该体系 n_ao = 6，基组近简并结构与 ghost 版本相同（预期仍恰好丢 2 列，m = 4），但 PySCF 参考能量需对**完全相同**的几何/电荷/cart=True 重新生成并钉入测试。

**新公开符号**：

- `ATV::cpu::linalg::CanonicalOrth`（新 struct）：`{ RealMatrix X /* n×m 列主序 */, int n_dropped, double s_min }`。
- `ATV::cpu::linalg::canonical_orthogonalization(const RealMatrix& S, double tau) -> CanonicalOrth`：`xt::linalg::eigh(S)` 本征分解（升序），保留 s_i ≥ τ 的列，X = Ũ diag(s̃^{−1/2})。τ ≤ 0 时不筛除（等价保留全部列）。遵循 "no fallback values"：若筛除后 m = 0 抛 `std::runtime_error`。所有中间量与返回值均 `double`（AccReal 恒 double）。
- `SCFConfig::lindep_threshold`（新字段，`double`，默认 `1e-8`）：S 本征值筛除阈值；`0` 表示禁用（保持 dsygvd 直解路径，行为与现状逐位一致）。Doxygen 注明单位（无量纲的 S 本征值）与 PySCF ≥ 2.13 默认 1e-6 的差异。

**求解器接线（最小行为差异原则）**：

- `SCFSolver::run()`（`src/cpu/scf/scf_solver.cpp:91`）现调用 `build_orthogonalization_matrix`（Löwdin，无阈值）产出 `X_` 仅供 DIIS 误差用；改为调用 `canonical_orthogonalization(S_cm, config_.lindep_threshold)`，`X_` 成为可能矩形的 n×m 矩阵（`scf_solver.h:99` 的成员注释同步更新），并记录 `n_dropped`。
- `internal::EigenSolver` 增加两步路径入口 `step_orthogonal(F, X)`（Xᵀ F X → eigh → C = X C′，复用既有 `linalg::diagonalize_fock`），保留 `step(F, S)`（dsygvd）。**路径选择**：仅当 `n_dropped > 0` 时 SCF 主循环（scf_solver.cpp:199-205）与初猜对角化（scf_solver.cpp:99 的 Hcore 解）切换到两步路径；`n_dropped == 0` 时一律走 dsygvd —— 健康基组的全部既有校验数值不得移动。
- DIIS 误差 `internal::diis_error(F, P, S, X_)` = Xᵀ(FPS−SPF)X：X 为 n×m 时结果为 m×m，维度自洽，无需改 DIIS 代码，但需测试覆盖。
- `SCFResult` 结构不变；文档注明筛除发生时 `eigenvalues_*` 长度为 m ≤ n_ao、`C_*` 为 n×m（n_mo ≤ n_ao），`P_*` 仍为 n×n。

**生命周期/所有权**：`CanonicalOrth` 按值返回；`X_` 仍由 `SCFSolver` 持有并每次 `run()` 重建（几何相关）。`eigensolve.h` 为 `cpu::scf` 与 `cpu::gpw` 共享，仅新增符号、不改既有符号签名，GPW 不受影响。

**教程**：在 `docs/tutorial/ch02-rhf.md` 的 RHF 部分之后新增小节 "Why HF becomes a matrix equation — and when the basis goes bad"（英文），Jensen 风格：Roothaan–Hall 推导（为何广义本征问题）、S^{−1/2} 与正则正交化、线性相关的物理来源与 τ 的作用；给出 H2 + H⁺ 演示的 S 最小本征值、丢弃列数、能量对比数字。教程引用的全部数字由新测试 `test_scf_lindep.cpp` 钉死（同一体系、同一代码路径），形成教程回归。

## Files to create or modify

- `include/atomiverse/cpu/linalg/eigensolve.h` — 新增 `CanonicalOrth` + `canonical_orthogonalization()`
- `include/atomiverse/cpu/scf/scf_solve.h` — `SCFConfig::lindep_threshold`；`SCFResult` 的 n_mo ≤ n_ao 文档
- `include/atomiverse/cpu/scf/internal/eigen_solver.h` — `EigenSolver::step_orthogonal(F, X)` 两步路径
- `include/atomiverse/cpu/scf/scf_solver.h` — `X_` 成员注释更新 + `n_dropped` 记录
- `src/cpu/scf/scf_solver.cpp` — 正交化构建切换、路径选择、初猜对角化
- `tests/cpu/linalg/test_canonical_orthogonalization.cpp` (new) — linalg 单元测试（新目录）
- `tests/cpu/scf/test_scf_lindep.cpp` (new) — 干净/近简并 RHF + 教程数字回归
- `tests/CMakeLists.txt` — 注册两个新测试（默认门内，非 E2E 段）
- `docs/tutorial/ch02-rhf.md` — 新教程小节

## Tasks

- [x] Write failing tests for canonical_orthogonalization (tests/cpu/linalg/test_canonical_orthogonalization.cpp, new dir; register in tests/CMakeLists.txt): synthetic near-singular S → X shape n×m, ‖XᵀSX − I_m‖ < 1e-12, exact n_dropped/s_min; tau=0 keeps all n columns; all-dropped case throws
- [x] Implement CanonicalOrth + canonical_orthogonalization in include/atomiverse/cpu/linalg/eigensolve.h with Doxygen (@brief/@param/@return, units noted dimensionless)
- [x] Regenerate PySCF 2.13.1 reference (cart=True) for the chosen H2 + H⁺/6-31G system (H at 0,0,0; 0,0,0.7414; 0,0,1e-5 Å; charge=+1) and record E_ref, s_min spectrum, expected m=4 as constants in tests/cpu/scf/test_scf_lindep.cpp
- [x] Write failing tests for SCF linear-dependence handling (tests/cpu/scf/test_scf_lindep.cpp; register in tests/CMakeLists.txt): (i) clean H2/6-31G RHF energy unchanged vs current value to 1e-12 with default lindep_threshold; (ii) H2 + H⁺ near-duplicate converges with exactly 2 dropped columns and |E − E_ref| < 1e-6 Ha; (iii) lindep_threshold=0 on the clean system reproduces the current dsygvd result bitwise; (iv) tutorial-quoted numbers pinned
- [x] Implement SCFConfig::lindep_threshold (default 1e-8, 0 disables) in include/atomiverse/cpu/scf/scf_solve.h and document n_mo ≤ n_ao on SCFResult fields
- [x] Implement EigenSolver::step_orthogonal in include/atomiverse/cpu/scf/internal/eigen_solver.h and wire src/cpu/scf/scf_solver.cpp + include/atomiverse/cpu/scf/scf_solver.h: build X_ via canonical_orthogonalization, select two-step path only when n_dropped > 0 (SCF loop and initial Hcore diagonalization), keep dsygvd otherwise; verify DIIS error dimensions with rectangular X_
- [x] Write tutorial section "Why HF becomes a matrix equation — and when the basis goes bad" in docs/tutorial/ch02-rhf.md (English, Jensen-style: Roothaan–Hall, S^{−1/2} vs canonical orthogonalization, thresholds incl. PySCF ≥ 2.13 = 1e-6), quoting only numbers produced by test_scf_lindep.cpp
- [x] Verify no-regression: full existing cpu/scf test set passes untouched and both new tests complete within the default ctest gate (< 1 s each)
- [x] Run full check + test suite

## Testing strategy

- **Happy path（单元）**：合成对称正定 S 注入一个 1e-12 量级本征方向 → `canonical_orthogonalization` 返回正确形状、XᵀSX = I_m 至 1e-12、`n_dropped`/`s_min` 精确匹配。
- **边界**：τ = 0（不筛除，X 为 n×n 且与 Löwdin 张同一空间）；τ 大到丢弃全部列 → 抛异常（no-fallback 约定）；τ 恰好等于某本征值（≥ 为保留）。
- **无回归**：干净 H2/6-31G 与既有全部 SCF 测试（H2O/STO-3G 等）在默认 `lindep_threshold=1e-8` 下能量零漂移（1e-12 级）——健康基组 n_dropped=0 → dsygvd 路径未被触碰；`lindep_threshold=0` 与现状逐位一致。
- **领域验证（science.required）**：H2 + H⁺/6-31G 近简并体系：dsygvd 直解路径在此体系失败或失稳（作为动机记录），筛除路径收敛、恰丢 2 列（m=4）、|E − E_PySCF| < 1e-6 Ha（PySCF 2.13.1，cart=True，同一几何与电荷重新生成）。
- **教程回归**：ch02 小节引用的每个数字（s_min、n_dropped、两个能量）与 `test_scf_lindep.cpp` 钉死的常量一致。
- **门控**：两个新测试均为毫秒级小体系，进默认 ctest 门（约 16 s 预算），不进 E2E 禁用列表。

## Out of scope

- CUDA SCF（`cuda::scf::SCFSolver`）、GPW、平面波路径的线性相关防护（`eigensolve.h` 新符号可被后续复用，但本 spec 不接线）。
- 零核电荷 ghost 原子基组支持（已用 H2 + H⁺ 方案规避；ghost 支持留待独立 spec）。
- PySCF 式 S 归一化阈值（Snorm = D S D）与 cond(S) 自动触发器 —— 本 spec 采用固定 τ 作用于原始 S 本征值，触发器作为后续增强。
- MP2 及其它 n_vir = m − n_occ 消费方（hf-robust 链后续环节）。
- SCF 梯度（`scf_gradient`）在筛除激活时的正确性验证（矩形 C 对梯度装配的影响单独立 spec 评估）。
