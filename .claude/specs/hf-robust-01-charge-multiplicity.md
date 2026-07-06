---
title: Gaussian SCF charge/multiplicity resolver and closed-shell validity guard
status: code-complete
created: 2026-07-02
---

# Gaussian SCF charge/multiplicity resolver and closed-shell validity guard

## Summary

当前 `include/atomiverse/cpu/scf/scf_solve.h:255` 的文档声称 `default_config`（Restricted）路径"Throws-via-solver if n_electrons is odd"，但代码中不存在任何此类抛出：奇数电子体系经 `rhf_config` → `n_occ() = n_electrons / 2` 被静默截断，得到 `n_alpha + n_beta = N_e − 1` 的物理错误配置并照常求解。本工作为 Gaussian SCF 配置层新增一个由 (电荷, 自旋多重度) 到 (n_alpha, n_beta) 的解析函数，并在 `SCFSolver::run` 入口加一道配置有效性守卫：电子数与占据数不一致、Restricted 下开壳层、以及非法 (N_e, M) 奇偶/符号组合，统一抛出一条指明用户可见参数的 `std::invalid_argument`（符合仓库 "no fallback values" 规则），使原本虚假的文档变为真实行为并同步修正。配套在教程 ch02（闭壳层有效性）与 ch03（电荷/多重度约定与解析器）中新增 Jensen 式（概念 → 数学 → 算法 → 局限）章节，其中的算例数字必须由新代码产生并由回归测试逐条锁定。不重构 `SCFConfig`，diff 保持最小。

## Domain basis

（scientist 输出，已用 PySCF 2.13.1 验证）

- 电子数与占据数解析：

  N_e = Σ_A Z_A − q（q 为分子电荷）；多重度 M = 2S + 1，故 n_α − n_β = M − 1；

  n_α = (N_e + M − 1) / 2，n_β = (N_e − M + 1) / 2。

- 有效性条件：
  - N_e 与 M − 1 同奇偶（偶 N_e ⇒ 奇 M；奇 N_e ⇒ 偶 M），否则 n_α/n_β 非整数；
  - n_β ≥ 0 ⇒ M ≤ N_e + 1；
  - RHF（闭壳层）要求 M = 1 且 N_e 为偶数。
- PySCF 约定对照：`mol.spin = 2S = n_α − n_β`（不是 M）；奇偶违例抛 `RuntimeError "Electron number 2 and spin 1 are not consistent"`（gto/mole.py:2386），n_β < 0 抛 `AssertionError`。Atomiverse 对两类违例统一抛**一条**描述性 `std::invalid_argument`，消息以用户可见的 (charge, multiplicity)/(N_e, M) 为键。测试断言条件本身，不比对 PySCF 消息文本。
- 验证目标（(charge, M) → (n_α, n_β)）：H2O(0,1)→(5,5)；OH⁻(−1,1)→(5,5)；H2O⁺(+1,2)→(5,4)；O 原子(0,3)→(5,3)。必须拒绝：H2O(0,2)、H2O(+1,1)、H(0,4)。
- 参考文献：F. Jensen, *Introduction to Computational Chemistry*, 3rd ed., Wiley (2017), ISBN 978-1-118-82599-0（ch. 3–4，闭壳层/开壳层与多重度约定）；A. Szabo, N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996), ISBN 978-0-486-69186-2；PySCF: Q. Sun *et al.*, WIREs Comput. Mol. Sci. **8**, e1340 (2018), DOI 10.1002/wcms.1340；Q. Sun *et al.*, J. Chem. Phys. **153**, 024109 (2020), DOI 10.1063/5.0006074。
- 单位：电子数/占据数为无量纲整数；SCF 能量内部单位 Hartree（本工作不触碰能量学）。

## Design

> **Design correction (2026-07-02, user directive — "no factory functions;
> instantiate and use structs directly")**: `resolve_occupations()` 自由函数与
> `uhf_config`/`uks_config` 的 multiplicity 重载全部取消。取而代之：
> `SpinResolution` 带校验构造函数（`SpinResolution(int n_electrons, int
> multiplicity)`，非法组合抛 `std::invalid_argument`），用户直接实例化后把
> `occ.n_alpha / occ.n_beta` 喂给**既有**的原始计数 config 构建器。以下原文
> 中涉及 `resolve_occupations` 与重载的段落按此修正解读。

全部改动落在 `cpu/scf` 配置/校验层与文档，不新建文件层级、不改 `SCFConfig` 结构体、不动求解循环。

新增公有符号（header-only，`include/atomiverse/cpu/scf/scf_solve.h`，`ATV::cpu::scf` 命名空间）：

- `struct SpinResolution { int n_alpha; int n_beta; };` — 简单聚合，无生命周期问题。
- `SpinResolution resolve_occupations(int n_electrons, int multiplicity);` — 纯函数（`inline`），实现 Domain basis 公式；奇偶违例或 `n_beta < 0` 或 `multiplicity < 1` 或 `n_electrons < 0` 时抛出一条 `std::invalid_argument`，消息中给出 N_e 与 M 的实际值及违反的条件。电荷在 `build_integrals`/`make_atom_list`（src/cpu/scf/basis_set.cpp:72,85）阶段已折入 `IntegralPack::n_electrons`，故解析器以 (N_e, M) 为输入即可覆盖 (charge, M) 语义。
- 重载 config 构建器（保留现有原始计数重载，遵循仓库既有"自由 config 构建器"风格，不引入工厂）：
  - `uhf_config(const IntegralPack& ints, int multiplicity, double tol = 1e-8)`；
  - `uks_config(const IntegralPack& ints, int multiplicity, XCType xc = XCType::LSDA, double tol = 1e-8)`。

  两者内部调用 `resolve_occupations(ints.n_electrons, multiplicity)` 填充 `n_alpha/n_beta`，其余字段与对应原始计数重载逐字段一致。重载消解经核实是确定的：`uks_config` 因 `XCType` 是 enum class（int 不隐式转换）无歧义；`uhf_config(ints, a, b)` 两个 int 实参精确匹配原始计数重载（int→double 是转换，精确匹配优先），行为可预期，测试中锁定。

修改的既有行为：

- `SCFSolver::run`（src/cpu/scf/scf_solver.cpp:64，唯一 choke point——`scf_solve()` 与 `rhf_solve`/`uhf_solve`/`uks_solve` 皆经此）在进入迭代前校验：
  1. `config_.n_alpha ≥ 0 && config_.n_beta ≥ 0`；
  2. `config_.n_alpha + config_.n_beta == integrals.n_electrons`；
  3. `spin == Restricted` ⇒ `n_alpha == n_beta`。

  任一违反抛 `std::invalid_argument`，消息包含 n_electrons、n_alpha、n_beta、自旋类型及"请检查 charge/multiplicity"指引。奇电子体系走 `rhf_config`（n_α=n_β=N/2，和为 N_e−1）恰好触发条件 2，`default_config` 文档的声明由此成真。
- 修正 `scf_solve.h` 文档：:255 行改述为真实抛出行为（`std::invalid_argument` from `SCFSolver::run`），:270 行 `@param ints ...（assumed even）` 同步改写；`rhf_config`（:180）补注奇电子体系将在求解入口被拒绝。新公有符号按仓库规则配 Doxygen `@brief`/`@param`（含单位/无量纲说明）/`@return`/`@throws`。

教程（英文正文，匹配现有 book 风格，Jensen 式 概念→数学→算法→局限）：

- `docs/tutorial/ch02-rhf.md` 新增一节"Closed-shell validity"：为什么 RHF 只对偶 N_e、M=1 成立；`n_occ = N_e/2` 的截断陷阱；新守卫的抛出示例（H2O⁺ 走 `rhf_config` 抛出）。
- `docs/tutorial/ch03-open-shell.md` 新增一节"Charge and multiplicity conventions"：N_e/M/n_α/n_β 公式；与 PySCF `mol.spin = 2S` 约定的对照表；`resolve_occupations` 与 `uhf_config(ints, multiplicity)` 用法；算例表（H2O、OH⁻、H2O⁺、O 原子）与三个被拒组合。教程中每一组数字均与回归测试断言一一对应（两侧以注释互相引用）。

不做：CUDA 后端镜像（`cuda::scf` 配置路径独立，见 Out of scope）、`SCFConfig` 字段增删、`build_integrals` 签名变更。

## Files to create or modify

- include/atomiverse/cpu/scf/scf_solve.h — `SpinResolution` + `resolve_occupations` + multiplicity 重载；修正 :255/:270/:180 文档
- src/cpu/scf/scf_solver.cpp — `SCFSolver::run` 入口配置校验守卫
- tests/cpu/scf/test_spin_resolution.cpp (new) — 解析器表 + 拒绝用例 + 守卫用例 + 教程回归表
- tests/CMakeLists.txt — 注册 `test_cpu_scf_spin_resolution`（`atv_add_cpp_test`，默认 gate，ms 级）
- docs/tutorial/ch02-rhf.md — "Closed-shell validity" 节
- docs/tutorial/ch03-open-shell.md — "Charge and multiplicity conventions" 节

## Tasks

- [x] Write failing tests for the SpinResolution constructor acceptance table and the three rejection cases (tests/cpu/scf/test_spin_resolution.cpp, new; register test_cpu_scf_spin_resolution in tests/CMakeLists.txt)
- [x] Write failing tests for the SCFSolver::run validity guard — odd-electron RHF (H2O⁺ via rhf_config/rhf_solve) throws; Restricted with n_alpha≠n_beta throws; occupation-sum mismatch throws; neutral-H2O RHF energy unchanged (tests/cpu/scf/test_spin_resolution.cpp)
- [x] Implement SpinResolution with a validating constructor (throws std::invalid_argument) and Doxygen @brief/@param/@throws in include/atomiverse/cpu/scf/scf_solve.h
- [x] Verify SpinResolution composes with the existing raw-count uhf_config/uks_config builders (no new overloads — user design directive; FeedsUhfConfigRawCounts / FeedsUksConfigRawCounts tests)
- [x] Implement the configuration validity guard at the top of SCFSolver::run in src/cpu/scf/scf_solver.cpp, throwing one descriptive std::invalid_argument
- [x] Fix the false docstrings in include/atomiverse/cpu/scf/scf_solve.h (default_config :255 "Throws-via-solver", :270 "assumed even", rhf_config :180) to describe the now-real behavior
- [x] Write the tutorial sections — "Closed-shell validity" in docs/tutorial/ch02-rhf.md and "Charge and multiplicity conventions" in docs/tutorial/ch03-open-shell.md, Jensen-style, worked numbers taken from the validation table
- [x] Verify tutorial-regression linkage: every worked (charge, M, n_alpha, n_beta) tuple and every rejected combination in the tutorial sections appears as an assertion in tests/cpu/scf/test_spin_resolution.cpp, with cross-referencing comments on both sides
- [x] Run full check + test suite (default ctest gate stays within ~16 s budget)

## Testing strategy

- **Happy path**：`resolve_occupations` 对四个验证目标逐一断言 (n_α, n_β)：H2O(0,1)→(5,5)、OH⁻(−1,1)→(5,5)、H2O⁺(+1,2)→(5,4)、O(0,3)→(5,3)；`uhf_config(ints_H2O⁺, /*M=*/2)` 产出 `n_alpha=5, n_beta=4` 且其余字段与 `uhf_config(ints, 5, 4)` 逐字段相等；中性 H2O 的 `rhf_solve` 能量与改动前逐位一致（守卫零开销于合法路径）。
- **Edge cases**：拒绝 H2O(0,2)（奇偶违例）、H2O(+1,1)（奇 N_e 配 M=1）、H(0,4)（n_β<0），均断言抛 `std::invalid_argument`（断言条件与异常类型，不比对消息全文）；`multiplicity<1`、`n_electrons<0` 拒绝；守卫三条件各自的独立触发；`uhf_config(ints, a, b)` 双 int 实参仍解析到原始计数重载（重载消解锁定）。
- **Domain validation**：验证表与 PySCF 2.13.1 行为对齐（同接受/同拒绝，见 Domain basis）；仅测条件不测消息文本。
- **教程回归**：测试文件中以注释标注 "tutorial ch02/ch03 worked example" 的断言块与教程数字一一对应；教程改数字必须同步改测试，反之亦然。
- **门控**：新测试为纯配置/解析逻辑加一次 STO-3G H2O 小体系求解，ms–亚秒级，进默认 ctest gate（~16 s 预算不破）；不新增 E2E 表尾条目。

## Out of scope

- CUDA 后端（`cuda::scf::SCFSolver` / `cuda::scf::SCFEvaluator`）的镜像守卫——考虑过随本 spec 一并做以满足"同一特性同一目录跨后端正交"，但 CUDA 配置路径（energy/density only，无 UHF 多重度入口）差异较大，留作后续 spec。
- `SCFConfig` 结构体重构（如把 charge/multiplicity 收进 config）——设计指引明确要求最小 diff。
- `build_integrals` / `make_atom_list` 签名或行为变更（电荷折算已正确）。
- `IntegralPack::n_occ()` 本身的语义变更——保留为"闭壳层双占据轨道数"，其误用由求解器守卫拦截。
- `uhf_solve` / `uks_solve` 自由函数新增 multiplicity 形参重载（config 构建器层已覆盖该需求）。
- PySCF 消息文本兼容。
