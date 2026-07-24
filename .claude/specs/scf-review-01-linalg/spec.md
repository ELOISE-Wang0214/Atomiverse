---
slug: scf-review-01-linalg
title: "SCF review 01 — hoist shared linear algebra (sym_power, nonsymmetric eig, PCG, dRPA solve)"
status: code-complete
created: 2026-07-24
revised: 2026-07-24
domain: scf / linalg
queue: scf-review
depends_on: []
grilled: true
---

# SCF review 01 — hoist shared linear algebra

## Summary
把 `cpu/scf` 里重复实现的四类稠密线性代数原语上移到欠填充的 `cpu/linalg` 基础设施层（以及 `OrbitalHessian`），使 rpa / g0w0 / bse / tdhf / adc / eom(ee,ea,ip) / cphf / ucphf 各处调用同一份代码。这是纯行为保持型重构：用户可见行为不变，由现有单元测试与 regression 门控，dRPA 相关能量保持 ≤1e-12 Ha 不变量。落地后 `cpu/linalg` 从 matrix/eigensolve/density 三件套补齐对称矩阵幂、非对称特征值、预条件共轭梯度三项通用原语，纠正当前把通用线代散落在 SCF 方法文件里的分层漂移。

## Domain basis
本子规范唯一带物理不变量的项是 dRPA（direct RPA 相关能）。其余三项（`sym_power`、`eig_nonsymmetric`、`solve_pcg`）是初等稠密线性代数原语，无需文献引用，其物理语境分别为 Casida/RPA 的 (A±B)^{1/2} 变换、EOM-CCSD/ADC/TDHF 的非厄米本征问题、CPHF/CP-KS 线性响应求解。

- **对称矩阵幂**（`sym_power`）：对实对称 M = U λ Uᵀ，M^p = U λ^p Uᵀ（λ 无量纲或 Hartree，视矩阵而定）。可选 `floor` 把近零本征值下截，避免 p<0 时的放大/发散——这正是 tdhf 的 4 参数变体所做的。
- **dRPA / Casida**（rpa.cpp:155 现有公式，保持不变）：
  以 A、B 为轨道旋转 Hessian 分块（Hartree），
  M = (A−B)^{1/2} (A+B) (A−B)^{1/2}，
  激发能 ω_n = sqrt(eig(M))（ω_n² 是 M 的本征值），
  相关能 E_c = ½ (Σ_n ω_n − Tr A)。
  其中 (A−B)^{±1/2} 即由 `sym_power` 给出，`kOmega2Tol = 1e-10` 守卫 ω²≥0。
  文献：Furche, *Phys. Rev. B* **64**, 195120 (2001), DOI 10.1103/PhysRevB.64.195120；Casida, in *Recent Advances in DFT* (World Scientific, 1995), 第 155 页公式一致。
  **不变量（scientist 项 4，解决其开放问题 (c)）：** 本次上移保持相同 LAPACK 例程、相同排序、相同单线程归约顺序 ⇒ 目标为逐位一致（bitwise）；验收以可复核的 **≤1e-12 Ha（ω_n 与 E_c）** 为诚实门槛（bitwise 是其子集），并在共享版本中原样保留 `kOmega2Tol=1e-10` 守卫。
- **单位约定**：A/B/M 与 E_c 内部以 Hartree 计（与 SCF 内部一致），regression 直接以 Ha 断言（与既有 `regressions/scf/direct-rpa-h2-sto3g` 相同口径）；公共 frame 能量的 eV 口径不在本子规范触及范围。

## Design
**触及实体**：`cpu/linalg`（infra 层，当前仅 matrix/eigensolve/density）新增两项自由函数原语和一个新文件；`cpu/scf/internal/OrbitalHessian` 增补一个方法；九个 SCF 方法 `.cpp` 收敛到共享符号。四项互不耦合，可独立测试。

- **A. `sym_power`（对称矩阵幂）**：在 `cpu/linalg/eigensolve.h` 内**紧邻** `build_orthogonalization_matrix` 增加 header-inline 自由函数 `RealMatrix sym_power(const RealMatrix& M, double p, double floor = 0.0, std::string_view label = "sym_power")`——`label` 进入错误消息（audit 核实：rpa/g0w0 的 3 参数版与 tdhf 的 4 参数版都带 `what`/上下文串，折入为尾参保持诊断信息）。header-inline 的无归属纯数学是工厂禁令下被明确认可的例外（librarian item 4）。`build_orthogonalization_matrix(S)` 与 `CanonicalOrth::lowdin` 改为调用 `sym_power(S, -0.5)`（对健康基集逐位一致）；tdhf 的 4 参数 floor 变体折进默认 `floor` 参数。
- **B. `eig_nonsymmetric`（非对称特征值 / dgeev 封装）**：当前树内**不存在**任何共享 dgeev 封装（librarian item 1）。在 `cpu/linalg/eigensolve.h` 新增 `linalg::eig_nonsymmetric`，严格沿用同文件 `diagonalize_fock_generalized`/`dsygvd_` 已确立的"显式工作区尺寸 + info≠0 抛出、不做 LWORK=-1 查询"模式（该模式的"捆绑 LAPACK 的 LWORK=-1 损坏"注释就在此文件），返回实/虚部特征值（供 EOM/ADC/TDHF 的实/共轭对使用）。命名、构造、错误处理与既有 `dsygvd_` 封装同构，读起来与现有代码一致。
- **C. `solve_pcg`（预条件共轭梯度）**：当前树内**不存在**任何 PCG（librarian item 6：minimize/cg 是非线性 PR+ 优化器，soscf/planewave 的 Davidson 是本征求解器，均不适配）。新建文件 `cpu/linalg/linear_solve.h`，提供模板化 `solve_pcg(matvec, diag, b, tol, max_iter)`——方法无关的求解器归属比 `scf/internal` 有更深复用面；cphf/ucphf 两个真实调用点已过"第二调用点"门槛。未在 `max_iter` 内收敛时抛出，保持 cphf 现有的非收敛抛出契约。
- **D. dRPA 求解上移到 `OrbitalHessian`**：`internal/orbital_hessian.h` 的 `OrbitalHessian` 已拥有 `dense_A(Channel, alpha_hf)`/`dense_B`（文档写明"可被 CPHF/TDHF 迭代求解器复用"），是 dRPA 稠密求解的自然归属。为其增加 `direct_rpa` 方法，封装 (A−B)^{1/2}(A+B)(A−B)^{1/2} 的构造、对称本征分解、E_c 求和与 X/Y 幅度组装，返回 `{omega, e_c, X, Y}`。**Audit 核实（g0w0.cpp:69-113）**：g0w0 的 `solve_drpa` 不止要 E_c——它由 T=(A−B)^{1/2}、Ti=(A−B)^{-1/2} 与本征向量 Z 组装完整的 X、Y 幅度矩阵供屏蔽库仑 W 使用；故共享方法的结果结构体必须导出 X/Y（rpa 只读 omega+e_c，忽略 X/Y；T/Ti 是方法内部量，不导出）。两条调用路径共用同一方法，不重复实现。

### Reuse decision
- **`sym_power` → reuse** `cpu/linalg/eigensolve.h`：作为 header-inline 自由数学**扩展**该文件，`build_orthogonalization_matrix` + `CanonicalOrth::lowdin` 成为其调用者；折入 tdhf floor 变体。
- **dRPA solve → generalize** `OrbitalHessian`（`internal/orbital_hessian.h`）：新增 `direct_rpa` 方法，同时服务 `rpa.cpp` 与 `g0w0.cpp`（对应"Generalize …"任务，不新建平行实现）。
- **`eig_nonsymmetric` → new**（无既有共享封装）——严格遵循同文件 `dsygvd_`/`diagonalize_fock_generalized` 的显式工作区 + 抛出模式。
- **`solve_pcg` → new**（无既有 PCG）——新建 `cpu/linalg/linear_solve.h`，按 eigensolve.h 的显式工作区/抛出风格书写。
- **`CCSDSolution`（librarian item 1 owner）→ 不触及**：item 1 的 make_tau/EomImds 部分属于 scf-review-03；本子规范只取其 dgeev 封装部分。

## Files to create or modify
- include/atomiverse/cpu/linalg/eigensolve.h — 新增 `sym_power`、`eig_nonsymmetric`；令 `build_orthogonalization_matrix` + `CanonicalOrth::lowdin` 走 `sym_power`
- include/atomiverse/cpu/linalg/linear_solve.h (new) — `solve_pcg`
- include/atomiverse/cpu/scf/internal/orbital_hessian.h — 声明 `OrbitalHessian::direct_rpa`
- src/cpu/scf/orbital_hessian.cpp — 实现 `direct_rpa`
- src/cpu/scf/rpa.cpp — 改用 `linalg::sym_power` + `OrbitalHessian::direct_rpa`，删本地副本
- src/cpu/scf/g0w0.cpp — 改用 `linalg::sym_power` + `OrbitalHessian::direct_rpa`，删本地副本
- src/cpu/scf/bse.cpp — 改用 `linalg::sym_power`，删本地副本
- src/cpu/scf/tdhf.cpp — 改用 `linalg::sym_power`（floor 变体）+ `linalg::eig_nonsymmetric`，删本地副本
- src/cpu/scf/adc.cpp — 改用 `linalg::eig_nonsymmetric`，删 dgeev 封装
- src/cpu/scf/eom_ee_ccsd.cpp — 改用 `linalg::eig_nonsymmetric`，删 dgeev 封装
- src/cpu/scf/eom_ea_ccsd.cpp — 改用 `linalg::eig_nonsymmetric`，删 dgeev 封装
- src/cpu/scf/eom_ip_ccsd.cpp — 改用 `linalg::eig_nonsymmetric`，删 dgeev 封装
- src/cpu/scf/cphf_pack.cpp — 改用 `linalg::solve_pcg`，删本地 PCG
- src/cpu/scf/ucphf_pack.cpp — 改用 `linalg::solve_pcg`，删本地 PCG
- tests/cpu/linalg/test_eigensolve.cpp (new) — `sym_power` + `eig_nonsymmetric`
- tests/cpu/linalg/test_linear_solve.cpp (new) — `solve_pcg`
- tests/CMakeLists.txt — 注册 `test_cpu_linalg_eigensolve`、`test_cpu_linalg_linear_solve`
- regressions/scf/scf-review-01-linalg/reg_scf-review-01-linalg.cpp (new) — dRPA E_c + CPHF 静态极化率，硬编码文献值
- regressions/CMakeLists.txt — `atv_add_regression(scf scf-review-01-linalg)`

## Tasks
- [x] Write failing unit tests for linalg::sym_power and linalg::eig_nonsymmetric (tests/cpu/linalg/test_eigensolve.cpp → TestEigensolve; register in tests/CMakeLists.txt)
- [x] Implement sym_power (defaulted floor) and eig_nonsymmetric in include/atomiverse/cpu/linalg/eigensolve.h with Doxygen @brief/@param(units)/@return; route build_orthogonalization_matrix + CanonicalOrth::lowdin through sym_power
- [x] Write failing unit test for linalg::solve_pcg (tests/cpu/linalg/test_linear_solve.cpp → TestLinearSolve; register in tests/CMakeLists.txt)
- [x] Implement solve_pcg in include/atomiverse/cpu/linalg/linear_solve.h (new) with Doxygen and throw-on-non-convergence
- [x] Generalize OrbitalHessian in include/atomiverse/cpu/scf/internal/orbital_hessian.h (+ src/cpu/scf/orbital_hessian.cpp) with a direct_rpa method serving rpa.cpp and g0w0.cpp (result exposes omega, e_c, and the X/Y amplitude matrices g0w0's screened Coulomb consumes)
- [x] Port sym_power and dRPA callers (rpa.cpp, g0w0.cpp, bse.cpp, tdhf.cpp) to linalg::sym_power + OrbitalHessian::direct_rpa; delete the in-file copies
- [x] Port dgeev callers (tdhf.cpp, adc.cpp, eom_ee_ccsd.cpp, eom_ea_ccsd.cpp, eom_ip_ccsd.cpp) to linalg::eig_nonsymmetric; delete the five wrappers
- [x] Port PCG callers (cphf_pack.cpp, ucphf_pack.cpp) to linalg::solve_pcg; delete the two copies
- [x] Add regression example regressions/scf/scf-review-01-linalg/reg_scf-review-01-linalg.cpp (public API: dRPA E_c + CPHF static polarizability; hard-coded cited goldens) and register in regressions/CMakeLists.txt
- [x] Run full check + test suite
- [x] Hygiene cleanup (simplify ran clean: dead includes rpa/g0w0, stale RED comments; manual handoffs: PD-guard dedup → refactor backlog, dgeev EigPair glue → scf-review-03)

## Testing strategy
单元测试仅在 `tests/` 下，路径镜像 `cpu/linalg/`，每个测试针对单一函数，硬编码期望值：

- tests/cpu/linalg/test_eigensolve.cpp → TestEigensolve
  - happy: `sym_power(S, 0.5)` 对 SPD 定值 `S=[[4,0],[0,9]]`（或经旋转的等价矩阵）满足 `sym_power(S,0.5)·sym_power(S,0.5) == S`，`sym_power(S,-0.5)` 复现 diag(1/2, 1/3)，容差 ≤1e-12。
  - edge: `sym_power` 对 λ_min=1e-14 的病态 S 施加 `floor=1e-10` 后返回有限值（无 inf/nan）。
  - happy: `eig_nonsymmetric` 对上三角 `A=[[2,1],[0,3]]` 返回实特征值 {2,3}（硬编码）；对含共轭复对的 2×2 返回正确实/虚部。
- tests/cpu/linalg/test_linear_solve.cpp → TestLinearSolve
  - happy: `solve_pcg` 对 SPD `A=[[4,1],[1,3]]`、`b=[1,2]`、对角预条件，解与稠密直解一致，容差 ≤1e-10，且在 max_iter 内收敛。
  - edge: `max_iter` 过小时抛出（对齐 CPHF 的非收敛抛出契约）。

dRPA 的物理数值不变量需要一次真实 SCF 求解（属 e2e），因此**不放进 `tests/`**，而由下述 regression 承担（保持单元套件秒级）。

**Regression example（本规范的永久门）**：regressions/scf/scf-review-01-linalg/reg_scf-review-01-linalg.cpp——库外消费者仅用公共 API，在 H2/STO-3G 上：(1) 走公共 RHF + direct-RPA 路径（与 regressions/scf/direct-rpa-h2-sto3g 同一公共面）复现 Furche PRB 64,195120 口径的 dRPA 相关能，容差 ≤1e-12 Ha；(2) 走 CPHF 静态极化率路径（与 regressions/scf/static-polariz-h2o 同一公共面）复现引用的各向同性 α，容差 ≤1e-6 a.u.。文献值以字面量嵌入（注释标注文献、日期），运行时不触任何第三方 oracle，≤60 s，失败非零退出。

## Out of scope
- natm-vs-n_atoms 命名统一（本轮 review 明确排除）。
- make_tau / EomImds / build_imds 去重（→ scf-review-03）；ao2mo / spin_orbital_eri 上移（→ scf-review-02）。
- FD-Hessian / GradientResult（→ scf-review-04）；所有性能项（→ scf-review-06/07）；CCSD 空间轨道重写（→ scf-review-08）。
- 泛化 `DIISMixer` 或 `vib::Hessian` 委托——本子规范不动。
- `.claude/notes` 的 5 条可复用规则捕获（units::pi、no-abort 等）——按流程属实现后 `/mol:note`，非本规范任务。
