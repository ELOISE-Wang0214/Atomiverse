---
title: AO→MO integral transformation and closed-shell MP2 correlation energy
status: approved
created: 2026-07-02
---

# AO→MO 积分变换与闭壳层 MP2 相关能

## Summary

为 Atomiverse 的高斯基 SCF 域新增第一个后 HF（post-SCF）相关方法：AO→MO 双电子积分变换（四步 O(N^5) 四分之一变换，只构建 MP2 所需的 (ia|jb) occ×vir×occ×vir 块）与闭壳层 MP2 相关能 `mp2_energy`。二者是 `ATV::cpu::scf` 命名空间下的自由函数，消费已收敛的 `SCFResult` 与显式物化的 `IntegralPack::ERI`——与 `scf_gradient` 同型的"SCF 结果的纯函数"，不引入任何 Driver/Component。同时把 docs/tutorial/ch07-correlation.md 从"MP2 未实现"的状态章改写为真正的实战章 "Electron correlation beyond Hartree–Fock"（单行列式 → 缺失的相关能 → 占据-虚轨道微扰激发 → 后 HF 为何依赖 HF 输出质量，Jensen 风格），配套新示例 `examples/13_mp2.cpp`。正确性以 PySCF 2.13.1 固定参考值为准：默认测试层验证 H2O/STO-3G 与 He/STO-3G 边界情形，E2E 层覆盖全部 10 个 (体系, 基组) 组合，容差 1e-8 Ha。本 spec 独立成链尾；若 hf-robust-02 先落地，则沿用其 n_mo 语义（软依赖，见 Design）。

## Domain basis

**MP2 相关能**（化学家记号，实系数正则空间轨道，闭壳层）：

$$ E_\mathrm{MP2} = \sum_{ij}^{\mathrm{occ}} \sum_{ab}^{\mathrm{vir}} \frac{(ia|jb)\,[\,2(ia|jb) - (ib|ja)\,]}{\varepsilon_i + \varepsilon_j - \varepsilon_a - \varepsilon_b} $$

分母恒为负（需断言 ε_a − ε_i > 0），故 E_MP2 < 0；闭壳层基态无需简并分母保护。

**AO→MO 变换**：四步顺序 O(N^5) 四分之一变换，占据指标优先以缩小中间体：
(μν|λσ) → (iν|λσ) → (ia|λσ) → (ia|jσ) → (ia|jb)。禁止 O(N^8) 直接缩并。

**约定**：全电子（对齐 PySCF MP2 默认 frozen=None，mp/mp2.py:522）；n_vir = n_mo − n_occ，其中 n_mo 可能因线性相关剔除而小于 n_ao。

**参考文献**：
- Møller, C.; Plesset, M. S. *Phys. Rev.* **46**, 618 (1934). DOI 10.1103/PhysRev.46.618
- Head-Gordon, M.; Pople, J. A.; Frisch, M. J. *Chem. Phys. Lett.* **153**, 503 (1988). DOI 10.1016/0009-2614(88)85250-3（四分之一变换）

**固定参考值**（PySCF 2.13.1，`cart=True`，SCF `conv_tol=1e-12`；E_corr，单位 Hartree；几何与 hf-koopmans spec 相同：H2 键长 0.7414 Å；H2O O(0,0,0.1173)/H(0,±0.7572,−0.4692) Å；NH3、CH4 标准几何；He 单原子）：

| 体系 | STO-3G | 6-31G |
|---|---|---|
| He | **0.0（严格为零——无虚轨道，必须返回 0 而非 NaN）** | −0.0112001229 |
| H2 | −0.0131707665 | −0.0173964438 |
| H2O | −0.0355456516 | −0.1288509172 |
| NH3 | −0.0472198510 | −0.1169596914 |
| CH4 | −0.0561492533 | −0.0996598545 |

对比容差 1e-8 Ha。参考值再生：`mol = gto.M(atom, basis, cart=True); mf = scf.RHF(mol); mf.conv_tol = 1e-12; mp.MP2(mf).kernel()[0]`。

## Design

**放置**：`cpu/scf` 域——新建 `include/atomiverse/cpu/scf/mp2.h` + `src/cpu/scf/mp2.cpp`（头/实现分离，镜像 `gradient.h`/`gradient.cpp`）。按仓库规则 "SCF 是内层求解器"，MP2 与 `scf_gradient` 同类，是后 SCF 分析纯函数——**不是** Driver 子类，**不是** Component，不注册进任何 Plan。

**新公共符号**（均在 `ATV::cpu::scf`）：

- `xt::xtensor<double, 4> ao2mo_ovov(const xt::xtensor<double, 4>& ERI, const xt::xtensor<double, 2>& C, int n_occ, int n_mo)` — 只返回 (ia|jb) 块，形状 (n_occ, n_vir, n_occ, n_vir)；**不**构建完整 MO 张量。每步四分之一变换 reshape 成矩阵后用 xt::linalg（BLAS-backed dot/tensordot）缩并——仓库有 `-march=native` + OpenBLAS，禁止裸标量四重循环。全程 `AccReal`=double。
- `struct MP2Result { double e_corr; double e_tot; };` — 单位 Hartree；`e_tot = scf.energy + e_corr`。
- `MP2Result mp2_energy(const IntegralPack& ints, const SCFResult& scf)` — 从 `scf.C_alpha`、`scf.eigenvalues_alpha`、`scf.n_alpha` 取输入；`n_occ = scf.n_alpha`，**`n_mo = C_alpha.shape(1)`**（不是 `ints.n_basis`——线性相关剔除后 n_mo < n_ao；这是对 hf-robust-02 n_mo 语义的软依赖，即便其未落地此写法也正确）。

**生命周期/所有权**：无状态自由函数；(ia|jb) 中间体由调用返回值持有，函数内不缓存。Incore-only：范围内体系极小（H2O/6-31G n_ao=13 → ERI 13⁴ double ≈ 0.2 MB）。

**守卫**（无回退值规则——一律 `std::invalid_argument` 抛出，不返回哨兵）：
- `ints.ERI` 未物化（`ERI.size() == 0`，即调用方忘记 `build_integrals(..., /*compute_eri=*/true)`）→ 抛出，消息指明需 `compute_eri=true`；
- `!scf.converged` → 抛出；
- 非限制性结果（`scf.n_alpha != scf.n_beta`）→ 抛出（本 spec 仅闭壳层）；
- 分母符号断言：所有 ε_a − ε_i > 0；
- **He/STO-3G 边界**：n_vir == 0 → 直接返回 `{0.0, scf.energy}`，不进入循环、不产生 NaN。

**示例**：`examples/13_mp2.cpp`（接续 12_delta_scf 之后；11/12 由 hf-koopmans/hf-delta-scf 姊妹链提供，若其尚未落地则暂存编号空洞，不阻塞本 spec），打印 E_HF、E_corr、E_MP2 三行；在 `examples/CMakeLists.txt` 的 `if(ATV_BUILD_RHF)` 块注册。

**教程**：重写 `docs/tutorial/ch07-correlation.md`（保留章号 7，英文内容）：替换现第 11–17、93–103、134–141 行等 "not implemented / status chapter" 叙述为真实章节——单行列式的局限 → 相关能定义（沿用现有开头）→ MP2 微扰推导概要与 (ia|jb) 变换 → 为何后 HF 质量取决于 HF 输出（轨道与轨道能是零阶量）→ 13_mp2 可运行示例与参考值表。章内引用的 E_corr 数字必须与测试中固定的常量逐位一致（教程回归）。

**测试**：`tests/cpu/scf/test_mp2.cpp`，按命名约定注册为 `test_cpu_scf_mp2`。默认 ~16 s 门槛层只跑 H2O/STO-3G 平价 + He 边界 + (ia|jb) 对称性 + 守卫抛出（毫秒级）；全部 10 组 (体系, 基组) 扫描用 GoogleTest 独立测试用例置于 `ATV_RUN_E2E_TESTS` 门后（沿用 tests/CMakeLists.txt 尾部的 DISABLED 列表机制）。

## Files to create or modify

- `include/atomiverse/cpu/scf/mp2.h` (new) — `ao2mo_ovov`、`MP2Result`、`mp2_energy` 声明 + Doxygen
- `src/cpu/scf/mp2.cpp` (new) — 实现
- `CMakeLists.txt` — 在 `ATV_BUILD_RHF` 源文件列表（约第 316 行 `src/cpu/scf/gradient.cpp` 旁）加入 `src/cpu/scf/mp2.cpp`
- `tests/cpu/scf/test_mp2.cpp` (new) — 默认层 + E2E 层测试
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_mp2`；E2E 扫描用例挂入 `ATV_RUN_E2E_TESTS` 门
- `examples/13_mp2.cpp` (new) — 打印 E_HF / E_corr / E_MP2
- `examples/CMakeLists.txt` — `if(ATV_BUILD_RHF)` 块内注册 `13_mp2`
- `docs/tutorial/ch07-correlation.md` — 由状态章改写为实战章

## Tasks

- [ ] Write failing tests for ao2mo_ovov + mp2_energy default tier (tests/cpu/scf/test_mp2.cpp): H2O/STO-3G E_corr = −0.0355456516 Ha @1e-8, He/STO-3G returns exactly 0.0 (finite, no NaN), (ia|jb)=(jb|ia) symmetry @1e-12, guard throws (ERI absent / unconverged / unrestricted); register `test_cpu_scf_mp2` in tests/CMakeLists.txt
- [ ] Write failing E2E-tier sweep tests in tests/cpu/scf/test_mp2.cpp covering all 10 (system, basis) pairs from the Domain-basis table @1e-8 Ha, gated behind `ATV_RUN_E2E_TESTS` in tests/CMakeLists.txt
- [ ] Implement ao2mo_ovov in include/atomiverse/cpu/scf/mp2.h + src/cpu/scf/mp2.cpp (four BLAS-shaped quarter transforms via xt::linalg, occupied indices first, ovov block only); add src/cpu/scf/mp2.cpp to the ATV_BUILD_RHF source list in CMakeLists.txt
- [ ] Implement MP2Result + mp2_energy in include/atomiverse/cpu/scf/mp2.h + src/cpu/scf/mp2.cpp (n_mo = C_alpha.shape(1), n_occ = scf.n_alpha, denominator sign assert, n_vir==0 → 0.0, invalid_argument guards)
- [ ] Add docstring per doxygen style with units (Hartree) for every public symbol in include/atomiverse/cpu/scf/mp2.h (`@brief`/`@param`/`@return`)
- [ ] Implement examples/13_mp2.cpp printing E_HF, E_corr, E_MP2 for H2O/STO-3G and H2O/6-31G; register in examples/CMakeLists.txt under ATV_BUILD_RHF
- [ ] Write the replacement chapter docs/tutorial/ch07-correlation.md ("Electron correlation beyond Hartree–Fock", Jensen-style, keeps chapter number, embeds the 13_mp2 walk-through, MP2 equation, both DOIs, and the pinned reference table)
- [ ] Verify tutorial-regression coupling: assert in tests/cpu/scf/test_mp2.cpp the exact E_corr constants quoted in ch07 (H2O/STO-3G and H2O/6-31G) so chapter numbers cannot drift from code
- [ ] Verify default ctest gate: measure `test_cpu_scf_mp2` default-tier runtime stays ms-scale and the default `ctest` wall time keeps the ~16 s budget; run E2E sweep once with `-DATV_RUN_E2E_TESTS=ON` (10/10 pass)
- [ ] Run full check + test suite

## Testing strategy

- **主路径（默认层）**：H2O/STO-3G — `build_integrals(frame, BasisSet::STO3G, 0, /*compute_eri=*/true)` + `rhf_solve(ints, 1e-12, 100)` + `mp2_energy`，E_corr 对 PySCF 固定值 −0.0355456516 Ha 容差 1e-8；`e_tot == scf.energy + e_corr`。
- **边界情形（默认层）**：He/STO-3G 无虚轨道 → E_corr 严格 == 0.0 且 `std::isfinite`（NaN 回归守卫）；ERI 未物化 / SCF 未收敛 / 非限制性输入 → `std::invalid_argument`。
- **单元不变量（默认层）**：H2O/STO-3G 的 (ia|jb) 张量满足 (ia|jb) = (jb|ia)，逐元素 ≤ 1e-12；E_MP2 < 0（H2 以上体系）。
- **领域验证（E2E 层，`ATV_RUN_E2E_TESTS=ON`）**：{He, H2, H2O, NH3, CH4} × {STO-3G, 6-31G} 全部 10 组对 Domain-basis 表中固定 PySCF 值 @1e-8 Ha。
- **教程回归**：测试断言 ch07 章内引用的 E_corr 数字（与测试固定常量同源），防止文档与代码漂移。
- **示例烟测**：`13_mp2` 目标在 `ATV_BUILD_RHF=ON` 下编译并运行，输出含 E_HF、E_corr、E_MP2。
- **门槛保护**：默认 `ctest` 总时长维持 ~16 s 预算；新增默认层用例为毫秒级。

## Out of scope

- 开壳层 MP2（UMP2 / ROMP2）——本 spec 仅闭壳层限制性结果，非限制性输入显式抛出。
- Frozen-core 近似（对齐 PySCF 默认全电子；冻芯是后续独立 spec）。
- Direct MP2 / DF-MP2（RI）/ 半直接积分变换——范围内体系极小，incore 足够；曾考虑 DF 路线，因需引入辅助基组基础设施而明确排除。
- 完整 MO 基 ERI 张量（`ao2mo_full`）——MP2 只需 ovov 块；考虑过通用变换，按"最窄接口"原则排除。
- MP2 解析梯度、弛豫密度、自然轨道；SCS-/SOS-MP2 标度变体。
- CUDA 后端 MP2（`cuda::scf` 无对应；CUDA SCFEvaluator 本身无梯度，平价路线未建立）。
- 任何 Driver/Component 包装（如 "MP2Evaluator"）——ch07 旧章曾勾画该形态，本 spec 明确采用 `scf_gradient` 式纯函数；组件化留待真实 AIMD-on-MP2 需求出现。
