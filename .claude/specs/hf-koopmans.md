---
title: Koopmans theorem tutorial — HOMO to vertical IP with OrbitalInfoObserver
status: approved
created: 2026-07-02
---

# Koopmans theorem tutorial — HOMO to vertical IP with OrbitalInfoObserver

## Summary

在现有 RHF 教程章节（`docs/tutorial/ch02-rhf.md`）中新增一节 Koopmans 定理内容，按 Jensen 式教学脉络（概念 → 数学推导 → 冻结轨道近似 → 算法 → 性质 → 局限性）讲解如何从正则 RHF 的 HOMO 轨道能得到垂直电离势（IP ≈ −ε_HOMO），包括误差抵消机制、基组敏感性与适用边界。配套交付一个新的编号示例 `examples/11_koopmans.cpp`（用现有的 `cpu::scf::OrbitalInfoObserver` 读取 HOMO，打印 He/H₂/H₂O/NH₃/CH₄ × {STO-3G, 6-31G} 的 Koopmans 表，Ha 与 eV 双单位，并与实验 IP 定性对照）和一套 PySCF 校验回归测试：教程正文中出现的每一个数值都由引擎算出并被测试钉死。唯一的引擎侧改动是给 6-31G 基组数据表补上 He（Z=2）一行——当前 `pople_631g.h` 明确标注 `Z=2 unsupported`，缺它则 He/6-31G 列无法计算。

## Domain basis

**Koopmans 定理**（T. Koopmans, *Physica* **1**, 104 (1934), DOI [10.1016/S0031-8914(34)90011-2](https://doi.org/10.1016/S0031-8914(34)90011-2)）：对正则 RHF 轨道，在冻结轨道近似下垂直电离势 IP_i ≈ −ε_i。推导核心（教程必须完整给出）：E(N) − E_frozen(N−1; hole i) = h_ii + Σ_j (2J_ij − K_ij) = ε_i。

**误差抵消**：Koopmans 忽略轨道弛豫（使 IP 被高估）又忽略电子相关（符号相反）；在接近 HF 极限的基组下两者部分抵消；最小基组下基组误差占主导。

**参考几何**（Å，逐字嵌入教程与测试）：
- He: `He 0 0 0`
- H2: `H 0 0 0; H 0 0 0.7414`
- H2O: `O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692`
- NH3: `N 0 0 0; H 0 -0.9377 -0.3816; H 0.8121 0.4689 -0.3816; H -0.8121 0.4689 -0.3816`
- CH4: `C 0 0 0; H 0.627 0.627 0.627; H 0.627 -0.627 -0.627; H -0.627 0.627 -0.627; H -0.627 -0.627 0.627`

**PySCF 2.13.1 基准（cart=True, conv_tol=1e-12），−ε_HOMO [Ha]**：

| 体系 | STO-3G | 6-31G |
|---|---|---|
| He | 0.87603551 | 0.91412663 |
| H2 | 0.57797481 | 0.59539276 |
| H2O | 0.39123677 | 0.50136813 |
| NH3 | 0.35252132 | 0.41503216 |
| CH4 | 0.51931445 | 0.54404220 |

**RHF 总能量锚点 [Ha]**：STO-3G: He −2.8077839575, H2 −1.1166843871, H2O −74.9630231385, NH3 −55.4540871348, CH4 −39.7268336369；6-31G: He −2.8551604262, H2 −1.1267339671, H2O −75.9839744727, NH3 −56.1610211797, CH4 −40.1805127401。

单位换算：1 Ha = 27.211386245988 eV。

**实验 IP（NIST；教程必须区分绝热/垂直）**：He 24.587 eV；H₂ 15.426 eV（绝热；真实垂直值更高，约 16 eV——谨慎引用或不给数字）；H₂O 12.621 绝热 / 12.60–12.62 垂直；NH₃ 10.07 绝热 / 10.8–10.85 垂直（阳离子趋于平面化）；CH₄ 12.61 绝热 / ~13.6 垂直（Jahn–Teller 畸变）。

**基组敏感性（教程必须展示）**：H₂O 的 −ε_HOMO 从 10.6 eV（STO-3G）→ 13.6 eV（6-31G），即 0.39123677 Ha → 0.50136813 Ha。

测试与 PySCF 比较 ε_HOMO 到 **1e-8 Ha**；与实验只做定性比较，绝不设数值容差。

## Design

**引擎面（近零改动）**：不新增任何公共类。复用已验证的链路：`cpu::scf::build_integrals(frame, BasisSet, charge=0)` + `rhf_solve(ints, 1e-10, 100)`（`res.eigenvalues_alpha` 供测试直接断言 Ha 值），以及 `cpu::Evaluate` + `cpu::scf::SCFEvaluator` + `cpu::scf::OrbitalInfoObserver`（示例展示 Observer 用法，绑定到 `evaluate::Phase::finalize`，经 `scf.result_handle()` 共享状态，读 `homo_eV()`）。

**唯一引擎改动**：`include/atomiverse/cpu/scf/basis_data/pople_631g.h` 的 `Z=2` 空行替换为 He 的 6-31G `AtomBasis`（两个 s 壳层）。数据必须取自 PySCF `gto.M(atom='He', basis='6-31g', cart=True)._basis`，沿用文件既有的未归一化收缩系数约定，与 PySCF 逐位一致；同时把文件头 `Supported elements: H, C, N, O, F` 更新为含 He。注意：文件头引用的 `tools/gen_basis_data.py` 在当前树中**不存在**，故本次以 oracle 脚本输出 + 注释注明出处（PySCF 版本、basis 名）替代"重新生成"，不重建整个生成器（超出范围）。

**已核实的边界情形**：He/STO-3G 只有 1 个基函数、无虚轨道——`OrbitalInfoObserver::launch` 已显式处理（LUMO 回退为 HOMO），示例中 He 行只呈现 HOMO/IP，不展示 gap。

**示例**（`examples/11_koopmans.cpp`，单一概念：Koopmans IP 表）：对 5 体系 × 2 基组循环，Evaluate + SCFEvaluator + OrbitalInfoObserver，打印表格：体系 | 基组 | −ε_HOMO [Ha] | IP_Koopmans [eV] | 实验 IP [eV]（注明绝热/垂直）。注册在 `ATV_BUILD_RHF` 门后；`examples/README.md` 补第 11 行。

**测试分层**（遵循 `tests/CMakeLists.txt` 现有 E2E 尾部 DISABLED 模式）：
- 默认 ~16 s 闸门内：`test_cpu_scf_koopmans`（`tests/cpu/scf/test_koopmans.cpp`）— 4 个锚点：He/H₂/H₂O @ STO-3G + He @ 6-31G（后者同时端到端验证新基组行）。
- E2E 层：`test_cpu_scf_koopmans_parity`（`tests/cpu/scf/test_koopmans_parity.cpp`）— 完整 5×2 扫描，−ε_HOMO 与 E_RHF 双断言到 1e-8 Ha；无条件加入 `ATV_E2E_TESTS` 尾部列表（编译保留、默认 DISABLED、CI 的 `-DATV_RUN_E2E_TESTS=ON` 启用）。
- 钉死的参考表以 `constexpr` 结构写在 parity 测试中（单一事实来源），教程表格从它转抄；两侧互加注释交叉引用。

**Oracle**：`tests/cpu/scf/pyscf_reference_koopmans.py`（新），仿照现有 `pyscf_reference.py` 模式：`scf.RHF`, `cart=True`, `conv_tol=1e-12`，输出 −ε_HOMO 与 E_RHF 的 constexpr 块。

**教程落点**：扩写 `docs/tutorial/ch02-rhf.md`（轨道能已在该章），在 "The Result That Convinces You" 之后、"What to Take Into Chapter 3" 之前插入新节（英文正文），避免动 ch00–ch15 章节编号（章节地图归 manual-system spec 管）。内容脉络：概念（电离即抽走一个电子）→ 数学（冻结轨道推导 ε_i）→ 近似（无弛豫、无相关、误差抵消）→ 算法（引擎里如何取 ε_HOMO：`rhf_solve` 特征值 / OrbitalInfoObserver）→ 性质（5×2 表 + H₂O 基组敏感性 10.6→13.6 eV + 绝热 vs 垂直辨析）→ 局限（最小基组、开壳层/Kohn–Sham 下的失效、Jahn–Teller/结构弛豫使"实验对照"必须定性）。

## Files to create or modify

- `docs/tutorial/ch02-rhf.md` — 新增 Koopmans 一节
- `examples/11_koopmans.cpp` (new) — Koopmans 表示例
- `examples/CMakeLists.txt` — 注册 `11_koopmans`（`ATV_BUILD_RHF` 门内）
- `examples/README.md` — 表格补第 11 行
- `include/atomiverse/cpu/scf/basis_data/pople_631g.h` — 补 He (Z=2) 行 + 文件头元素列表
- `tests/cpu/scf/pyscf_reference_koopmans.py` (new) — PySCF oracle 脚本
- `tests/cpu/scf/test_koopmans.cpp` (new) — 默认闸门锚点测试
- `tests/cpu/scf/test_koopmans_parity.cpp` (new) — E2E 全扫描 + 教程数值回归
- `tests/CMakeLists.txt` — 注册两个测试；parity 加入 `ATV_E2E_TESTS` 尾部列表

## Tasks

- [ ] Write oracle script `tests/cpu/scf/pyscf_reference_koopmans.py` (RHF, cart=True, conv_tol=1e-12; emits −ε_HOMO + E_RHF constexpr block and the He 6-31G `_basis` rows for the basis-table task)
- [ ] Write failing anchor tests in `tests/cpu/scf/test_koopmans.cpp` (He/H2/H2O @ STO-3G + He @ 6-31G, −ε_HOMO vs pinned PySCF to 1e-8 Ha via `rhf_solve`) and register `test_cpu_scf_koopmans` in `tests/CMakeLists.txt`
- [ ] Write failing full-sweep parity test `tests/cpu/scf/test_koopmans_parity.cpp` (5 systems × 2 bases: −ε_HOMO and E_RHF to 1e-8 Ha; pinned table as single `constexpr` source of truth with tutorial cross-reference comment); register `test_cpu_scf_koopmans_parity` and append it to the `ATV_E2E_TESTS` tail list in `tests/CMakeLists.txt`
- [ ] Implement He (Z=2) row in `include/atomiverse/cpu/scf/basis_data/pople_631g.h` from PySCF `_basis` data (unnormalized-coefficient convention; update the Supported-elements doc line; note PySCF-2.13.1 provenance in the comment since `tools/gen_basis_data.py` is absent from the tree)
- [ ] Implement `examples/11_koopmans.cpp` (Evaluate + SCFEvaluator + OrbitalInfoObserver at `evaluate::Phase::finalize`; prints system | basis | −ε_HOMO [Ha] | IP [eV] | experimental IP [eV] with adiabatic/vertical labels; He/STO-3G row shows no gap; Doxygen `@brief` file header with units); register in `examples/CMakeLists.txt` behind `ATV_BUILD_RHF` and add row 11 to `examples/README.md`
- [ ] Write the Koopmans section in `docs/tutorial/ch02-rhf.md` (English; Jensen progression concept → math → approximation → algorithm → property → limitations; frozen-orbital derivation, error cancellation, adiabatic-vs-vertical distinction, H2O 10.6→13.6 eV basis sensitivity, Koopmans 1934 DOI; every number transcribed from the pinned table in `test_koopmans_parity.cpp`)
- [ ] Verify tutorial-number regression and gate budget: every −ε_HOMO/IP value in the ch02 section appears in the pinned parity table; default `ctest` shows `test_cpu_scf_koopmans_parity` as Disabled without `-DATV_RUN_E2E_TESTS=ON` and the anchor test adds well under ~2 s
- [ ] Run full check + test suite

## Testing strategy

- **正路**：4 个锚点（He/H₂/H₂O @ STO-3G，He @ 6-31G）在默认闸门中断言 −ε_HOMO 与钉死的 PySCF 值差 < 1e-8 Ha；E2E 层做完整 5×2 扫描并同时断言 RHF 总能量（护住"教程数字来自引擎"的锚）。
- **边界**：He/STO-3G 无虚轨道 → 断言 `OrbitalInfoObserver::lumo_eV() == homo_eV()`（已实现的回退语义）且示例不打印该行 gap；He/6-31G 在基组行落地前必须抛出（无回退值规则）——即锚点测试先 RED。
- **领域校验**（`science.required`）：eV 换算用 1 Ha = 27.211386245988 eV 显式断言一处；与实验 IP 仅定性对照（教程文字），测试不含实验数值容差；oracle 脚本可重跑再生全部钉死值（cart=True 与引擎笛卡尔 GTO 约定一致）。
- **闸门纪律**：parity 测试编译保留但默认 DISABLED（`ATV_E2E_TESTS` 尾部机制），CI（`ATV_RUN_E2E_TESTS=ON`）启用。

## Out of scope

- 重建缺失的 `tools/gen_basis_data.py` 生成器（本次仅手补 He 一行并注明出处；生成器复活是独立的工具债）。
- `examples/README.md` 中缺失的 09/10 两行（属于既有文档滞后，不在本 spec 顺手修）。
- ΔSCF / EOM / 电子亲和能（−ε_LUMO 对 EA 的更差近似只在"局限性"一段提及，不做数值交付）；UHF/开壳层与 Kohn–Sham 轨道能的 Koopmans 失效同样只述不算。
- 新教程章节或 ch00–ch15 重编号（章节地图归 manual-system spec）。
- CUDA 侧任何改动（`cuda::scf::SCFEvaluator` 无梯度亦无此教程需求）。
