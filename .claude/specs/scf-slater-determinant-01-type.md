---
title: SCF Slater 行列式 · 01 类型（cpu::scf::SlaterDeterminant 值类型 + 单测）
slug: scf-slater-determinant-01-type
status: code-complete
chain: scf-slater-determinant
link: 01/04
created: 2026-07-20
---

# SCF Slater 行列式 · 01 类型

## Summary
新增公共值类型 `cpu::scf::SlaterDeterminant`（`include/atomiverse/cpu/scf/slater_determinant.h`，头文件实现），承载均场波函数的最小完备表示：每自旋道 {MO 系数 C、轨道能 ε、占据向量 f}，受限（RHF/RKS）单道复用（β 访问器即 α 视图）、非受限（UHF/UKS）α/β 双道。按自旋道存储是拍板约定——P^σ 的幂等不变量最干净，factor-2 闭壳层总密度只是派生值。成员方法给出密度矩阵（委托 `linalg::build_density_matrix`，GPW 共用、本轮不改名）、单电子期望值、HOMO/LUMO/轨道能访问器（吸收 `orbital_info_observer.cpp:43-60` 的手工逻辑，observer 迁移在 03），以及 ΔSCF/MOM 入口 `with_occupation`。构造即校验不变量（违规 `throw std::invalid_argument` 带上下文、无 fallback，`internal::SpinResolution` 先例）；需要重叠矩阵 S 的不变量（Tr(P^σS)=N^σ、幂等）作 `validate(S)` 方法而非构造强制。**本环零消费者接线：现有生产代码一行不改**，仅新头 + 单测 + 测试注册；02/03/04 依次消费本类型。

## Domain basis
（scientist 已验证，逐字采纳）

1. 闭壳层密度 P_μν = 2Σ_a^{N/2} C_μa C_νa；自旋分辨（Pople–Nesbet）P^α = C^α_occ (C^α_occ)ᵀ、P_total = P^α + P^β。Szabo & Ostlund, *Modern Quantum Chemistry* (Dover 1996), Ch. 3。与 `linalg/density.h`、`mean_field_method.cpp:388-401` 现行代码零偏差。
2. ⟨Φ|Ô|Φ⟩ = Tr(P·O) **仅对单电子算符成立**；双电子 ⟨Φ|Σg|Φ⟩ = ½Σ_ij(⟨ij|ij⟩−⟨ij|ji⟩) 不可写为固定 Tr(P·G)——`expect()` 文档必须写明。
3. E_HF = ½Tr[P(H_core+F)] + E_nn；Σε_i = Tr(PF) 双计电子互斥——**均场能量不是 Fock 期望值**（`expect()` 文档警告的依据）。
4. Koopmans：IP ≈ −ε_HOMO（冻结轨道近似）。T. Koopmans, *Physica* 1, 104–113 (1934), DOI 10.1016/S0031-8914(34)90011-2。
5. ΔSCF/MOM：非 Aufbau 占据 + 重收敛 = ΔSCF；MOM 按与上一迭代占据轨道最大重叠（O = C_oldᵀ S C_new）选占据防塌缩。Gilbert, Besley, Gill, *JPCA* 112(50), 13164–13171 (2008), DOI 10.1021/jp801738f。`with_occupation` 即标准入口。
6. 不变量：Tr(P^σS) = N^σ（广义为 Σf^σ）；P^σSP^σ = P^σ 当且仅当 f∈{0,1}（单行列式判别）；Pauli 界 0 ≤ f ≤ 1 每自旋道。factor-2 约定下 PSP = 2P——按道存储回避此歧义。

单位：ε 存 Hartree（SCF 内部恒 Hartree，eV 换算属消费方——observer 现即如此）；C、f、AO 基 P 无量纲；`expect(O)` 返回值单位随 O。精度域 AccReal/double（SCF 恒 double）。

## Design

**放置（librarian 裁定）：** `cpu/scf/` 顶层公共头——`internal/` 会经 `MeanFieldOutputs` 泄漏进 recorder 公共头（禁止），`stages/` 只放 Component。值类型先例 `cpu::Box`（`cpu/frame/box.h:46`：值类型 + 不变量缓存 + mutate-through-method）；命名值类型名词、裸名词访问器。矩阵类型直接复用 `dense_matrix.h` 的 `scf::RealMatrix`/`RealVector`（= `linalg::` 别名）。

**存储：** 内嵌 `enum class Spin { alpha, beta }`（树内无现成自旋道枚举，`scf_result.h:24` 的 `SpinType` 是受限/非受限判别、语义不同，不复用不混淆）。每道一个 `Channel { RealMatrix C; RealVector epsilon; RealVector occupation; }`。受限形态只存 α 道，β 访问全部转发 α（视图语义，不复制）；非受限形态存双道。`restricted()` 访问器暴露形态。

**构造与校验（`internal::SpinResolution` 先例：校验构造 + throw + 上下文）：**
- 受限构造 `SlaterDeterminant(C, epsilon, occupation, n_electrons_per_spin)`；非受限构造 `SlaterDeterminant(C_a, eps_a, f_a, n_alpha, C_b, eps_b, f_b, n_beta)`。
- 构造即校验（违规 `throw std::invalid_argument` 带具体数值上下文，无 fallback）：维度一致（C 列数 = ε 长 = f 长，C 行数 = AO 维 n）；Pauli 界 0 ≤ f_i ≤ 1（自旋轨道）；Σf_i^σ = N^σ（|Σf−N| ≤ 1e-12·max(1,N)）。
- 带 S 的不变量**不在构造强制**（构造时无 S）：`validate(const RealMatrix& S, double tol = 1e-8) const` 校验 Tr(P^σS) = N^σ（每道）；当该道 f 全为 {0,1} 时另校验幂等 P^σSP^σ = P^σ（分数占据道跳过幂等——Domain basis 第 6 条，幂等当且仅当整占据）。违规 throw 带道别与偏差量。

**方法：**
- `density(Spin)` → P^σ = Σ_i f_i^σ c_i c_iᵀ。实现：收集 f_i>0 的列、按 √f_i 加权拼为 C_occ'，委托 `linalg::build_density_matrix(C_occ', n_occ', n, 1.0)`（`include/atomiverse/cpu/linalg/density.h`，GPW 共用，本轮不改名——04-methods 里门禁临时豁免它；`internal/density_builder.h::DensityBuilder::step` 现在就是这么委托的）。Aufbau 整占据时数值上与现行 `C_occ C_occᵀ` GEMM 逐位一致；非 Aufbau（`with_occupation` 之后）与分数占据同一条路径覆盖。
- `density_total()` → P^α + P^β；受限即 2·C_occ C_occᵀ（factor-2 是派生值，不是存储约定）。
- `expect(const RealMatrix& O)` → Tr(P_total·O)。Doxygen 必须写明：仅单电子算符成立（双电子需二粒子密度矩阵）；均场总能量不是 `expect(F)`（Σε 双计互斥，E_HF = ½Tr[P(H_core+F)]+E_nn）。
- `orbital_energy(Spin, i)` / `homo(Spin)` / `lumo(Spin)` [Ha]：吸收 `orbital_info_observer.cpp:43-60` 的手工逻辑——以占据向量定义（HOMO = 最高占据轨道 ε，LUMO = 最低空轨道 ε；ε 升序 Aufbau 情形与 observer `eps(n_occ-1)`/`eps(n_occ)` 逐位一致）；无空轨道时 `lumo()` 返回 `homo()`（observer:57-59 现行语义）；越界 i throw。道参数显式（按道存储的诚实 API；受限两道同值）。本 spec 只实现方法，observer 迁移在 03。
- `with_occupation(Spin, from, to)`：ΔSCF/MOM 入口——要求 f_from = 1 且 f_to = 0（整占据转移，违规 throw），返回修改占据后的**新值**并重新过构造校验。ref-qualified &/&& 链式对（`mean_field_method.h:180-187` `with_gradient` 先例）：`const &` 重载拷贝后修改返回新值，`&&` 重载原地修改后移动返回。

**Reuse decision：**
- `reuse` `linalg::build_density_matrix`（density 委托，不改名不改签名）；`scf::RealMatrix`/`RealVector` 别名；`internal::SpinResolution` 的校验构造风格。
- `new` `slater_determinant.h` + `test_slater_determinant.cpp`；无 `make_*` 工厂、无 `*Config`、无新 `.sh` 门。
- **不动**：`scf_result.h`（02/03）、`orbital_info_observer.*`（03）、`mean_field_method.*` 与一切 `build_*` 转换（04）、`linalg/density.h`（04 豁免后处理）。

## Files to create or modify
- `include/atomiverse/cpu/scf/slater_determinant.h` (new) — 值类型全部实现（头文件内联，值类型 + 校验构造，`cpu::Box` 先例）。
- `tests/cpu/scf/test_slater_determinant.cpp` (new) — 每方法一测。
- `tests/CMakeLists.txt` — 注册 `atv_add_cpp_test(test_cpu_scf_slater_determinant cpu/scf/test_slater_determinant.cpp)`（HF 测试段，`test_cpu_scf_spin_resolution` 旁）。

## Tasks
- [x] Write failing tests for SlaterDeterminant (tests/cpu/scf/test_slater_determinant.cpp): ctor validation throws, density vs direct GEMM, density_total both forms, expect vs manual trace, homo/lumo/orbital_energy on known ε, with_occupation &/&& pair, validate(S) on orthonormal toy basis; register test_cpu_scf_slater_determinant in tests/CMakeLists.txt
- [x] Implement SlaterDeterminant storage + validating constructors in include/atomiverse/cpu/scf/slater_determinant.h (per-spin Channel, restricted β-as-α view, nested Spin enum, throw std::invalid_argument with context)
- [x] Implement density(Spin) and density_total() delegating to linalg::build_density_matrix (√f-weighted occupied-column gather; Aufbau path bit-identical to current GEMM)
- [x] Implement expect(O), orbital_energy(Spin,i), homo(Spin), lumo(Spin) (occupation-defined; no-virtual edge returns homo, observer:57-59 semantics)
- [x] Implement with_occupation(Spin,from,to) as ref-qualified const&/&& pair returning a re-validated new value (with_gradient precedent)
- [x] Implement validate(S,tol): Tr(P^σS)=N^σ per channel; idempotency P^σSP^σ=P^σ only for integer occupations
- [x] Add docstring per doxygen with units: @brief/@param/@return on every public symbol; ε [Ha]; expect() one-electron-only + energy-is-not-expectation warnings
- [x] Verify zero-consumer wiring: grep SlaterDeterminant over include/ and src/ matches only slater_determinant.h; no existing production file modified
- [x] Run full check + test suite
- [x] Hygiene cleanup — /mol:simplify ran clean (janitor: zero findings, 2026-07-20)

## Testing strategy
- **构造校验（happy + throw）**：合法受限/非受限构造成功；f<0、f>1、Σf≠N、维度不一致各自 `std::invalid_argument` 且消息含具体数值。
- **density 一致性**：3×3 玩具 C 上 `density(alpha)` 与手写 `xt::linalg::dot(C_occ, transpose(C_occ))` 逐元素 ≤1e-12；非 Aufbau 占据（1,0,1）与分数占据（√f 加权）各一例对手工 Σf_i c_i c_iᵀ。
- **density_total 两形态**：受限 = 2·density(alpha)；非受限 = P^α+P^β（≤1e-12）。
- **expect**：`expect(O)` = 手工 Tr(P_total·O)（≤1e-12）。
- **homo/lumo/orbital_energy**：已知 ε 序列（升序）Aufbau 下与 `eps(n_occ-1)`/`eps(n_occ)` 一致；无空轨道时 lumo==homo；越界 throw。
- **with_occupation**：换道后新值 f 正确、原值不变（const& 路径）；&& 链式路径可用；f_from≠1 或 f_to≠0 throw。
- **受限/非受限形态**：受限 β 访问器与 α 同值（density(beta)==density(alpha)、homo(beta)==homo(alpha)）。
- **validate(S)（domain validation）**：正交归一玩具基（S=I）下 Tr(PS)=N 与幂等 PSP=P 通过（整占据）；分数占据跳过幂等仍通过 Tr(PS)=N；非正交归一 C 配 S=I 时 throw。即 Domain basis 第 6 条不变量的直接验证。
- **速度**：全部玩具矩阵（n≤4），无 SCF 全解——整套秒级，E2E 禁入 `tests/`。
- **Regression example**：本环不新增 `regressions/` 案例——纯类型引入、零消费者，行为面在 03/04 接线后由既有 SCF 回归守卫。

## Out of scope
- 一切消费者接线（observer/recorder/stages 迁移）——03。
- `scf_result.h` 改动（`SlaterDeterminant` 入 result）——02/03。
- `build_*` 转换与 `linalg::build_density_matrix` 改名/门禁豁免处理——04-methods。
- `Wavefunction` 基类——单实现不归纳（已定案）。
- 算符代数/表达式模板——性能 + 物理双重否决（已在讨论定案）。
- MOM 的最大重叠选择算法本体（O=C_oldᵀSC_new 逐迭代重排）——`with_occupation` 只是占据入口，重收敛与 MOM 循环属后续方法层。
