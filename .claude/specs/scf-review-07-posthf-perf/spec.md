---
slug: scf-review-07-posthf-perf
title: "post-HF 积分/变换层性能重构（eri_dense 块装配 + RI-MP2 GEMM + ao2mo 入口拷贝消除）"
status: approved
created: 2026-07-24
domain: scf / post-hf performance
queue: scf-review
depends_on: [scf-review-02-eri-hoists]
grilled: true
---

# post-HF 积分/变换层性能重构（eri_dense 块装配 + RI-MP2 GEMM + ao2mo 入口拷贝消除）

## Summary
本 spec 消除所有 post-HF 方法共享的积分/变换层里三处 review 认定的热点：稠密 N⁴ AO ERI 张量（被 MP2/MP3/MP4/CCSD/CISD/FCI/ADC/EOM/CASSCF/GF2/G0W0/BSE 全部消费）当前用逐-AO 标量 `native::eri` 逐元素装配，对每个 (pp|pp) 四元组的 81 个 AO 积分反复重建 Boys/R-table 且无 OpenMP；RI-MP2 能量在 O(no²nv²) 幅度循环内用标量内层 P-求和并逐 aux 逐元素重建稠密 M[μν]；以及 UMP2/UCCSD 变换入口 `cur = eri_ao` 的整张 N⁴ 拷贝。重构后这三处分别落到既有壳块引擎 `native::eri_block`（摊销 Boys/R-table + OpenMP 并行）、单次 BLAS GEMM（`xt::linalg::dot`）、以及首个 quarter-transform 直接消费输入张量。对用户可见行为为零改变：稠密 ERI 张量逐元素 ≤1e-12 不变，RI-MP2 能量 ≤1e-10 Ha 不变，全部 post-HF regression 复现其引用值不变；收益是这些共享路径的墙钟时间下降（咨询性证据，不作门禁）。

## Domain basis
- 稠密 AO ERI 张量 ERI(μν|λσ) 采用 chemist 记法与 8 重置换对称性；块装配把逐-AO 标量核换成壳块核 `native::eri_block`，二者均已逐元素对标 libcint 为基准（`native_integrals.h`、`rys_eri.h`），因此两条独立实现（McMurchie–Davidson 标量路径 vs Rys/xsimd 壳块路径）在 s+p 基组上须一致到 ≤1e-12/元素——这正是等价单测的判据。
- RI-MP2：(ia|jb) ≈ Σ_P B_ia^P B_jb^P，其中 B = (μν|Q)·Wᵀ。**度量核实（scientist 开放问题 (b)）**：读 `src/cpu/scf/internal/density_fitting.cpp:182-215`，当前实现对 2-中心 Coulomb 度量 V_PQ=(P|Q) 做 Cholesky 分解 V=L·Lᵀ 并以白化因子 W=L⁻¹ 生成 B（`b_tensor()`）。由 W·Wᵀ = L⁻ᵀL⁻¹ = (L·Lᵀ)⁻¹ = V⁻¹，拟合能量中内嵌的正是 **Coulomb 度量 V⁻¹**——与对称 J^{-1/2}=V^{-1/2} 形式在能量层面等价（两者都实现同一 robust Coulomb-metric 拟合，差异仅在 FP 重结合层级）。结论：度量为 Coulomb 度量、经 Cholesky 白化实现 V⁻¹，是 Weigend–Häser 稳健拟合的正确形式。
- **冻结不变量（禁止在本 spec 改动）**：2-中心拟合度量必须保持 Coulomb 度量 V=(P|Q) 并经 W=L⁻¹ 施加 V⁻¹；把度量换成 overlap 或任何其它核属于模型改变（数量级 ~1e-3 Ha），禁止。GEMM 重写只重排 (ia|P) 收缩与 B_ov 的 BLAS 化装配，**直接复用 `df.b_tensor()` 输出的 B 张量不变**。参考：Weigend, Häser, *Theor. Chem. Acc.* 97, 331 (1997)，DOI 10.1007/s002140050269。
- RI-MP2 能量的容差为 ≤1e-10 Ha（**非机器精度**）：GEMM 与标量求和只是浮点重结合顺序不同，误差界约 N·ε ≈ 1e-10，而非 1e-15。

## Design
触及的实体全部在 `cpu/scf` 域内、均为既有函数的内部重实现，**不新增任何公共符号、不新增 `*Config` 包、不新增 factory**。

- **item 17 — 稠密 ERI 块装配**（`IntegralPack::eri_dense`，`src/cpu/scf/two_electron.cpp:17`，经 `integrals.cpp:106` 在 `compute_eri=true` 时到达）：把逐-AO `native::eri` 双重壳对循环替换为：`blocks_from_shells(shells.shells)` 分块 → 遍历规范化壳块四元组 `(A B | C D)`（外层壳块对 `left`，内层 `right ≤ left`）→ 每四元组调用一次 `native::eri_block(A,B,C,D,blk)` → 按 8 重对称把 `blk` 散射进 N⁴ 张量。构造/命名/缓冲区（`kMaxNCartOrbital^4` 栈缓冲、`n_threads()`/`thread_id()` 助手、`#pragma omp parallel for schedule(dynamic)`）严格照抄既有 Closest pattern `ShellBlockPairList::g_rhf`（`src/cpu/scf/operators.cpp:506`），使新代码读起来与既有 Fock 块构建一致。
  - **并发正确性（关键设计点，区别于 g_rhf）**：g_rhf 因是 J/K **累加**才用 thread-local G 归约；`eri_dense` 是把积分值**赋值**进互斥单元。因 `right ≤ left` 的规范枚举使每个无序壳块对 {left,right} 的 8 重对称轨道互不相交（内层 ket 壳块对索引恒 ≤ 外层 bra 索引，镜像单元 (kl|ij) 落在别的 `left` 从不覆盖的区域），故对外层 `left` 做 OpenMP 并行时各线程写入张量的**互斥单元**，无需线程局部缓冲、原子操作或归约。
  - 张量元素类型 `double`（AccReal 恒 double，符合 ERI 精度层）。s+p 之外的 l（Cartesian d/f，def2 族先在 Cartesian 引擎装配再经 `eri_to_spherical`）由 `eri_block` 的 Rys 路径原生支持，是当前标量 `native::eri`（仅 s+p）的严格超集；现役 post-HF regression 均为 STO-3G（纯 s+p），等价性判据在其上生效。
- **item 18 — RI-MP2 GEMM**（`ri_mp2_energy`，`src/cpu/scf/mp2.cpp:408`）：把 `:456-462` 的逐 aux 逐元素稠密 M[μν] 重建与 `:488-491` 的标量内层 P-求和，替换为：一次性把 B 装配成 BLAS 形状的 B_ov（(n_occ·n_vir)×n_aux），随后 `(ia|jb) = B_ov · B_ovᵀ` 用单次 `xt::linalg::dot`（dgemm）得到 (nov×nov) 的稠密收缩，OS/SS 分裂与冻结核逻辑不变。B_ov 的装配沿用既有 `xt::view` + `xt::linalg::dot`（Coᵀ·M·Cv）成语，但把逐 P 的 element-by-element M 重建替换为对 `df.b_tensor()` 列切片的 BLAS 化 reshape/视图。**`df.b_tensor()` 输出的 B 张量按位复用、拟合度量不改**。
- **item 21 折叠 + 拷贝消除**（`ao2mo_mixed_`，`src/cpu/scf/mp2.cpp:29`；`ao2mo_uhf_block`，`src/cpu/scf/ccsd.cpp:972`）：两函数首行 `xt::xtensor<double,4> cur = eri_ao;` 是整张 N⁴ 拷贝（UMP2 调 3×，UCCSD 调 4×）。重构：首个 quarter-transform（round 0）直接 `reshape_view(eri_ao,…)` 读入并写出全新工作张量，后续 round 才在 `cur` 上迭代，删去入口整拷贝。变换的数值语义与四轮收缩顺序完全不变。
- **共享路径收敛的范围界定（显式）**：mp3/mp4/cisd/fci/adc/gf2/g0w0/bse 各自的 `ao2mo_full` N⁴ 物化与标量反对称化循环（item 21 其余部分）**不在本 spec 逐个重设计**——它们已在 dep `scf-review-02-eri-hoists` 的 hoists 中收敛。本 spec 仅保证它们调用的**共享** dense-ERI（`IntegralPack::eri_dense`）与 RI 路径落在摊销/BLAS 版本；即 item-17/18 的修复落在这些消费者共享的入口处（若 02 已把入口拷贝的 `ao2mo_mixed_`/`ao2mo_uhf_block` 提升进 `internal/ao2mo.h`，则同一处编辑随其迁移到该文件，语义不变）。
- **CCSD 边界**：链上 08 负责 CCSD 空间轨道迭代重写。本 spec 只触碰 `ccsd.cpp` 的**共享入口拷贝** `ao2mo_uhf_block`，**不动**任何 CCSD 迭代内核。

### Reuse decision
（本次无 librarian_report 传入；以下候选来自 Closest-pattern 阅读，逐个裁决。）
- reuse `native::eri_block`（`include/atomiverse/cpu/scf/native_integrals.h:207`）——稠密装配调用它，不重实现壳块 ERI。
- reuse `blocks_from_shells`（`include/atomiverse/cpu/scf/internal/shell_pair.h:184`）——AO 壳分块，不新写分块逻辑。
- reuse `xt::linalg::dot`（xtensor-blas，BLAS/dgemm）——RI-MP2 收缩，不手写 GEMM。
- reuse `internal::DensityFitting::b_tensor()`（`density_fitting.h`）——B 张量与拟合度量按位复用，禁止改度量。
- pattern `ShellBlockPairList::g_rhf`（`src/cpu/scf/operators.cpp:506`）——OpenMP 结构、栈缓冲尺寸、8 重散射与助手命名照抄，使新代码风格一致。
- new — 无。本 spec 不新增任何公共或内部符号，仅重实现既有函数体。

## Files to create or modify
- `src/cpu/scf/two_electron.cpp` — 重实现 `IntegralPack::eri_dense`（块装配 + OpenMP，item 17）
- `src/cpu/scf/mp2.cpp` — 重写 `ri_mp2_energy` 为 GEMM（item 18）+ `ao2mo_mixed_` 入口拷贝消除（item 21）
- `src/cpu/scf/ccsd.cpp` — `ao2mo_uhf_block` 入口拷贝消除（item 21，仅共享入口，不动迭代内核）
- `include/atomiverse/cpu/scf/internal/integrals.h` — `IntegralPack::eri_dense` 声明处的 Doxygen（单位）
- `include/atomiverse/cpu/scf/mp2.h` — `ri_mp2_energy` 声明处的 Doxygen（单位）
- `tests/cpu/scf/test_two_electron.cpp` (new) — 块装配 ERI 张量等价单测
- `tests/cpu/scf/test_ri_mp2.cpp` — 追加 GEMM 路径 golden 单测
- `tests/CMakeLists.txt` — 注册 `test_cpu_scf_two_electron`
- `regressions/scf/posthf-eri-block-mp2/reg_posthf-eri-block-mp2.cpp` (new) — 经块装配稠密 ERI 的 H2O/STO-3G MP2 regression
- `regressions/scf/posthf-eri-block-mp2/REFERENCE.md` (new) — PySCF 固定的 golden 数值（工具/版本/命令/日期）
- `regressions/CMakeLists.txt` — 注册 `atv_add_regression(scf posthf-eri-block-mp2)`

## Tasks

- [ ] Write failing unit test `TestEriDenseBlockEquivalence` in tests/cpu/scf/test_two_electron.cpp (register test_cpu_scf_two_electron in tests/CMakeLists.txt) asserting IntegralPack::eri_dense matches a per-AO native::eri reference tensor ≤1e-12 element-wise on an STO-3G fixture containing p shells
- [ ] Reimplement IntegralPack::eri_dense in src/cpu/scf/two_electron.cpp over blocks_from_shells + native::eri_block, OpenMP over the canonical outer block-pair index, 8-fold-symmetry scatter into disjoint tensor cells (no thread-local reduction)
- [ ] Write failing unit test for the GEMM RI-MP2 path in tests/cpu/scf/test_ri_mp2.cpp asserting ri_mp2_energy reproduces the embedded golden e_corr ≤1e-10 Ha
- [ ] Rewrite ri_mp2_energy in src/cpu/scf/mp2.cpp to assemble B_ov BLAS-shaped and contract (ia|jb)=B_ov·B_ovᵀ with a single xt::linalg::dot, removing the scalar inner-P loop and per-P dense M reconstruction (df.b_tensor()/metric unchanged)
- [ ] Eliminate the entry N⁴ copy `cur = eri_ao` in ao2mo_mixed_ (src/cpu/scf/mp2.cpp) and ao2mo_uhf_block (src/cpu/scf/ccsd.cpp) by reading eri_ao directly in the first quarter-transform
- [ ] Add Doxygen for IntegralPack::eri_dense (include/atomiverse/cpu/scf/internal/integrals.h) and ri_mp2_energy (include/atomiverse/cpu/scf/mp2.h) with @brief/@param(units, Hartree/Bohr)/@return
- [ ] Add regression example regressions/scf/posthf-eri-block-mp2/reg_posthf-eri-block-mp2.cpp (+ REFERENCE.md; register in regressions/CMakeLists.txt) reproducing a PySCF-pinned H2O/STO-3G MP2 correlation energy through the block-assembled dense ERI (hard-coded golden literals, no live PySCF at runtime)
- [ ] Run full check + test suite and the post-HF `regressions` target

## Testing strategy
单测只在 `tests/`，路径镜像源码（`src/cpu/scf/two_electron.cpp` → `tests/cpu/scf/test_two_electron.cpp`；`FooClass`→`TestFooClass`），每个单测只针对单一函数/方法，禁止 e2e。
- `tests/cpu/scf/test_two_electron.cpp`（新，针对 `IntegralPack::eri_dense`）：
  - happy path — 在含 p 壳的小 STO-3G fixture（如 H2O 或单个 C/N/O 原子）上，块装配张量与本地逐-AO `native::eri` 参考张量逐元素 `max|Δ| ≤ 1e-12`（参考实现即调用仍保留的、已对标 libcint 的 `native::eri`；方法说明：以 `native::eri` 现构参考张量，不依赖已被替换的旧 `eri_dense`）。
  - edge cases — 纯 s 壳 fixture（H2/STO-3G，走 `eri_blk_ss_ss` 路径）；含 (pp|pp) 四元组（触发摊销与 8 重散射）；对称性抽查若干 (ij|kl)=(kl|ij) 等价单元。
  - domain validation — 几个硬编码 ERI 数值锚点（如 H2/STO-3G 的 (00|00)、(00|11)），防两实现同向漂移。
- `tests/cpu/scf/test_ri_mp2.cpp`（追加，针对 `ri_mp2_energy`）：GEMM 路径的 e_corr 复现嵌入的 golden `≤1e-10 Ha`（golden 从当前实现于小 fixture 上一次性记录并硬编码，守卫 item 18 仅为重结合/BLAS 化的等价改写）。
- copy-elision（item 21）为行为保持型改写，由既有 `test_cpu_scf_mp2`（UMP2 路径）与 `reg_ccsd-h2-sto3g` 守卫，不新增专用 RED 单测。
- 回归示例（唯一，`regressions/`，仓库根，不入 `tests/`）：`regressions/scf/posthf-eri-block-mp2/` —— 最小 H2O/STO-3G MP2 公共 API 脚本，经块装配稠密 ERI 复现 PySCF 固定的 MP2 相关能（literal 值内嵌，注明工具 PySCF + 版本 + 命令 + 日期；运行期**不**调用/子进程 PySCF）。
- **咨询性性能证据（非门禁）**：在一条 post-HF regression 上测量重构前后墙钟时间，写入 impl 总结，不注册为 ctest/acceptance 门禁。

## Out of scope
- 不逐个重设计 mp3/mp4/cisd/fci/adc/gf2/g0w0/bse 各自的 `ao2mo_full` N⁴ 物化与标量反对称化循环（item 21 其余部分）——它们在 dep 02 收敛；本 spec 仅保证共享 `eri_dense`/RI 路径为摊销/BLAS 版本。
- 不改 `ccsd.cpp` 的 CCSD 迭代内核（链上 08 负责空间轨道重写），仅触碰共享入口拷贝 `ao2mo_uhf_block`。
- 不改 RI 的 2-中心拟合度量（Coulomb 度量 V=(P|Q)，经 Cholesky 白化 W=L⁻¹ 实现 V⁻¹）；度量替换为模型改变，禁止。
- 不改 CUDA 后端、不改 CCSD/其它方法的公共 API 与签名、不新增公共符号。
- 性能墙钟前后对比不作 ctest 门禁（仅 impl 总结中的咨询性证据）。
- 备选考量：曾考虑对 `eri_dense` 加 Cauchy–Schwarz 屏蔽（如 g_rhf），但稠密 post-HF 消费者需要完整 N⁴ 张量（含被屏蔽的近零元），故本 spec 只做块装配 + 并行、不做屏蔽。
