---
slug: scf-review-06-scf-perf
title: "scf-review-06-scf-perf — KS/杂化 SCF 热路径去冗余（Phi 缓存 · 融合 JK · ERI 微分配）"
status: done
created: 2026-07-24
domain: scf / performance
queue: scf-review
depends_on: [scf-review-05-api-hygiene]
grilled: true
---

# scf-review-06-scf-perf — KS/杂化 SCF 热路径去冗余（Phi 缓存 · 融合 JK · ERI 微分配）

## Summary

本 spec 针对均场/DFT 的 Fock 构建热路径，消除四类跨 SCF 迭代重复的冗余工作，且**保持数值逐位不变**：（1）把几何固定的 AO-on-grid 矩阵 `Phi[G×N]`（GGA 另含 `dPhi_x/y/z`）从"每次 KS 迭代重建"改为"每几何构建一次、跨迭代复用"；（2）把受限杂化（B3LYP/PBE0）中"先 `g_rhf` 全扫一遍再 `k` 全扫一遍"的两遍 ERI 扫描融合为**一遍**同时产出 J 与 K；（3）修掉 ERI 内层三处微分配/串行热点（整块 80KB 清零、串行 K 收缩、每四元组堆分配）。对用户可见的结果（`RHF/UHF/RKS/UKS` 能量、力、Vxc）不变，仅缩短 def2-SVP 苯量级 KS/杂化单点与 AIMD 的每步壁钟时间。所有改动落在 `cpu/scf` 单一 domain 内，依赖 `scf-review-05-api-hygiene` 先序列化 rys 表面改动。

## Domain basis

受限全局杂化 KS 的 Fock 矩阵（内部单位 Hartree，公开面 eV）：

- `F = h_core + J[P] − ½·a_x·K[P] + V_xc[ρ]`
- Coulomb：`J_{μν} = Σ_{λσ} P_{λσ} (μν|λσ)`
- 精确交换：`K_{μν} = Σ_{λσ} P_{λσ} (μλ|νσ)`

J 与 K 由**同一批**两电子积分 `(μν|λσ)` 收缩得到，故一遍扫描同时累加 J 与 K 在代数上与两遍扫描**严格相等**——线性性保证结果一致，仅浮点加法次序不同（因此验收用逐元素 `‖F_fused − F_twopass‖∞ ≤ 1e-12` 而非等式）。杂化系数 `a_x`：B3LYP `a_x = 0.20`，PBE0 `a_x = 0.25`。两电子积分经 Rys 求积块内评估（`eri_block`）。KS 的 `V_xc` 由分子求积网格上的 `ρ = rowsum((Phi·P)⊙Phi)` 与 GEMM 组装，网格随几何固定，仅密度逐迭代变化。

参考：
- Almlöf, Faegri, Korsell, direct SCF，J. Comput. Chem. 3, 385 (1982), DOI 10.1002/jcc.540030314（一遍 JK 的直接-SCF 依据）。
- Dupuis, Rys, King, J. Chem. Phys. 65, 111 (1976), DOI 10.1063/1.432807（Rys 求积 ERI）。
- Becke, J. Chem. Phys. 98, 5648 (1993), DOI 10.1063/1.464913；Stephens 等, J. Phys. Chem. 98, 11623 (1994), DOI 10.1021/j100096a001（B3LYP，a_x=0.20）。
- Adamo, Barone, J. Chem. Phys. 110, 6158 (1999), DOI 10.1063/1.478522（PBE0，a_x=0.25）。

## Design

涉及实体：`evaluate_xc_batched`（`xc_grid_engine.h` 内联引擎，被 `hamiltonian_builder.cpp`、`tdhf.cpp`、gradient 共享）、`internal::HamiltonianBuilder`、`ShellBlockPairList`、`IntegralPack::coulomb_exchange_fused`、`rys_eri_block_canonical` / `rys_eri_block`。全部为方法/成员（无自由函数工厂、无 `*Config` 包），符合默认 OOP 形状。

**1. 几何固定的 AO-on-grid 缓存（Item 15）。**
- 新增内部工作区类型 `XCGridBasis`（置于 `xc_grid_engine.h`，非 `*Config` 包、非"构建计算"的中间世界模型；带动词方法 `build`）：持有 `xt::xtensor<double,2> Phi, dPhi_x, dPhi_y, dPhi_z;` 与 `bool gga=false, built=false;`，`build(const MolecularGrid&, const std::vector<native::Shell>&, int n_basis, bool is_gga)` 复用现有 `eval_basis` / `eval_basis_val_grad`（与 cube 写出同一码路），把当前 `evaluate_xc_batched` 第 129–176 行的逐点 Phi/dPhi 构建搬入。
- `evaluate_xc_batched` 拆为两个重载：新增 `evaluate_xc_batched(const XCGridBasis& aog, const MolecularGrid& grid, P_alpha, P_beta, xc_type)` 执行原步骤 2–5（密度 GEMM → XC batch → Vxc GEMM）；保留原签名作**薄包装**（内部构造临时 `XCGridBasis` 再委托），使 `tdhf.cpp` 与 gradient 调用点逐位不变。
- `HamiltonianBuilder` 增私有 `mutable XCGridBasis xc_basis_;`。`prepare()`（每几何一次、重建 `grid_`）令 `xc_basis_.built = false` 使其失效；`step()`（`const`）在首个 KS 迭代按 `!xc_basis_.built || xc_basis_.gga != is_gga(xc_type)` 惰性 `build`，其后各迭代直接复用，调用新重载。
- **内存策略（无条件缓存，不设阈值）**：def2-SVP 苯（C6H6，n_ao_cart ≈ 120，中等网格 G ≈ 1.2×10^5）下 `Phi` ≈ 115 MB，GGA 的 `Phi+3·dPhi` ≈ 460 MB。理由：原 `evaluate_xc_batched` 每次调用本就把 `Phi`（+3 `dPhi`）连同 `T_a,T_b`（2·G·N）作为局部量分配，单次 Fock 构建的峰值工作集已 ≥ 此；缓存只是让同一批缓冲跨 ~15–30 次 SCF 迭代常驻，省去每次 G×N 次 `exp/pow` 超越函数评估。稳态净增 1×G·N（LDA）/ 4×G·N（GGA），与已存在的 `T_a/T_b` 瞬态同阶。在头文件注释中记录 GGA ×4 成本；若未来超大网格 OOM，逃生口是按网格点数设阈值——**当前不引入任何投机性配置（禁 `*Config` 包）**。
- 同文件两处 🟢：删除只写不读的 `std::vector<char> active`（原 :142/:168/:173，最小修复即删除死写；筛选接线属越界优化）；把 `P→double`（:180–185）与 `Vxc→double`（:351–355）逐元素拷贝用 `if constexpr (std::is_same_v<RealMatrix::value_type, double>)` 守卫——当元素类型已是 `double`（默认 `AccReal`）时跳过 O(N²) 拷贝、直接引用/移动，否则保留转换（为将来 `StorageReal=float` 留路）。

**2. 融合一遍 JK（Item 16）。**
- `ShellBlockPairList` 新增成员 `void jk_rhf(const RealMatrix& P_total, RealMatrix& J, RealMatrix& K) const;`（声明于 `shell_pair.h`，实现于 `operators.cpp`）：在**同一** `eri_block` 四元组扫描内，用 `g_rhf` 的 Coulomb 累加（含 `fab`/`fcd` 的对角 2 倍多重度）产出裸 J，用 `k` 的 8 项交换累加（无 0.5）产出裸 K。**RSH 兼容约束（绑定，来自延后的 rsh-eri-cam-b3lyp）**：K 累加与 J 分离、K 不吸收任何 `a_x` 权重，保持"可再叠加一个不同算子（erf 衰减、omega）的第二个 K"的形状——不得把 J+K 焊死成不透明单例。
- `HamiltonianBuilder::step` 的受限杂化 `Direct`（原 :335–337）分支改为单次 `jk_rhf(P_total_J, J_eng, K)`，随后 `F = Hcore + J_eng − 0.5·alpha_hf·K + Vxc`（`E_K = −0.25·Tr[P·K]`），代数上等于原 `G + 0.5·(1−alpha_hf)·K`。CFMM 分支的 J 仍来自 `cfmm_->j`，故**仅 Direct 分支做 J+K 融合**，CFMM 分支保持调 `k()` 单遍即可，不引入第二遍、不新增开销亦不回退。
- **命名风险（在 spec 显式命名）**：8 倍 ERI 置换多重度对 J 与 K 目标的解包不同（J 用 `fab/fcd` 对角因子，K 用逐项镜像）——这是经典的一遍-JK bug；逐元素 `≤1e-12` 门（ac-003）正是为它设立。

**3. ERI 内层微分配（Item 20）。**
- (a) `ShellBlockPairList::k`（原 :385）与 `g_rhf`（原 :541）里 `double blk[kMaxNCartOrbital⁴] = {}` 每四元组整块清零 ~80KB；改为按实际块大小 `na*nb*nc*nd` 显式循环清零（对齐 `shell_pair.cpp:154–155`，`eri_block` 需预清零后累加）。
- (b) `IntegralPack::coulomb_exchange_fused`（原 :64–77）：J 用线程化 BLAS，K 收缩串行 O(n⁴)；给 K 的外层 `i` 循环加 `#pragma omp parallel for schedule(dynamic)`。每个输出元素 `(i,j≥j')` 及其镜像 `(j,i)` 由单一 `i` 迭代独占、按升序 `k,l` 完整归约，故对线程数**逐位不变**。
- (c) `rys_eri_block_canonical`（原 :416）每四元组 `std::vector<AngIdx> ang_idx(n_out)` 堆分配、以及非规范包装 `rys_eri_block`（原 :1074）的 `std::vector<double> scratch`——均转为 `static thread_local` 复用（沿用同文件 `tls_g_buf`（:380）的 resize-if-smaller 模式）；`ang_idx` 全量覆盖无需清零，`scratch` 使用前 `std::fill` 清零其使用前缀。`ang_idx` 仅依赖 `(la,lb,lc,ld)` 签名、纯整数索引预计算，无浮点，故确定且线程安全。

### Reuse decision

未提供 `librarian_report`（调用方注明蓝图刷新延后），故直接对在树符号裁定：
- `reuse eval_basis` / `eval_basis_val_grad` —— `XCGridBasis::build` 直接调用，不重实现基函数求值。
- `reuse ShellBlockPairList::g_rhf` / `ShellBlockPairList::k` 的累加逻辑 —— `jk_rhf` 复用二者的 Coulomb/交换收缩模式（合并为一遍），不新造 ERI 求值。
- `reuse tls_g_buf` 线程局部模式 —— 新 `tls_ang_idx` / `tls_swap_scratch` 沿用同一 resize-if-smaller、每 worker 独享形状。
- `new — XCGridBasis`：现无持有跨迭代 AO-on-grid 缓存的类型；作为 `HamiltonianBuilder` 私有工作区新增（非公开用户 API，非 `*Config` 包）。
- `new — ShellBlockPairList::jk_rhf`：现无"一遍同出 J 与 K"的成员；`g_rhf` 产出的是 `J−0.5K` 且丢弃 K，无法满足 RSH 可分离 K 的约束，故新增成员，构造/错误处理沿用同类 `g_rhf`/`k` 风格。

## Files to create or modify

- `include/atomiverse/cpu/scf/xc_grid_engine.h` — 新增 `XCGridBasis`、`evaluate_xc_batched` 预算基重载、删除 `active[]`、`if constexpr` 守卫拷贝。
- `include/atomiverse/cpu/scf/internal/shell_pair.h` — 新增 `ShellBlockPairList::jk_rhf` 声明。
- `include/atomiverse/cpu/scf/internal/hamiltonian_builder.h` — 新增 `mutable XCGridBasis xc_basis_;` 私有成员。
- `src/cpu/scf/operators.cpp` — 实现 `jk_rhf`；`k`/`g_rhf` 的 `blk` 按块大小清零；`coulomb_exchange_fused` 的 K 加 OpenMP 归约。
- `src/cpu/scf/internal/hamiltonian_builder.cpp` — XC 基缓存生命周期（`prepare` 失效、首迭代 `build`、复用重载）+ 受限杂化 Direct 分支融合 `jk_rhf`。
- `src/cpu/scf/rys_eri.cpp` — `ang_idx` 与非规范 `scratch` 转 `thread_local` 复用。
- `tests/cpu/scf/test_xc_grid_engine.cpp` (new) — XC 基缓存等价性单测。
- `tests/cpu/scf/test_operators.cpp` (new) — `jk_rhf` 融合等价 + `coulomb_exchange_fused` 线程不变单测。
- `tests/cpu/scf/test_rys_eri.cpp` — 扩充 `thread_local` 复用确定性单测。
- `tests/CMakeLists.txt` — 注册两个新测试目标（`atv_add_cpp_test`）。
- `regressions/scf/b3lyp-fused-jk-h2o/reg_b3lyp_fused_jk_h2o.cpp` (new) — 受限 B3LYP/def2-SVP H2O 单点走融合路径。
- `regressions/scf/b3lyp-fused-jk-h2o/REFERENCE.md` (new) — PySCF 硬编码金标（工具/版本/命令/日期）。
- `regressions/CMakeLists.txt` — 注册新回归用例。

## Tasks

- [x] Write failing unit tests for the geometry-fixed XC-basis cache (tests/cpu/scf/test_xc_grid_engine.cpp → TestXcGridEngine): cached vs rebuild Vxc_alpha/Vxc_beta/E_xc bitwise-equal for LDA and PBE on a stored H2O fixture
- [x] Implement XCGridBasis + evaluate_xc_batched precomputed-basis overload in include/atomiverse/cpu/scf/xc_grid_engine.h; delete dead active[]; if constexpr-guard the P/Vxc double copies
- [x] Implement the per-geometry mutable xc_basis_ cache in HamiltonianBuilder (include/atomiverse/cpu/scf/internal/hamiltonian_builder.h + src/cpu/scf/internal/hamiltonian_builder.cpp): invalidate in prepare(), lazy build on first KS step, reuse across iterations
- [x] Write failing unit tests for fused jk_rhf and thread-invariant stored-ERI K (tests/cpu/scf/test_operators.cpp → TestShellBlockPairList / TestIntegralPack): jk_rhf J,K == g_rhf/k, per-element ‖F_fused − F_twopass‖∞ ≤ 1e-12; coulomb_exchange_fused output identical for 1 vs N threads; and thread_local ERI-reuse determinism (tests/cpu/scf/test_rys_eri.cpp)
- [x] Implement ShellBlockPairList::jk_rhf one-sweep J+K in include/atomiverse/cpu/scf/internal/shell_pair.h + src/cpu/scf/operators.cpp, keeping K separable (unscaled) for RSH extensibility
- [x] Wire fused jk_rhf into the restricted-hybrid Direct branch of HamiltonianBuilder::step in src/cpu/scf/internal/hamiltonian_builder.cpp (F = Hcore + J − 0.5·alpha_hf·K + Vxc), leaving the CFMM branch single-sweep
- [x] Optimize ERI inner-loop allocations in src/cpu/scf/operators.cpp: zero only na·nb·nc·nd of blk[] in ShellBlockPairList::k and g_rhf; add OpenMP reduction to coulomb_exchange_fused K contraction
- [x] Convert per-quartet ang_idx and the non-canonical scratch to thread_local reuse in src/cpu/scf/rys_eri.cpp (tls_g_buf resize-if-smaller pattern; fill scratch prefix before accumulation)
- [x] Add regression example regressions/scf/b3lyp-fused-jk-h2o/ (restricted B3LYP/def2-SVP H2O single point via the fused path; hard-coded PySCF golden energy, no third-party runtime) and register it in regressions/CMakeLists.txt
- [x] Run full check + test suite (ATV_BUILD_HF=ON): full unit ctest + the DFT/SCF regression families reproduce cited values

## Testing strategy

单测仅在 `tests/`（镜像布局、单函数/方法）：

- `tests/cpu/scf/test_xc_grid_engine.cpp` → `TestXcGridEngine`：在一个小型存储网格 fixture（H2O，少量壳）上，断言"先 `XCGridBasis::build` 再走预算基重载"得到的 `Vxc_alpha`、`Vxc_beta`、`E_xc`、`n_*_integrated` 与旧的原地重建路径**逐位相等**（LDA 与 PBE 两分支）；grep 断言 `active` 死数组已删除。
- `tests/cpu/scf/test_operators.cpp` → `TestShellBlockPairList`：在存储的 H2O/def2-SVP（含 d 壳）shell-block fixture 上，`jk_rhf` 产出的 J、K 与 `g_rhf`（换算 `J=G+0.5K`）、`k` 分别相等；对 `a_x ∈ {0.20, 0.25}` 组装 `F_fused` 与两遍 `F_twopass`，硬编码门 `max_abs(F_fused − F_twopass) ≤ 1e-12` Ha（域不变量，守 8 倍多重度解包 bug）。→ `TestIntegralPack`：`coulomb_exchange_fused` 在 `omp_set_num_threads(1)` 与 `=N` 下输出逐位相同。
- `tests/cpu/scf/test_rys_eri.cpp`：扩充——同一 `(la,lb,lc,ld)` 块四元组在单线程与多线程（`#pragma omp parallel`）下 `eri_block` 值不变，验证 `thread_local ang_idx`/`scratch` 复用无跨线程串扰；grep 断言规范核内不再有局部 `std::vector<AngIdx>`。

回归（`regressions/`，公开 API、引用文献/硬编码金标、≤60s）：
- 本 spec 回归示例：`regressions/scf/b3lyp-fused-jk-h2o/` —— 用 `cpu::scf::RKS`（B3LYP、def2-SVP、H2O）单点，经受限杂化融合 `jk_rhf` 路径，复现 **硬编码 PySCF B3LYP/def2-SVP H2O 总能量**（`REFERENCE.md` 注明 PySCF 版本、命令、日期；不在测试期调用/子进程任何第三方）。
- 既有 DFT/SCF 回归族（数值不变守卫，按 `regressions/scf/` 枚举）：`h2_full_grid_alda`（LSDA/ALDA）、`blyp-water-sto3g`（BLYP GGA）、`full-lr-tddft-h2-pbe` 与 `uks-full-lr-o2-lsda`、`uks-tda-o2-lsda`（PBE/LSDA + TDA）、`h2-b3lyp-hybrid-tda`、`h2-b3lyp-rpa-tddft`、`b2plyp-h2o`（杂化/双杂化交换）、`d3bj-molecular-pbe`、`dft-d4-h2o-pbe`（PBE+色散）、`rks-hessian-h2o-sto3g`、`uks-hessian-o2-sto3g`（梯度/Hessian）、`sf-tddft-o2`、`spin-flip-tda-o2`（SF-TDA）——全部复现既有引用值不变。

## Out of scope

- CUDA/torch 后端（本 spec 仅 CPU `cpu/scf`；CUDA SCF 为已知 debt，且 dev 上不编译）。
- 把 `active[]` 接入实际网格筛选（grid screening）——本 spec 仅删死写；筛选优化另立。
- 将 XC 基缓存下推到 `tdhf.cpp` / gradient 的迭代路径——二者保留旧薄包装签名、逐位不变，本 spec 不改其内部循环。
- RSH（CAM-B3LYP、omega）的第二个 erf 衰减 K 求值——仅**保留可扩展形状**，不实现（属延后的 rsh-eri-cam-b3lyp）。
- deriv 核（`rys_eri_block_deriv_canonical` 的 `ang_idx`/`scratch`，:742/:1025）的 thread_local 化——本 spec 仅按 scope 处理 :416/:1074 的能量路径；deriv 孪生留作后续（明确记录、非默默遗留）。
- 任何 `*Config` 包 / 工厂 / 网格阈值配置项（内存策略为无条件缓存，逃生口仅在注释中说明，不落代码）。
