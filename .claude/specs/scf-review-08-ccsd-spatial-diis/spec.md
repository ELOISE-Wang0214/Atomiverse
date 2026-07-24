---
slug: scf-review-08-ccsd-spatial-diis
title: "Closed-shell spatial CCSD + CC amplitude DIIS"
status: done
created: 2026-07-24
domain: scf / coupled-cluster
queue: scf-review
depends_on: [scf-review-02-eri-hoists, scf-review-03-posthf-dedup, scf-review-07-posthf-perf]
grilled: true
---

# Closed-shell spatial CCSD + CC amplitude DIIS

## Summary

针对 review 条目 19（当前最大的单项性能工作）：`src/cpu/scf/ccsd.cpp` 的闭壳层 CCSD 完全在**自旋轨道基**（`no=2·n_occ`、`nv=2·n_vir`，见 :463-464）中求解，FLOP 与内存约为闭壳层空间轨道公式的 16 倍，约一半的收缩在与保证为零的自旋块相乘（`so_eri`），且每个中间量都是手写的标量四重循环（限速的 `Wabef·τ` 在 :655-674，`vvvv` 全张量被物化两次于 :489/:655，`ring_at` 每个振幅调用四次于 :735-745）；收敛用纯 Jacobi（:913-914），无 DIIS，迭代数多 2-3 倍。本 spec 交付两件事、**严格按此顺序**：先做低风险、立即惠及现有自旋轨道代码的 **CC 振幅 DIIS**，再做闭壳层 **空间轨道 CCSD 重写**。用户可见的结果是：`RHF::ccsd(frame)` 在闭壳层参考下自动走空间轨道快路径，收敛更快、内存更省，而 EOM/(T)/梯度/CC2/UCCSD 的数值结果逐位不变。

## Domain basis

- **自旋轨道 CCSD（保留为验证 oracle）**：Stanton, Gauss, Watts, Bartlett, J. Chem. Phys. **94**, 4334 (1991)，DOI 10.1063/1.460620；Gauss & Stanton, J. Chem. Phys. **103**, 3561 (1995)。现有 `so_ccsd_kernel` 即此公式，须保持可调用。
- **闭壳层空间轨道 CCSD（本 spec 新增默认路径）**：Scuseria, Janssen, Schaefer, J. Chem. Phys. **89**, 7382 (1988)，DOI 10.1063/1.455269。空间振幅 `t_i^a`（T1）、`t_ij^ab`（T2），存储量 `n_o²·n_v²`。库仑-交换组合 `L_iajb = 2·(ia|jb) − (ib|ja)`（化学家记号）；相关能
  `E_corr = Σ_ijab L_iajb · ( t_ij^ab + t_i^a·t_j^b )`。
  自旋求和后的中间量以融合复合下标（occ·vir）的 GEMM 形式经 `xt::linalg` 收缩。
- **CC 振幅 DIIS**：Scuseria, Lee, Schaefer, Chem. Phys. Lett. **130**, 236 (1986)，DOI 10.1016/0009-2614(86)80461-4；Pulay, Chem. Phys. Lett. **73**, 393 (1980)。以扁平化向量 `T = [t1; t2]` 为振幅、误差向量 `e = Δt = R/D`（即 Jacobi 增量）。增广正规方程
  `[[B, −1], [−1ᵀ, 0]]·[c; λ] = [0; −1]`，其中 `B_ij = e_i·e_j`；外推 `T = Σ_i c_i·T_i`。`B` 近奇异时正则化（Tikhonov）或丢弃最旧向量。
- **不变量（永久回归门，非一次性）**：
  1. 空间轨道 `E_corr` 与自旋轨道 `E_corr` 在 H2O/STO-3G 尺度上相差 ≤ 1e-9 Ha，两者均收敛至 `max|Δt| < 1e-8`。
  2. 相同残差容差下，`E_CCSD`（含 DIIS）与 `E_CCSD`（纯 Jacobi）相差 ≤ 1e-9 Ha（同一定点）。
- **单位**：CC 相关能与总能均为 Hartree（公开面按既有约定，`frame.meta.get(kEnergy)` 仍走 SCF 侧）。迭代数收益是 impl summary 中的实测证据，**不在数值上编码约 2-3× 的加速倍数**。

## Design

**风险命名（须显式验证）**：空间轨道重写的经典 bug 是自旋求和系数错误（`2J−K` 权重、exchange/ring 项里 `τ` 与 `t2` 混淆）。oracle 交叉校验回归正是为此存在。

**入口与兼容性（明确的选择点）**：
- `ccsd_energy(ints, det, e_scf_ha, tol, max_iter)`（即公开的 `MeanFieldMethod::ccsd → RHF::ccsd`）在闭壳层 RHF 参考下**默认路由到新的 `ccsd_energy_spatial`**（空间轨道快路径）。这是唯一切换空间/自旋轨道的选择点。
- `ccsd_solve(...)`（返回带自旋轨道块的 `CCSDSolution`）**保持自旋轨道、签名与语义完全不变**：EOM-EE/EA/IP-CCSD 的 `build_imds` 依赖 `CCSDSolution` 的 SO 块，`triples_energy`/`ccsd_t_energy`（(T)）依赖 SO 振幅，`ccsd_gradient`（`internal::ccsd_analytic_forces`）与 `cc2_energy` 亦然——它们继续原样消费 `CCSDSolution`。
- `uccsd_solve(...)`（UHF/开壳层 UCCSD）**保持走 `so_ccsd_kernel` 自旋轨道路径不变**，因此开壳层与 EOM-SF 不受影响。
- 结论：空间轨道成为**闭壳层能量入口的默认**，自旋轨道保留给开壳层 + oracle + 所有 `CCSDSolution` 消费者。

**新符号与所有权（默认 OOP 形态）**：
- `CcAmplitudeDiis`（`internal/cc_diis.h`）——CC 局部 DIIS 加速器，OOP 类。构造参数 `int max_vectors = 8`（构造参数直接赋成员，**无 `*Config` 包、无工厂函数**）。主方法 `std::vector<double> step(const std::vector<double>& amp, const std::vector<double>& err)`：存入一对 `(amp, err)`，解增广正规方程并返回外推振幅；辅以 `void reset()`、`int size() const`。内部 `std::deque<std::vector<double>>` 存历史。命名、构造、错误处理刻意贴合 `internal::DIISMixer`（`max_vectors` 构造参、`step`/`reset`/`size`、deque 历史），使新代码读起来与既有 SCF DIIS 一致。
- `ccsd_energy_spatial(const IntegralPack&, const SlaterDeterminant&, double e_scf_ha, double tol, int max_iter)`（`internal/ccsd_spatial.h`，返回 `CCSDResult`）——闭壳层空间轨道核；复用既有 `ao2mo_full` 的空间 MO ERI（不做自旋展开，直接切片 spatial 块）与 `CcAmplitudeDiis`。这是计算函数（与既有 `ccsd_energy`/`ccsd_solve` 同类），非工厂。

**Reuse decision**（解决 librarian advisory 中的每个候选）：
- `internal::DIISMixer`（`internal/diis_mixer.h`）→ **new — DIISMixer 面向 Fock（自旋分辨的 `RealMatrix` 对、`deque<RealMatrix>` 二维历史、`Tr(eᵢ·eⱼ)` 度量、FPS−SPF 伴随），不适用于扁平 `[t1;t2]` 向量；且泛化它会触及活跃的 SCF mix 阶段**。新 `CcAmplitudeDiis` 仅复用 Pulay B-矩阵最小二乘的形状，并沿用其命名/构造风格。**不得** generalize，不得触碰 `diis_mixer.*` 或 SCF mix 阶段。
- `ao2mo_full`（既有 AO→MO 变换）→ **reuse**：空间路径直接复用同一 MO ERI 变换，不重实现。
- `so_ccsd_kernel`（既有自旋轨道核）→ **reuse/keep**：作为 oracle 与开壳层路径保留，仅在其振幅循环内接入 DIIS。

## Files to create or modify

- `include/atomiverse/cpu/scf/internal/cc_diis.h` (new) — `CcAmplitudeDiis` 声明。
- `src/cpu/scf/internal/cc_diis.cpp` (new) — B-矩阵增广正规方程求解（`xtensor-blas`）。
- `include/atomiverse/cpu/scf/internal/ccsd_spatial.h` (new) — `ccsd_energy_spatial` 声明。
- `src/cpu/scf/ccsd_spatial.cpp` (new) — 闭壳层空间轨道 CCSD 核实现。
- `src/cpu/scf/ccsd.cpp` (modify) — 将 `CcAmplitudeDiis` 接入 `so_ccsd_kernel` 振幅循环；`ccsd_energy` 闭壳层路由到 `ccsd_energy_spatial`；`ccsd_solve`/`uccsd_solve` 不变。
- `CMakeLists.txt` (modify) — 在 `ATV_BUILD_HF` 块内加入 `src/cpu/scf/internal/cc_diis.cpp` 与 `src/cpu/scf/ccsd_spatial.cpp`。
- `tests/cpu/scf/test_cc_diis.cpp` (new) — `CcAmplitudeDiis` 单元测试。
- `tests/cpu/scf/test_ccsd_spatial.cpp` (new) — `ccsd_energy_spatial` 单元测试。
- `tests/CMakeLists.txt` (modify) — 在 `ATV_BUILD_HF` 块内注册 `test_cpu_scf_cc_diis`、`test_cpu_scf_ccsd_spatial`。
- `regressions/scf/ccsd-spatial-so-parity/reg_ccsd-spatial-so-parity.cpp` (new) — 空间 vs 自旋轨道、DIIS vs Jacobi 双不变量 + PySCF golden。
- `regressions/scf/ccsd-spatial-so-parity/REFERENCE.md` (new) — 引用值出处（PySCF 版本/命令/日期）。
- `regressions/CMakeLists.txt` (modify) — 在 `ATV_BUILD_HF` 块注册 `atv_add_regression(scf ccsd-spatial-so-parity)`。

## Tasks

- [x] Write failing unit tests for CcAmplitudeDiis in tests/cpu/scf/test_cc_diis.cpp and register test_cpu_scf_cc_diis in tests/CMakeLists.txt
- [x] Implement CcAmplitudeDiis in include/atomiverse/cpu/scf/internal/cc_diis.h + src/cpu/scf/internal/cc_diis.cpp and register the source in CMakeLists.txt (ATV_BUILD_HF)
- [x] Wire CcAmplitudeDiis into the spin-orbital CCSD amplitude loop in src/cpu/scf/ccsd.cpp (so_ccsd_kernel; error e = Δt = R/D)
- [x] Write failing unit tests for ccsd_energy_spatial in tests/cpu/scf/test_ccsd_spatial.cpp and register test_cpu_scf_ccsd_spatial in tests/CMakeLists.txt
- [x] Implement closed-shell spatial CCSD in include/atomiverse/cpu/scf/internal/ccsd_spatial.h + src/cpu/scf/ccsd_spatial.cpp (reuse ao2mo_full + CcAmplitudeDiis) and register the source in CMakeLists.txt (ATV_BUILD_HF)
- [x] Route closed-shell ccsd_energy to ccsd_energy_spatial in src/cpu/scf/ccsd.cpp, keeping ccsd_solve/uccsd_solve spin-orbital and unchanged for EOM/(T)/gradient/CC2/UCCSD
- [x] Add Doxygen @brief/@param(units)/@return to CcAmplitudeDiis and ccsd_energy_spatial
- [x] Add regression example regressions/scf/ccsd-spatial-so-parity (spatial-vs-spin-orbital E_corr and DIIS-vs-Jacobi fixed point on H2O/STO-3G; hard-coded PySCF goldens, no third-party runtime) and register it in regressions/CMakeLists.txt
- [x] Verify existing CCSD/EOM/CCSD-gradient/CC2 regressions and SCF unit tests are unchanged and green
- [x] Run full check + test suite

## Testing strategy

- **Unit — `tests/cpu/scf/test_cc_diis.cpp` → `CcAmplitudeDiis` suite** (targets `CcAmplitudeDiis::step` in isolation, pure algebra, no SCF):
  - happy path: 在一个已知的线性混合玩具序列上，DIIS 外推使误差范数单调下降，几步后残差 ≤ 硬编码阈值（hard-coded 期望向量）。
  - edge: 单个历史向量时 `step` 返回输入本身；近奇异 `B`（近平行误差向量）触发丢弃最旧/正则化，返回有限值、无 NaN。
  - state: `reset()` 清空历史，`size()` 反映当前历史长度（硬编码计数）。
- **Unit — `tests/cpu/scf/test_ccsd_spatial.cpp` → `CcsdSpatial` suite** (targets `ccsd_energy_spatial`；沿用 `test_mp2.cpp` 的 He/H2 STO-3G 亚毫秒参考惯例):
  - domain validation: H2/STO-3G 闭壳层 `E_corr = −0.020524527092106 Ha`（与 `regressions/scf/ccsd-h2-sto3g` 同一 PySCF pin，硬编码），`converged == true`。
  - edge: He/STO-3G（零虚轨道）`E_corr == 0`。
- **Regression（永久，非一次性）— `regressions/scf/ccsd-spatial-so-parity/`**：库外消费者（仅公开头 + 链接 atomiverse），单个 `main()` 在 H2O/STO-3G 上：
  1. 跑闭壳层空间轨道路径（含 DIIS）得 `E_corr_spatial`；
  2. 跑自旋轨道 oracle 路径得 `E_corr_so`，两者均收敛至 `max|Δt| < 1e-8`，断言 `|E_corr_spatial − E_corr_so| ≤ 1e-9 Ha`；
  3. 跑自旋轨道 + 关闭 DIIS（纯 Jacobi）得 `E_corr_jacobi`，断言 `|E_corr_so − E_corr_jacobi| ≤ 1e-9 Ha`（相同残差容差、同一定点）；
  4. 断言 `E_corr_spatial` 与嵌入的 PySCF `cc.CCSD` H2O/STO-3G golden 相符 ≤ 1e-6 Ha（literal，注释注明 tool/version/command/date；ctest 时不调用任何第三方）。
- **Perf 证据（advisory）**：impl summary 记录 `RHF::ccsd` 在回归案例上的 before/after 墙钟时间与含/不含 DIIS 的 CCSD 迭代数；不编码加速倍数、不设数值门。

## Out of scope

- 开壳层/UHF 空间轨道（UCCSD）重写——保留自旋轨道 `uccsd_solve`。
- CCSD(T) 空间轨道化——(T) 继续消费自旋轨道 `CCSDSolution`（`triples_energy`/`ccsd_t_energy` 不变）。
- EOM-EE/EA/IP/SF 空间轨道中间量——EOM 继续消费自旋轨道 `CCSDSolution`（`build_imds` 不变）。
- CCSD 解析梯度空间轨道化——`ccsd_gradient` 继续走自旋轨道路径。
- CC2 空间轨道化——`cc2_energy` 不变。
- CUDA CCSD。
- 泛化或修改 `internal::diis_mixer.*` 与 SCF mix 阶段（librarian 明确禁止）。
- （已考虑替代：可让 `ccsd_solve` 也切空间轨道并重建 SO 块供 EOM——被否，因会破坏 EOM/(T)/梯度的 SO 契约且收益仅在能量项，风险远高。）
