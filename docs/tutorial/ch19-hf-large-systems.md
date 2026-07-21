# 第 19 章　大分子与库仑瓶颈

Hartree–Fock 的 Fock 矩阵只有 \(O(N^2)\) 个元素，但构造它所需的双电子
积分 \((\mu\nu|\lambda\sigma)\) 在形式上是 \(O(N^4)\) 的。大分子计算的
核心困难在此。本章按“如何降低库仑代价”展开：直接积分与筛选、密度拟合、
周期体系的 Ewald/PME、开边界多极方法，以及连续高斯分布上的快速多极
（CFMM）。

---

## 19.1 四中心积分与直接 SCF

把四中心张量存入内存通常不可行，因此采用 **direct SCF**：每步迭代重算
所需积分。Cauchy–Schwarz 上界

\[
|(\mu\nu|\lambda\sigma)|
  \le \sqrt{(\mu\nu|\mu\nu)}\,\sqrt{(\lambda\sigma|\lambda\sigma)}
\]

允许在壳对层面丢弃指数小的贡献，把实际代价压到大约 \(O(N^2)\)（固定
密度、固定精度时）。本引擎的壳块 Schwarz 筛选与 Rys 求积即此路径。

筛选不能改变渐近阶的“常数因子”：远距离壳对仍可能被反复触及。进一步
加速需要换一种表示——用低秩近似压缩库仑，或把远场交给多极/网格方法。

---

## 19.2 密度拟合

乘积密度 \(\rho_{\mu\nu}(\mathbf r)=\chi_\mu(\mathbf r)\chi_\nu(\mathbf r)\)
在空间上高度线性相关。用一组 \(O(N)\) 规模的辅助基 \(\{P\}\) 展开，并在
**库仑度规**下拟合（Vahtras–Almlöf–Feyereisen），得到

\[
(\mu\nu|\lambda\sigma)
  \approx \sum_{PQ}
  (\mu\nu|P)\,[V^{-1}]_{PQ}\,(Q|\lambda\sigma),
\qquad
V_{PQ}=(P|Q).
\]

四中心张量不再出现；存储的是三中心量 \((\mu\nu|P)\)。库仑度规使能量误差
对拟合残差为二阶小量。常用辅助基分为只拟合 \(J\) 的 JFIT 与同时拟合
\(J\) 与 \(K\) 的 JKFIT；交换对辅助基质量更敏感。本引擎仅提供
**def2 通用 JKFIT** 表（含 \(g\) 壳），`RI::J` 与 `RI::JK` 共用该表。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2TZVP;
rhf.ri = cpu::scf::RI::JK; // None | J | JK
rhf.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;
rhf.run(frame);
```

`RI::None` 为精确四中心路径。未接线的组合（如部分开壳、解析力与 RI 同开）
会明确拒绝，而不是静默降级。

水分子 def2-TZVP 上，RI-JK 与精确 ERI 的总能量差约 \(5\,\mu\mathrm{Ha}\)
（约 \(2\,\mu\mathrm{Ha}/\)原子），与 PySCF `density_fit()` 在同一拟合
空间内对齐到 \(10^{-12}\,\mathrm{Ha}\) 量级。拟合误差是系统偏移，对
反应能差往往部分抵消。辅助基须按球谐约定计数，否则与标准 def2-JKFIT
定义不一致。

---

## 19.3 周期边界：Ewald 与 PME

周期库仑和是条件收敛的：裸求和 \(\sum q_i q_j/r_{ij}\) 依赖于求和顺序。
Ewald 将 \(1/r\) 拆成短程与长程，

\[
\frac1r
  = \frac{\mathrm{erfc}(\alpha r)}{r}
  + \frac{\mathrm{erf}(\alpha r)}{r},
\]

实空间用截断的 erfc 项（与邻居列表配合），倒空间解平滑电荷分布的
Poisson 方程。两侧必须使用**同一个** \(\alpha\)。

直接 Ewald 对每个 \(\mathbf k\) 累加结构因子，固定精度下代价随体系增大
接近 \(O(N^2)\)。光滑粒子网格 Ewald（SPME；Essmann 等，1995）把电荷用
B 样条铺到规则网格，以三维 FFT 完成卷积，再把力插回粒子，总标度
\(O(N\log N)\)。

| 方法 | 标度 | 适用 |
|---|---|---|
| 直接 Ewald | \(\sim O(N^2)\) | 小盒、对照、无 FFT |
| SPME | \(O(N\log N)\) | 常规周期 MD |
| FMM / CFMM | \(\sim O(N)\) | 开边界分子库仑（下两节） |

实现分为两层：纯数学在 `include/atomiverse/longrange/`（无 Frame、无
MPI、无 FFT）；CPU 上的 `KSpaceEwald` / `KSpacePME` 负责网格、通信与
FFTW。开关名仍为 `-DATV_BUILD_KSPACE`。电子结构侧若需长程算子，应只
依赖 `atomiverse/longrange/*`。

```cpp
const double alpha = 0.3;  // 1/Å
const double cutoff = 8.0; // Å
sim.add_potential(cpu::CoulLong{alpha, cutoff});
sim.add(cpu::KSpacePME{alpha, /*grid_spacing=*/1.0, /*order=*/6});
```

岩盐结构的 Madelung 常数 \(M\approx 1.7476\)（每离子对）可作为倒空间
实现的文献检验（`tests/cpu/longrange/` 与相关回归）。

---

## 19.4 开边界：快速多极方法

无周期像时，远场可用多极级数压缩。一电荷团在充分分离的观察点处，

\[
\Phi(\mathbf r)
  = \sum_{n=0}^{p}\sum_{m=-n}^{n}
    M_n^m\, I_n^m(\mathbf r-\mathbf c)
  + O\!\Big(\frac1{c-1}\Big(\frac1c\Big)^p\Big),
\]

截断阶 \(p\) 的误差随间隔比 \(c\) 几何下降（Greengard–Rokhlin）。自适应
八叉树上的标准流程为：叶上 P2M，向上 M2M，对充分分离的胞对做 M2L，
向下 L2L，叶上 L2P，相邻叶用直接求和补近场。M2M 与 L2L 在阶 \(p\) 内
精确；截断误差仅来自 M2L。

这与 PME 互补：PME 处理**周期**长程，FMM 处理**开边界**远场。二者边界
条件不同，不能互换。

```cpp
#include "atomiverse/longrange/fmm_tree.h"

ATV::longrange::FmmTree tree(/*order=*/8, /*max_per_leaf=*/8);
double e = tree.energy(positions, charges); // 无量纲库仑和
```

对随机与晶格点电荷集，相对直接 \(O(N^2)\) 求和的能量误差随 \(p\) 下降；
\(p=16\) 时相对误差可低于 \(10^{-10}\)（见
`regressions/longrange/fmm_point_charges`）。将 \(N\) 增至 \(4N\) 时，
墙钟比应明显低于二次标度的期望值 16。

---

## 19.5 连续快速多极与 HF 库仑矩阵

分子积分中的库仑矩阵 \(J\) 来自高斯乘积分布，不是点电荷。White、
Johnson、Gill 与 Head-Gordon（1994）给每个乘积分布一个延伸半径

\[
r_{\mathrm{ext}}
  = \sqrt{2/\zeta}\,\mathrm{erfc}^{-1}(\varepsilon_{\mathrm{thr}}),
\qquad
\zeta=\zeta_a+\zeta_b,
\]

仅当两分布中心距大于延伸半径之和时才允许多极远场，否则走解析 ERI
近场。该判据在 `longrange/gaussian_extent.h`。

在 HF 中，CFMM 作为 **\(J\) 的第三条路径**（精确四中心、密度拟合、CFMM）
接入；**交换 \(K\) 仍走精确或 RI 路径，不经多极。**

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.coulomb = cpu::scf::Coulomb::CFMM;
rhf.cfmm_order = 10;
rhf.cfmm_threshold = 1e-9;
rhf.run(frame);
```

默认 `Coulomb::Direct` 与历史精确路径一致。`Coulomb::CFMM` 与
`ri = RI::J` 或 `RI::JK` 同时设置会被拒绝。

正确性标准必须先于任何速度比较：CFMM 的 \(J\) 应在给定容差内复现
Schwarz 直接 \(J\)。PySCF 分子代码不提供 FMM/CFMM，因此外部对照只能
钉直接路径；CFMM 的 oracle 是本库的 Direct。水/STO-3G 上单次 \(J\) 的
元素最大偏差 \(\le 10^{-8}\,\mathrm{Ha}\)；正癸烷 C₁₀H₂₂/STO-3G 上
收敛总能偏差 \(\le 10^{-7}\,\mathrm{Ha}\)。

`cfmm_order` 控制多极截断；`cfmm_threshold` 控制延伸半径与近远场划分。
小分子在默认阈值下可能几乎全部落在近场，此时 CFMM 与 Direct 数值重合，
墙钟优势要到更大体系才出现。标度测量方法见
`benchmarks/scf/validation/cfmm_scaling/`。

---

## 数值对照摘要

| 体系 | 量 | 说明 |
|---|---|---|
| H₂O def2-TZVP | RI-JK 与精确 ERI | \(\Delta E\sim 5\,\mu\mathrm{Ha}\) |
| 岩盐 | Madelung | \(\approx 1.7476\) / 离子对 |
| 点电荷 FMM | 相对直接和 | \(p=16\) 时 \(\lesssim 10^{-10}\) |
| C₁₀H₂₂ STO-3G | CFMM vs Direct 总能 | \(\le 10^{-7}\,\mathrm{Ha}\) |

更完整的 PySCF 能量对照见章内回归与 [附录](appendix-pyscf-oracle.md)。

### 进一步阅读

- Jensen, *Introduction to Computational Chemistry*，SCF 瓶颈与线性标度。
- Weigend, *Phys. Chem. Chem. Phys.* **8**, 1057 (2006)（辅助基）。
- Essmann 等, *J. Chem. Phys.* **103**, 8577 (1995)（SPME）。
- Greengard & Rokhlin, *J. Comput. Phys.* **73**, 325 (1987)（FMM）。
- White 等, *Chem. Phys. Lett.* **230**, 8 (1994)（CFMM）。
