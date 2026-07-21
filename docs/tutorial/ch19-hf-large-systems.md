# Chapter 19 — 大体系 HF:从 Schwarz 筛选到线性标度库仑

> **状态:部分成文**(19.1–19.4 已填充;19.5 仍为骨架)。各节由标注的链环填充;
> 本章对应 Jensen《Introduction to Computational Chemistry》§3.8 后半
> ("the SCF bottleneck" / linear scaling):方法为主干,把 direct SCF、密度
> 拟合与多极方法编成一条"库仑瓶颈如何被逐级击破"的叙事。

## 19.1 SCF 的瓶颈:四指标 ERI

Jensen(§3.8)把大体系 HF 的困境总结为一句话:Fock 矩阵本身只有
O(N²) 个元素,但喂给它的双电子积分张量 (μν|λσ) 有 **O(N⁴)** 个。存不
下,就每步重算(direct SCF);重算太贵,就用 Cauchy–Schwarz 上界
√(μν|μν)·√(λσ|λσ) 在壳对层面把指数衰减的贡献整批筛掉——这正是
Atomiverse 的现役路径(`ShellBlockPairList::g_rhf`,壳块级 Schwarz 筛
选 + Rys 求积),也是 PySCF 的分子路线(诚实注明:PySCF 分子代码同样
没有 FMM,它的快来自积分内核与迭代数,不是更好的标度)。

Schwarz 筛选把渐进标度压到 ~O(N²)(远壳对指数消失),但**前置系数**
仍由每步全量重算的四指标积分决定。本章逐级击破这个库仑瓶颈:先用密度
拟合把 4-指标压成 3-指标(本节),再把远场交给多极展开(19.3–19.5)。

## 19.2 密度拟合 RI

### 4-指标 → 3-指标:为什么能压缩

四指标张量贵,是因为它把**乘积密度** ρ_μν(r) = χ_μ(r)χ_ν(r) 的所有
O(N²) 个成员两两配对。密度拟合(RI,resolution of the identity)的观
察是:这 O(N²) 个乘积密度挤在同一段分子空间里,线性相关得厉害——用一
组 O(N) 的**辅助基** {P} 就能把它们展开:

```
ρ_μν(r) ≈ Σ_P c^P_μν P(r)
```

在**库仑度规**下求拟合系数(Vahtras–Almlöf–Feyereisen 1993,
DOI 10.1016/0009-2614(93)89151-7),得到经典的三中心分解:

```
(μν|λσ) ≈ Σ_PQ (μν|P) [V⁻¹]_PQ (Q|λσ),   V_PQ = (P|Q)
```

四指标张量从此不再出现:存的是 O(N²·N_aux) 的 (μν|P)。库仑度规的妙处
在于它让 J 的能量误差是拟合残差的**二阶小量**——残差 δρ 对能量的一阶
项 (δρ|ρ) 恰好被拟合条件消掉,剩下 (δρ|δρ)。配对好的辅助基下,这是每
原子几十 μHa 的量级(Weigend 2006,DOI 10.1039/b515623h)。

工程实现(`internal::DensityFitting`):V = L·Lᵀ 做 Cholesky,取白化
因子 W = L⁻¹(W·V·Wᵀ = I),把 B_{μν,P} = Σ_Q (μν|Q)·W_QP 一次算好。
之后每步 SCF:

- **RI-J(γ_P 路线)**:γ = Bᵀ·vec(P),J = B·γ——两次 GEMV;
- **RI-K(B 半变换路线,Weigend 2002,DOI 10.1039/B204199P)**:
  B^P_μi = Σ_ν B^P_μν C_νi,K_μν = Σ_{P,i} B^P_μi B^P_νi。

### 单表取舍:只有 JKFIT,没有 JFIT

文献里辅助基分两族:只拟 J 的 JFIT 和拟 J+K 的 JKFIT。Atomiverse 只
收 **def2 万能 JKFIT** 一张表(`basis_data/def2_jkfit.h`,含 g 壳,
l = 4):交换拟合的误差比库仑大一个量级,**必须** JKFIT——JFIT 不够;
反过来用 JKFIT 兼任纯 J 拟合只是辅助维度略大,精度不降。一表两用,少
一张表少一类失配错误。配套的陷阱也直接编码:轨道基元素在 JKFIT 表无
行时**抛异常**,没有静默回退。

### 动手:一行成员开关

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2TZVP;
rhf.tol   = 1e-12;
rhf.ri    = cpu::scf::RI::JK;                     // None | J | JK
rhf.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;    // 默认值,可省略
rhf.run(frame);          // frame.meta[kEnergy];rhf.energy() [Ha]
```

`RI::None`(默认)与旧路径逐位一致;`RI::J` 只拟库仑(交换仍走精确壳
块 K);`RI::JK` 全拟。未接线的组合(ROHF/GHF/SOSCF/解析力)一律抛异
常说明原因。

### 诚实的误差预算(μHa)

对照 regressions/scf/h2o_ri_jk 与 benchmarks/scf/validation/ri_speed
的实测(def2 万能 JKFIT,球谐辅助跨度):

| 体系 | E(RI-JK) − E(精确 ERI) | 每原子 | 与 PySCF `density_fit()` 之差 |
|---|---|---|---|
| H2O / def2-TZVP | +5.41 μHa | 1.8 μHa | 3.3e-13 Ha(同一拟合空间) |
| 甘氨酸 / def2-SVP | +258 μHa | 26 μHa | ~1e-9 Ha |

两点提醒:(1) 误差随基组配对质量浮动——SVP 配万能 JKFIT 比 TZVP 粗;
(2) 拟合误差是**能量的系统偏移**,不是噪声,对反应能差还会进一步抵
消。注意辅助跨度必须取**球谐**(与 def2 家族的定义一致):笛卡尔 g 壳
会多张出低阶函数,拟合虽仍变分,但数值就对不上 PySCF 的
def2-universal-jkfit 了。

## 19.3 长程库仑 I:Ewald 与 PME

周期体系中的库仑和是**条件收敛**的:裸的 ∑_{i<j} q_i q_j / r_{ij} 在
PBC 下按求和顺序发散(正/负表面偶极层)。Ewald 求和把 1/r 劈成短程
(在实空间用 erfc 截断)与长程(在倒空间用高斯平滑的 Poisson 解):

\[
\frac{1}{r}
  = \frac{\mathrm{erfc}(\alpha r)}{r}
  + \frac{\mathrm{erf}(\alpha r)}{r}
\]

实空间半截 `CoulLong(α, cutoff)` 只走邻居列表;倒空间半截由
`KSpaceEwald` 或 `KSpacePME` 提供。两边**必须共用同一个 α**
(计划在 `build()` 时 fail-closed 校验,见
[`docs/developer/kspace-long-range.md`](../developer/kspace-long-range.md))。

### 为何 PME 是 O(N log N)

直接 Ewald 对每个 k 向量累加结构因子 S(k)=∑_j q_j e^{ik·r_j},
k 点数随精度与体积增长 → 整体 O(N²)(固定精度、密度时对 N)。
**光滑粒子网格 Ewald**(SPME; Essmann *et al.*, *J. Chem. Phys.*
**103**, 8577 (1995), DOI: [10.1063/1.470117](https://doi.org/10.1063/1.470117))
把电荷用 B 样条铺到规则网格,一次 3D FFT 完成所有 k 模的卷积,再
样条插回原子力:

| 步骤 | 复杂度 |
|---|---|
| 电荷铺展(B 样条 stencil) | O(N · n³)(n = 样条阶,常数) |
| 正/逆 FFT | O(K log K), K ∝ N(固定网格密度) |
| 力插值 | O(N · n³) |

于是总标度 O(N log N)。Ewald 1921 原始公式见 *Ann. Phys.* **369**,
253, DOI: [10.1002/andp.19213690304](https://doi.org/10.1002/andp.19213690304)。

### 代码居所(facade + 薄包装)

| 层 | 路径 | 内容 |
|---|---|---|
| 纯算子数学 | `include/atomiverse/longrange/` | `EwaldSum`, `bspline` / `spline_moduli` / `good_fft_size`, 影响函数与力插值 — **无** Frame / Component / MPI / FFT |
| CPU 包装 | `include/atomiverse/cpu/longrange/` | `KSpaceEwald` / `KSpacePME` — Frame span、Allreduce、FFTW 计划、meta 写回 |
| 开关 | `-DATV_BUILD_KSPACE=ON` | 门控 FFTW3 / PME(名字保留) |

规则(为后续 FMM/CFMM 立此存照):**SCF 侧只许 include
`atomiverse/longrange/*`,永不 `cpu/longrange/*`。**

```cpp
// Classical MD: real-space CoulLong + reciprocal KSpacePME
cpu::md::MD sim;
const double alpha = 0.3;   // 1/Å — must match both halves
const double cutoff = 8.0;  // Å
sim.add_potential(cpu::CoulLong{alpha, cutoff});
sim.add(cpu::KSpacePME{alpha, /*grid_spacing=*/1.0, /*order=*/6});
// … integrator, neighbor, build, run …
```

### 动手:NaCl Madelung 常数

岩盐结构(交替 ±1 电荷的简单立方)每离子对的 Madelung 常数

\[
M = 1.74756459\ldots
\quad\Rightarrow\quad
E = -\tfrac12 N\, M\, k_C / d
\]

(`k_C = units::coulomb`, d = 最近邻间距)。回归
`regressions/md/npt_pme_nacl` 与单元
`tests/cpu/longrange/test_{ewald,pme}_coul.cpp` 钉住 M ≈ 1.7476
(相对 1e-4 级);构建需 `-DATV_BUILD_KSPACE=ON`。

### 章末标量表(长程 I)

| 方法 | 标度 | 边界 | 何时用 |
|---|---|---|---|
| 直接 Ewald (`KSpaceEwald`) | O(N²) | PBC | 小盒、对照、无 FFTW |
| SPME (`KSpacePME`) | O(N log N) | PBC | 生产 MD / NPT |
| FMM / CFMM | O(N) | 开边界或 PBC | 下两节(环 20/21) |

## 19.4 长程库仑 II:快速多极方法

周期边界的长程库仑用 Ewald/PME(§19.3)在倒空间卷积;开边界分子体系
没有周期像,远场压缩换一条路——**快速多极方法**(FMM)。Greengard &
Rokhlin(*J. Comput. Phys.* **73**, 325 (1987),
DOI:[10.1016/0021-9991(87)90140-9](https://doi.org/10.1016/0021-9991(87)90140-9))
把一球电荷在远处的势写成多极级数;White & Head-Gordon
(*J. Chem. Phys.* **101**, 6593 (1994),
DOI:[10.1063/1.468354](https://doi.org/10.1063/1.468354))给出实固体谐波
形式的平移算子,使层级八叉树上的游走可以按 O(N) 组装全部远场。

### 多极压缩远场

一点电荷 q 在 |r| 大处的势是 q/|r|。一**团**落在中心 c 附近的电荷,
在 |r−c| 大于团半径时可以写成

\[
\Phi(\mathbf r)
  = \sum_{n=0}^{p}\sum_{m=-n}^{n}
      M_n^m\, I_n^m(\mathbf r-\mathbf c)
      \;+\; O\!\left(\frac{1}{c-1}\Big(\frac{1}{c}\Big)^{p}\right)
\]

其中 \(M_n^m\) 是多极系数(对每个源做 P2M 累加),\(I_n^m\) 是非正则
固体谐波,截断阶 p 的误差按间隔比 c 几何衰减——这就是验收时"误差表
逐阶几何衰减"的来源。对**接收**侧,同一势可以在无源球内再展成局部
(Taylor)级数 \(L_n^m R_n^m\),求值是 L2P。

### M2M / M2L / L2L 游走

自适应八叉树罩住包围立方体,叶容量为成员 `max_per_leaf`:

| 阶段 | 算子 | 作用 |
|---|---|---|
| 上行 | P2M → M2M | 叶上累加多极,向父胞精确平移(阶内 exact) |
| 交互 | M2L | 对 well-separated 胞对(≥1 胞宽间隔)把源多极译成目标局部展开;唯一引入截断误差的一步 |
| 下行 | L2L | 父局部展开精确下传到子胞 |
| 叶层 | L2P + 近场直和 | 远场用局部展开求值;相邻叶用 O(n_leaf²) 直和补齐 |

M2M 与 L2L 在截断阶 p 内是**恒等式**(单元测试钉 1e-12);只有 M2L 带
几何截断。这与 PME 的分工对偶:**PME 用网格 FFT 处理周期长程,FMM 用
层级多极处理开边界远场**——两者都把 O(N²) 的库仑和压到近线性,但边界
条件不同,不能互相替代。

### 代码居所

| 层 | 路径 | 内容 |
|---|---|---|
| 展开 | `longrange/multipole_expansion.h` | `MultipoleExpansion` / `LocalExpansion`(P2M、M2P、L2P) |
| 平移 | `longrange/fmm_operators.h` | `M2M` / `M2L` / `L2L` |
| 树 | `longrange/fmm_tree.h` | `FmmTree{order, max_per_leaf}` → `energy` / `potentials` |
| 核 | 无量纲 Σ q_i q_j / r_ij | e²/Å → eV 的换算留给链环 21(CFMM/HF 接线) |

```cpp
#include "atomiverse/longrange/fmm_tree.h"

ATV::longrange::FmmTree tree(/*order=*/8, /*max_per_leaf=*/8);
// positions: std::vector<ATV::Vec3>, charges: std::vector<ATV::AccReal>
double e = tree.energy(positions, charges);   // 无量纲;× units::coulomb → eV
std::vector<ATV::AccReal> phi;
tree.potentials(positions, charges, phi);     // Φ_i = Σ_{j≠i} q_j / r_ij
```

### 误差-阶数收敛(标量表)

随机与晶格点电荷集对 O(N²) 直和的相对能量误差(回归
`regressions/longrange/fmm_point_charges`,叶容量 8,MAC c≥2;完整表见该目录
`REFERENCE.md`)。每两阶取一点看包络;文档阶 **p★ = 16** 处相对误差 ≤ 1e-10。

| p | 双团簇随机(2×32) rel.err | 双晶格块(2×3³) rel.err |
|---:|---:|---:|
| 2 | 6.4e-6 | 1.0e-8 |
| 4 | 1.2e-7 | 1.6e-9 |
| 8 | 3.4e-10 | 7.4e-11 |
| 12 | 3.3e-13 | 2.3e-12 |
| **16** | **1.7e-15** | **2.0e-14** |

包络误差比 < 1 → 几何衰减,与 Greengard–Rokhlin 界
\((1/(c-1))(1/c)^p\) 同型(间隔比 c≥2 由 well-separated MAC 保证)。

### 标度:近线性(曲线)

同一回归对 N 与 4N 两次求值(相同 p、叶容量)的墙钟比 < 8(O(N²)
期望 16)。下图是示意性的 FMM vs 直和标度(对数坐标;实测点随机器浮动,
断言只钉比值上界)。

```molplot preset="molplot" theme="auto" width="100%" aspect="16:9"
mark: {type: line, point: true, tooltip: true}
data:
  values:
    - {N: 64, series: FMM, t: 1.0}
    - {N: 128, series: FMM, t: 2.1}
    - {N: 256, series: FMM, t: 4.3}
    - {N: 512, series: FMM, t: 8.8}
    - {N: 64, series: direct, t: 1.0}
    - {N: 128, series: direct, t: 4.0}
    - {N: 256, series: direct, t: 16.0}
    - {N: 512, series: direct, t: 64.0}
encoding:
  x: {field: N, type: quantitative, title: "N (particles)", scale: {type: log}}
  y: {field: t, type: quantitative, title: "Relative wall time", scale: {type: log}}
  color: {field: series, type: nominal, title: null}
```

FMM 近似斜率 1(近线性);直和斜率 2。链环 21 把同一算子核接到
高斯分布(CFMM)并写入 HF 的 J 矩阵。

## 19.5 长程库仑 III:CFMM 加速 HF

Jensen 线性标度叙事的收官:分子库仑从

1. **Direct + Schwarz**(精确 ERI,O(N²) 筛选后仍超线性)——**这正是 PySCF 的做法**;PySCF **不含** FMM/CFMM(可选 DF 是另一条路),故 CFMM 正确性的 oracle 是我们自家已被 PySCF 钉定的 Direct J,而非外部 FMM 数;
2. **RI-J / RI-JK**(链环 12,拟合辅助基,中等体系的速度赢家);
3. **CFMM**(本环):连续高斯分布的快速多极,线性标度选项。

### 延伸分类(WJGH 1994)

乘积高斯(总指数 ζ = ζ_a+ζ_b)不是点电荷。赋予延伸

\[
r_\mathrm{ext} = \sqrt{2/\zeta}\,\mathrm{erfc}^{-1}(\varepsilon_\mathrm{thr})
\quad[\mathrm{Bohr}]
\]

两分布仅当 \(|R_A-R_B| > r_\mathrm{ext}(A)+r_\mathrm{ext}(B)\) 才允许走多极,否则归近场解析 ERI。公式与判据在门面头 `longrange/gaussian_extent.h`(header-only)。

### 接入(成员驱动,第三条 J 路)

```cpp
#include "atomiverse/cpu/scf/methods.h"

cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.coulomb = cpu::scf::Coulomb::CFMM; // Direct | RI | CFMM
rhf.cfmm_order = 10;                  // 多极截断阶 p
rhf.cfmm_threshold = 1e-9;            // ε_thr → r_ext
rhf.run(frame);                       // K 永远不走 CFMM
// rhf.coulomb = CFMM 且 rhf.ri = RI::J/JK 会抛异常(J 路二选一)
```

默认 `Coulomb::Direct` 下现路**位级不变**。SCF 只 include 顶层
`atomiverse/longrange/*`,永不 include `cpu/longrange/`。

### 精度表(硬门:先于任何速度论断)

| 体系 | 量 | 容差 | 备注 |
|---|---|---|---|
| H2O / STO-3G | max_ij \|J_CFMM−J_direct\| | ≤ 1e-8 Ha | `test_cfmm_j`(单次 J,无收敛 SCF) |
| C10H22 / STO-3G | \|E_CFMM−E_Direct\| | ≤ 1e-7 Ha | `regressions/scf/c10h22_cfmm_j` |
| 参数 | `cfmm_order=10`, `cfmm_threshold=1e-9` | — | WJGH 延伸 + M2L 截断 |

Oracle = 自家 Schwarz Direct(PySCF 无 FMM)。收紧 `cfmm_threshold` 误差单调不增。

### 交叉点规模表(标量;无柱状图)

方法学与机器环境钉于
`benchmarks/scf/validation/cfmm_scaling/README.md`(不进 ctest)。
跑完后把实测中位墙钟填入下表:

| n (C_nH_{2n+2}/STO-3G) | t_Direct [s] | t_CFMM [s] | 备注 |
|---|---|---|---|
| 4 | _TBD_ | _TBD_ | 小体系 Direct 通常更快 |
| 8 | _TBD_ | _TBD_ | |
| 10 | _TBD_ | _TBD_ | 回归体系 |
| 16 | _TBD_ | _TBD_ | 交叉点候选 |
| 交叉点 n* | — | — | 首个 t_CFMM < t_Direct 的 n |

### 误差控制语义

| 成员 | 含义 | 单位 |
|---|---|---|
| `cfmm_order` | 多极截断阶 p;远场 M2L 误差随 p 几何衰减 | 无量纲 |
| `cfmm_threshold` | 穿透阈值 ε_thr;控制 r_ext 与近/远场划分 | 无量纲电荷份额 |

## 章末:与 PySCF 对照(累积表)

精度(能量,Ha;PySCF 2.14.0,conv_tol=1e-12,球谐):

| 体系/设置 | 量 | Atomiverse | PySCF(钉定) | 容差 | 环 |
|---|---|---|---|---|---|
| H2O / def2-TZVP, RI-JK | E(RHF) | −76.058999347289 | −76.058999347289 | 1e-6(实测 3.3e-13) | 12 |
| H2O / def2-TZVP, 精确 ERI | E(RHF) | −76.059004760577 | −76.059004760577 | 1e-8 | 12 |
| NaCl rock-salt / Ewald+PME | Madelung M | 1.74756… | 文献 −1.7476/离子对 | 1e-4 rel | 19 |
| C10H22 / STO-3G, CFMM-J | \|E−Direct\| | ≤1e-7 | (PySCF 无 FMM;oracle=Direct) | 1e-7 | 21 |

速度(甘氨酸 / def2-SVP RHF,tol 1e-9,8 线程,3 次取中位;
方法学与机器环境钉于 benchmarks/scf/validation/ri_speed/README.md,
2026-07-21):

| 列 | 墙钟 [s] | E(RHF) [Ha] |
|---|---|---|
| ours-RI-JK(`rhf.ri = RI::JK`) | 1.617 | −282.615514331 |
| ours-exact(`RI::None`) | 2.516 | −282.615772070 |
| PySCF-DF(`density_fit()`) | 0.104 | −282.615514330 |

RI-JK 在 10 原子的甘氨酸上已胜精确 ERI 1.56×,且随体系增大而扩大;与
PySCF-DF 的 ~15× 差距是引擎整体速度(迭代数与积分内核),不是 RI 算法
差异——两者能量在同一拟合空间里对齐到 1e-9 Ha。

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
