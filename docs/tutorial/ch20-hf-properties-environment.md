# 第 20 章　分子性质与环境

收敛的单行列式波函数已经给出能量与轨道。本章讨论如何从中读出可观测
量，以及如何把分子置于外场、溶剂或相对论核哈密顿中。次序与物理层次
一致：单电子期望值 → 布居与静电势 → 响应对核位移与光场 → 静电环境 →
标量相对论修正。

---

## 20.1 多极矩

对单行列式，任一单电子算符 \(\hat O\) 的期望值为

\[
\langle \Phi | \hat O | \Phi \rangle = \operatorname{Tr}(P O).
\]

电偶极与无迹四极（Buckingham 约定）因此写成

\[
\boldsymbol\mu = \sum_A Z_A \mathbf R_A - \operatorname{Tr}(P\,\mathbf r),
\qquad
\Theta_{ab} = \tfrac12 \sum q\,(3 r_a r_b - r^2 \delta_{ab}).
\]

原点平移 \(\mathbf d\) 时，\(\boldsymbol\mu' = \boldsymbol\mu - q_{\mathrm{tot}}\mathbf d\)。
中性体系偶极与原点无关；带电体系则必须同时报告原点。本引擎取**核电荷
中心**为原点，并写入 `kMultipoleOrigin`，使整体平移下数值不变。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.run(frame);

cpu::MultipoleObserver obs{rhf.result_handle()};
// 也可挂在宿主驱动的 end_of_step；结果在 kDipole、kQuadrupole、kMultipoleOrigin
```

一电子多极积分采用 Obara–Saika 移位中心公式（*J. Chem. Phys.* **84**, 3963
(1986)）。平移恒等式 \(M^1(\mathbf C') = M^1(\mathbf C) + (\mathbf C-\mathbf C')S\)
在单元测试中核对到 \(10^{-13}\)。

水分子 RHF/6-31G* 的 \(\mu_z\) 与 PySCF 对照见章末表；HF 对水偶极的高估
（约 2.23 D，实验约 1.85 D）是方法本身的系统误差，不是实现问题。

---

## 20.2 布居与静电势电荷

原子电荷不是可观测量：全部分配方案都是约定。下面三条路线目标不同，
数值不可直接互换。

### Mulliken 与 Löwdin

Mulliken 布居把重叠矩阵对半分给两个中心，对基组选择敏感。Löwdin 先做
对称正交化 \(S^{1/2}\)，再按原子块求和：

\[
q_A = Z_A - \sum_{\mu\in A} (S^{1/2} P S^{1/2})_{\mu\mu}.
\]

对称正交化在最小二乘意义下最接近原基组，但仍是轨道空间划分；弥散基下
“函数属于哪个原子”本身就不明确。实现上，本引擎在 **\(S_{ii}=1\)** 的
归一化 AO 基中做 Löwdin，以避免笛卡尔 d/f 工程归一化约定的差异。

```cpp
auto q  = cpu::lowdin_charge(rhf.result_handle(), frame); // [e]
auto sz = cpu::lowdin_spin(rhf.result_handle(), frame);   // 自旋布居
```

### ESP 拟合：MK 与 CHELPG

力场需要的电荷往往以**再现分子外静电势**为目标。在 vdW 面外采样

\[
V(\mathbf r_k)
  = \sum_A \frac{Z_A}{|\mathbf R_A-\mathbf r_k|}
  - \sum_{\mu\nu} P_{\mu\nu}
    \langle\mu|\,1/|\mathbf r-\mathbf r'|\,|\nu\rangle,
\]

在 \(\sum q = q_{\mathrm{mol}}\) 约束下解最小二乘（Lagrange 乘子 / KKT）。

- **MK**（Besler–Merz–Kollman）：Bondi 半径的若干倍壳层上采样。
- **CHELPG**（Breneman–Wiberg）：立方格点，落在 vdW 面外的薄壳内。

```cpp
cpu::MKCharges mk;
auto q_mk = mk.fit(rhf.result_handle(), frame);

cpu::ChelpgCharges chelpg;
auto q_chelpg = chelpg.fit(rhf.result_handle(), frame);
```

埋藏原子使 \(1/r\) 列近乎线性相关，最小二乘病态。奇异值过小时采用 SVD
最小范数解：电荷总和仍可靠，共线子空间内的划分则无唯一物理解释。对称
等价原子上的拟合电荷仅在格点离散误差量级（约 \(10^{-3}\,e\)）内相等，
这是采样约定的一部分。

水分子 RHF/6-31G*：Löwdin 氧电荷约 \(-0.73\,e\)，MK/CHELPG 约 \(-0.81\,e\)。
前者来自轨道划分，后者来自外部势拟合，二者目标不同。

---

## 20.3 解析二阶导数与振动

Hellmann–Feynman 定理使**一阶**能量导数不必知道密度对核坐标的响应。
**二阶**导数没有这一便利：必须求解轨道响应 \(U^x\)。耦合微扰
Hartree–Fock（CPHF；Pople 等，1979）对每个核坐标 \(x\) 给出

\[
{}^1(A+B)\,U^x = -B^x.
\]

左端的 \({}^1(A+B)\) 与轨道稳定性分析、SOSCF 中的 Hessian 是同一算符的
不同用法。右端含骨架导数 Fock 与重叠导数；实现中每一项都有有限差分
核对。装配后的解析 Hessian 先与“解析梯度的有限差分”一致，再与外部
程序对照。有限差分 Hessian 生产线仍保留，二者汇入同一
`vib::HessianResult`。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.tol = 1e-12;
cpu::vib::HessianResult H = rhf.hessian(frame); // eV/Å²

// 质量加权、平动转动投影、对角化 → 频率 [cm⁻¹]
// 见 regressions/scf/h2o_rhf_freq
```

`hessian()` 是方法输出，与 `forces()` 并列；振动分析属于 vib 组件。
当前解析 Hessian 限于闭壳层 RHF；无虚轨道的体系（如 He/STO-3G）会拒绝
计算。水分子 RHF/6-31G* 的三个谐振动频率与 PySCF 偏差约 \(10^{-4}\,\mathrm{cm}^{-1}\)
量级（章末表）。

---

## 20.4 时间相关 HF 与激发能

基态期望值回答“分子处于何种态”。激发能来自另一类响应：对频率 \(\omega\)
的振荡场，极化率 \(\alpha(\omega)\) 在激发能处出现极点。时间相关
Hartree–Fock（TDHF，亦称 RPA）把该自持密度振荡写成

\[
\begin{pmatrix} A & B \\ -B & -A \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \omega
\begin{pmatrix} X \\ Y \end{pmatrix}.
\]

\(X\) 为占据→虚轨道激发振幅，\(Y\) 为退激发振幅。矩阵 \(A\)、\(B\) 与
轨道 Hessian 内核相同，与稳定性、SOSCF、CPHF 共用一份 \(\sigma\)-向量
实现。

实轨道下可约化为对称本征问题

\[
(A-B)^{1/2}(A+B)(A-B)^{1/2} Z = \omega^2 Z.
\]

**Tamm–Dancoff 近似（TDA）** 丢掉 \(B\)（即丢掉 \(Y\)），得到 \(A X = \omega X\)，
在 HF 参考上即 CIS。变分性质保证每个对应根 \(\omega_{\mathrm{RPA}} \le
\omega_{\mathrm{TDA}}\)。

若参考态在三重态方向不稳定，则 \({}^3(A+B)\) 出现负本征值，RPA 的
\(\omega^2\) 可为负——这不是数值故障，而是微扰参考失效的信号。TDA 的
\(A\) 保持厄米，谱恒为实数。因此 `tdhf()` 在计算前做稳定性检查；不稳
定时抛出说明，并建议改用 `tda()` 或沿稳定性给出的方向重收敛。

跃迁偶极（长度规范）与振子强度

\[
\boldsymbol\mu_{0n}
  = \sqrt2 \sum_{ia}(X+Y)_{ia}\langle\phi_i|\mathbf r|\phi_a\rangle,
\qquad
f = \tfrac23\,\omega[\mathrm{Ha}]\,|\boldsymbol\mu_{0n}|^2[\mathrm{a.u.}].
\]

三重态通道自旋禁阻，\(f=0\) 是严格结果。

甲醛是常用的教学分子：最低单重激发 \(n\to\pi^*\)（\(^1A_2\)）偶极禁阻，
\(f=0\)；\(\pi\to\pi^*\) 为沿 C=O 轴的强吸收。单激发理论的能级通常偏高
约 0.5 eV，适用于定性归属与半定量比较，不宜当作实验光谱的精确预测。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.tol = 1e-12;
rhf.run(frame); // 甲醛几何见 regressions/scf/ch2o_tdhf

auto tda  = rhf.tda(3);
auto tdhf = rhf.tdhf(3);
```

---

## 20.5 外电场与点电荷嵌入

静电环境不引入新的相互作用类型，而是改写**单电子**哈密顿。点电荷嵌入
对应 QM/MM 中 QM 侧的静电势；均匀外场对应偶极耦合。双电子积分与 Frame
布局不变。

点电荷（位置以 Å 输入，电荷以 \(e\) 计）：

\[
h_{\mu\nu} \leftarrow h_{\mu\nu}
  - \sum_C q_C
    \langle\mu|\,1/|\mathbf r-\mathbf R_C|\,|\nu\rangle,
\]

经典核–电荷与电荷–电荷项计入核排斥能。积分与核吸引原语相同。

均匀电场 \(\mathbf F\)（原子单位），原点 \(\mathbf R_0\) 取核电荷中心
（与 §20.1 一致）：

\[
h_{\mu\nu} \leftarrow h_{\mu\nu}
  + \mathbf F\cdot\langle\mu|(\mathbf r-\mathbf R_0)|\nu\rangle,
\qquad
E_{\mathrm{nuc}} \leftarrow E_{\mathrm{nuc}}
  - \sum_A Z_A \mathbf F\cdot(\mathbf R_A-\mathbf R_0).
\]

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.point_charges = {
    {0.0, 0.0,  3.0,  0.5},
    {0.0, 0.0, -3.0, -0.5},
};
rhf.external_field = {0.0, 0.0, 1e-4};
rhf.run(frame);
// 含环境时 forces() 尚未实现，会抛出异常
```

一阶关系 \(E(\mathbf F)=E(0)-\boldsymbol\mu\cdot\mathbf F\) 给出交叉检验：
小场中心差分 \(\partial E/\partial F\) 应等于解析偶极。水/STO-3G、
\(F_z=10^{-4}\) a.u. 时二者在 \(10^{-5}\) a.u. 内一致。

静电场下能量无下界：场过强或基组过弥散时，电子可沿 \(-\mathbf F\) 方向
“电离”，SCF 失去物理意义。有限场计算应保持 \(|\mathbf F|\sim 10^{-4}\) a.u.
量级。

与 PySCF `qmmm.mm_charge` 对照时需注意：PySCF 默认**不含** MM–MM 自能，
本引擎将其计入核排斥；比较总能量时应把该项单独加上（见
`regressions/scf/h2o_embedding_field`）。

---

## 20.6 连续介质溶剂：C-PCM

溶剂可近似为介电常数 \(\varepsilon\) 的连续介质，溶质置于原子中心球
组成的空腔内。空腔表面的表观电荷响应溶质静电势，并在每次 SCF 迭代中
与密度一起更新——这与 §20.5 的固定点电荷不同。

在表面元（tessera）离散下，

\[
S\mathbf q = -f(\varepsilon)\,\mathbf V,
\qquad
f(\varepsilon)=\frac{\varepsilon-1}{\varepsilon+x},
\qquad
\Delta G_{\mathrm{el}}=\tfrac12\sum_k q_k V_k.
\]

\(x=0\) 为 C-PCM（Barone–Cossi；PySCF 默认），\(x=1/2\) 为原始 COSMO
（Klamt）。低 \(\varepsilon\) 时两者可差百分之几；水（\(\varepsilon=78.4\)）
下几乎重合。空腔半径取 Bondi 半径乘以缩放因子（默认 1.2）；落在其他
球内的格点剔除。表面在核坐标上不连续，几何扫描可能出现能量折角；
平滑开关（SWIG）不在本实现范围内。

本实现只含**静电** \(\Delta G_{\mathrm{el}}\)。空穴化、色散与排斥等非
静电能未计入，因此不能与实验水合自由能直接数值等同。溶剂打开时
`forces()` 未实现。

```cpp
cpu::scf::RHF sol;
sol.basis = cpu::scf::BasisSet::Pople631Gs;
sol.solvent_eps = 78.4;
sol.solvent_radii_scale = 1.2;
sol.solvent_x = 0.0; // C-PCM
sol.run(frame);
double dg_kcal = sol.e_solvent() * 627.509474;
```

`solvent_eps = 0`（默认）与气相路径一致。水 RHF/6-31G* 的 \(\Delta G_{\mathrm{el}}\)
约 \(-8.0\,\mathrm{kcal\,mol}^{-1}\)，与 PySCF C-PCM 相差约 \(0.06\,\mathrm{kcal\,mol}^{-1}\)
（章末表；完整设置见 `regressions/scf/h2o_rhf_631gs_cpcm`）。

---

## 20.7 标量相对论：x2c-1e

重核附近电子局域动量大，Dirac 理论给出 s/p 轨道径向收缩与能量降低，
d/f 轨道则因屏蔽变化而膨胀。对第三周期及更重元素的全电子计算，忽略
标量相对论会系统抬高内层能量。自旋–轨道耦合是另一通道，此处不讨论。

精确二分量（x2c）将四分量 Dirac 问题一步退耦到正能二分量空间
（Kutzelnigg–Liu；Peng–Reiher）。**x2c-1e** 只改造一电子哈密顿：在收缩
AO 基上构造 \(\{S,T,V,W\}\)，其中

\[
W_{\mu\nu}=\langle\nabla\mu|V|\nabla\nu\rangle,
\]

组装动能平衡基下的 Dirac 矩阵，取正能解得退耦矩阵 \(X\) 与重整化 \(R\)，
得到 \(h^{\mathrm{X2C}}\) 并替换 SCF 的核哈密顿。重叠、双电子积分与
Fock 迭代形式不变（二电子相对论修正不在 x2c-1e 定义内）。光速取
\(c=137.03599967994\) a.u.，与 PySCF 一致。

| | 有效核势 (ECP) | x2c-1e |
|---|---|---|
| 适用 | 重元素价层（如 def2-ECP） | 全电子基组覆盖的元素 |
| 机制 | 相对论效应参数化进有效势 | 显式改造一电子哈密顿 |
| 互斥 | 不得与 x2c 同时使用 | 不得与 ECP 同时使用 |

相对论 ECP 已含标量相对论；再叠 x2c 会双重计数。引擎在入口检查并拒绝
该组合。含环境或 x2c 时解析力尚未推导，`forces()` 会拒绝执行。

用 x2c 轨道计算偶极、ESP 等期望值时，若属性算符仍用非相对论形式，则
含 **picture-change** 误差。本实现不变换属性算符；此类数值仅可粗估，
不能当作 picture-change 修正后的结果。

```cpp
cpu::scf::RHF x2c;
x2c.basis = cpu::scf::BasisSet::Def2TZVP;
x2c.x2c = true;
x2c.run(ar_frame);
```

氩原子 RHF/def2-TZVP 的相对论能量降低约 \(1.84\,\mathrm{Ha}\)，与
PySCF `x2c1e()`（`xuncontract=False`）在 \(10^{-5}\,\mathrm{Ha}\) 内一致。

---

## 与外部程序的数值对照

下列数值来自本库回归与预先写入的 PySCF 2.14.0 参考（方法学见
[附录](appendix-pyscf-oracle.md)）。构建与测试不依赖 Python。

| 体系 | 量 | 本库 | PySCF | 容差 |
|---|---|---|---|---|
| H₂O RHF/6-31G* | \(\mu_z\) / a.u. | −0.87578587 | −0.87578587 | \(10^{-6}\) |
| H₂O RHF/6-31G* | Löwdin \(q_\mathrm{O}\) / e | −0.7284 | −0.7284 | \(10^{-6}\) |
| H₂O RHF/6-31G* | MK \(q_\mathrm{O}\) / e | −0.8071 | −0.8071 | \(10^{-3}\) |
| H₂O RHF/6-31G* | 谐频 / cm⁻¹ | 1843.7, 4036.2, 4146.4 | 同 | 1 |
| CH₂O RHF/6-31G* | TDA \(n\to\pi^*\) / eV | 4.6690 | 4.6690 | \(10^{-4}\) |
| CH₂O RHF/6-31G* | TDHF \(n\to\pi^*\) / eV | 4.4999 | 4.4999 | \(10^{-4}\) |
| CH₂O RHF/6-31G* | TDHF \(\pi\to\pi^*\) \(f\) | 0.1796 | 0.1796 | \(10^{-3}\) |
| H₂O STO-3G + 点电荷 | \(E\)（含 MM–MM）/ Ha | 见回归 | 见 REFERENCE | \(10^{-8}\) |
| H₂O 6-31G* C-PCM | \(\Delta G_\mathrm{el}\) / kcal mol⁻¹ | −8.02 | −7.96 | 0.1 |
| Ar def2-TZVP x2c-1e | \(E\) / Ha | 见回归 | −528.64422317 | \(10^{-5}\) |

### 进一步阅读

- Jensen, *Introduction to Computational Chemistry*，分子性质与溶剂化、相对论章节。
- Szabo & Ostlund, *Modern Quantum Chemistry*，§3 期望值；§7 时间相关方法。
- Dreuw & Head-Gordon, *Chem. Rev.* **105**, 4009 (2005)（激发能综述）。
- Barone & Cossi, *J. Phys. Chem. A* **102**, 1995 (1998)（C-PCM）。
- Peng & Reiher, *Theor. Chem. Acc.* **131**, 1081 (2012)（x2c 综述）。
