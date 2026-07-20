# Chapter 20 — HF 性质与环境:从偶极矩到相对论

> **状态:骨架章**(hf-capability-gaps 链 01 环创建)。各节由标注的链环填充;
> 本章把"从收敛的波函数读出物理"讲成一条弧:期望值 → 布居 → 响应(频率/激发)
> → 环境(外场/嵌入/溶剂)→ 相对论,对应 Jensen 性质诸章的 HF 侧。

## 20.1 分子多极矩

收敛的波函数里能读出的第一个物理量是**偶极矩**。对单行列式,任何单电子
算符的期望值都是一个迹:⟨Φ|Ô|Φ⟩ = Tr(P·O) — 这正是 `SlaterDeterminant::expect`
的契约(双电子算符没有这种形式;把 Fock 矩阵塞进去得到的 Σε 会双计电子互斥)。
于是

μ = Σ_A Z_A R_A − Tr(P·r),  Θ_ab = ½ Σ q (3 r_a r_b − r² δ_ab)(无迹 Buckingham)。

**原点问题**:把原点平移 d,偶极变换为 μ′ = μ − q_tot·d — 中性分子原点无关,
带电体系则不然;四极矩同理带上偶极交叉项。所以约定必须随数值一起报告:
Atomiverse 取**核电荷中心**为原点并写入 `kMultipoleOrigin`,这样全局平移下
报告值不变(原点随体系移动)。

### 上手

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.run(frame);

cpu::MultipoleObserver obs{rhf.result_handle()};
// 生产用法:挂到宿主驱动的 end_of_step;单点后也可直接 launch 一次。
// frame.meta 里出现 kDipole [Debye]、kQuadrupole [D·Å]、kMultipoleOrigin [Å]。
```

molvis 渲染:把 `kDipole` 矢量画在 `kMultipoleOrigin` 上即得教科书的
"偶极箭头"图(场景 API 见 molvis 文档)。

积分内核是移位中心闭式(Obara & Saika 1986, DOI 10.1063/1.450106),
平移恒等式 M¹(C′) = M¹(C) + (C−C′)·S 由单元测试在 1e-13 钉住。

## 20.2 布居分析与 ESP 电荷

"这个原子带多少电荷"不是可观测量:密度矩阵是整体的,任何按原子的划分都
是一种**约定**。本节把三种约定排成一条弧 —— Mulliken(重叠对半分)→
Löwdin(先正交化再划分)→ ESP 拟合(放弃轨道空间,直接拟合分子外的
静电势)—— 各自的适用域和陷阱都不同。

### Mulliken 回顾与 Löwdin

Mulliken 布居 N_A = Σ_{μ∈A}(PS)_{μμ} 把重叠布居对半劈给两个中心 —— 任意,
且对基组出了名地敏感(ch13 已见)。Löwdin 的修正是先把基组对称正交化:

q_A = Z_A − Σ_{μ∈A} (S^{1/2} P S^{1/2})_{μμ}   (Löwdin, 1950)

S^{1/2} 把每个 AO 变成"最像原来自己"的正交组合(对称正交化在最小二乘意
义下离原基组最近),重叠布居随之自然归属。两个警告:

- **Löwdin 对逐 AO 缩放不不变**(Mulliken/ESP 则不变)。笛卡尔 d/f 函数的
  工程归一化各码不同(PySCF cart 裸 libcint 与 Atomiverse 引擎的 d 自重叠
  差一个常数因子),所以 Atomiverse 在**单位归一化的 AO 基**(S_ii = 1)里
  做 Löwdin —— 数值与工程约定无关,神谕脚本同一步骤。
- 它仍是"轨道空间划分":基组越弥散,"这个函数属于哪个原子"越没意义。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Pople631Gs;
rhf.tol = 1e-12;
rhf.run(frame);                                   // H2O

auto q  = cpu::lowdin_charge(rhf.result_handle(), frame);   // [e] 按原子
auto sz = cpu::lowdin_spin(rhf.result_handle(), frame);     // P^α−P^β 版本
```

### ESP 拟合电荷:MK 与 CHELPG

力场要的电荷有一个明确的物理目标:**再现分子外的静电势**。ESP 拟合直接
以此为目标函数 —— 在分子 vdW 面外采样 QM 势

V(r_k) = Σ_A Z_A/|R_A−r_k| − Σ_{μν} P_μν ⟨μ|1/|r−r′||ν⟩

(电子项即核吸引积分内核,单位点电荷放在格点上),再解带 Lagrange 乘子的
约束最小二乘(Σq = 分子净电荷,KKT 线性系统)。两种采样约定:

- **MK**(Besler-Merz-Kollman 1990):1.4/1.6/1.8/2.0 × Bondi 半径的嵌套壳
  层,每壳黄金角 Fibonacci 球面 ~1 点/Å²,落入他原子壳内的点弃掉。
- **CHELPG**(Breneman-Wiberg 1990):0.3 Å 立方格点,vdW 面外、距任一原子
  vdW 面 ≤2.8 Å 的带内取点。点多(H₂O ~8500 vs MK ~360),对分子取向更钝。

```cpp
cpu::MKCharges mk;                 // 参数即成员:mk.shell_scales, mk.density_per_ang2
auto q_mk = mk.fit(rhf.result_handle(), frame);        // [e]

cpu::ChelpgCharges chelpg;         // chelpg.spacing_ang = 0.3, chelpg.band_ang = 2.8
auto q_chelpg = chelpg.fit(rhf.result_handle(), frame);
```

**埋藏原子陷阱**:ESP 只在分子外采样,深埋原子(大分子内部的 sp³ 碳)的
1/r 列几乎与近邻的列共线 —— 最小二乘病态,电荷数值大而无物理意义,微小
的密度变化引起电荷大摆。Atomiverse 在 KKT 矩阵奇异值坍缩时自动切 SVD 最小
范数解(不发散,但**共线组内的电荷划分是任意的**——总和仍可靠)。RESP 的
分阶段拟合与等价约束就是为此发明的(本环未实现)。对称等价原子(H₂O 的两
个 H)的拟合电荷只在格点离散化水平相等(~1e-3 e):确定性格点不严格镜像
对称,这是约定的一部分,不是 bug。

molvis 可视化:`mk.grid(frame)` / `chelpg.grid(frame)` 返回 Å 单位的采样
点,直接作为点云叠在分子上,按 V(r_k) 上色即得教科书的 ESP 采样图(场景
API 见 molvis 文档)。

三法对照(H₂O RHF/6-31G*,PySCF 2.14.0 钉定,`regressions/scf/h2o_charges_pyscf`):
Löwdin 给 q_O ≈ −0.73 e,MK/CHELPG 都给 ≈ −0.81 e —— 轨道划分与"再现外部
静电势的有效电荷"本来就是两个目标,数值见章末表。

## 20.3 解析 Hessian 与振动频率(骨架)

<!-- 填充:13-cphf — CPHF 轨道响应、稳定性内核第三次复用、H2O 频率表、
     molvis 简正模可视化 -->

## 20.4 激发态:TDHF/RPA 与 TDA(骨架)

<!-- 填充:14-tdhf — 极化率极点、A/B 内核第四配置、triplet instability
     回链 ch03、甲醛 n→π* -->

## 20.5 外场与点电荷嵌入(骨架)

<!-- 填充:15-embed — 环境项进 Fock、有限场 vs 解析偶极恒等式、坍缩风险 -->

## 20.6 隐式溶剂(骨架)

<!-- 填充:16-solvent — 导体近似、空腔、x 约定分歧、非静电项诚实声明 -->

## 20.7 标量相对论:x2c(骨架)

<!-- 填充:18-x2c — 轨道收缩叙事、x2c-1e、ECP vs x2c 决策、picture-change 警告 -->

## 章末:与 PySCF 对照(累积表占位)

| 体系/设置 | 量 | Atomiverse | PySCF(钉定) | 容差 | 环 |
|---|---|---|---|---|---|
| H2O RHF/6-31G* | μ_z [a.u.] | −0.8757858700 | −0.8757858708 | 1e-6 a.u.(实差 8e-10) | 03 |
| H2O RHF/6-31G* | μ_z [D] | −2.226026 | −2.2260255980 | —(HF 高估 vs 1.85 D 实验) | 03 |
| H2O RHF/6-31G* | Löwdin q_O [e](归一化 AO) | −0.72840153 | −0.7284015337 | 1e-6 e(实差 5e-9) | 04 |
| H2O RHF/6-31G* | MK q_O [e](文档化格点) | −0.80709214 | −0.8070921421 | 1e-3 e(实差 8e-10) | 04 |
| H2O RHF/6-31G* | CHELPG q_O [e](文档化格点) | −0.80500612 | −0.8050061239 | 1e-3 e(实差 9e-10) | 04 |
| N2 RHF/6-31G* | −ε_HOMO [eV](Koopmans IP) | 16.646730 | 16.64672961 | 1e-3 eV(实差 4e-7) | 04 |
| N2(+) UHF/6-31G* | ΔSCF IP [eV] | 15.692600 | 15.69256822 | 1e-3 eV(实差 3e-5) | 04 |
| _待 13 填充_ | H2O 频率 | — | — | 1 cm⁻¹ | 13 |
| _待 14 填充_ | CH2O 激发能 | — | — | 1e-4 eV | 14 |

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
