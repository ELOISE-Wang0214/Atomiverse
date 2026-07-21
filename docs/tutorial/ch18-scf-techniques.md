# 第 18 章　自洽场的收敛

Roothaan–Hall 过程

\[
F(P_k)\,C = SC\,\varepsilon
\quad\rightarrow\quad
P_{k+1}
\]

是不动点迭代，不保证收敛。本章讨论三类相互衔接的手段：改善初猜、
改善迭代混合，以及在近自洽区改用二阶轨道更新。叙述顺序与实际推荐
流程一致：先把起点放进合理盆地，再在远区用稳健混合器，最后在近区
用高阶方法收尾。

---

## 18.1 为何会发散

常见失败模式有两类。

**电荷振荡。** 两个近简并构型 \(A\)、\(B\) 中，\(P_A\) 产生的 Fock 偏爱
\(B\)，\(P_B\) 又偏爱 \(A\)，密度做周期为 2 的振荡。小 HOMO–LUMO 能隙、
拉长键、双自由基倾向的体系最容易出现。

**过冲。** 每步更新过大，新密度比旧密度更远离不动点。裸核哈密顿初猜的
第一步往往过猛：没有任何电子屏蔽时，电子被过度压向核。

两种经典对策在前面章节已出现，机理简述如下。

- **密度阻尼** \(P\leftarrow(1-\alpha)P_{\mathrm{new}}+\alpha P_{\mathrm{old}}\)
  （成员 `damping_factor`；\(\alpha\) 为保留旧密度的比例）。抑制振荡，
  同时降低收敛速度。
- **能级移位**（Saunders–Hillier）：抬高虚轨道，人为增大能隙。\(\lambda\)
  越大越稳、越慢；收敛后引擎用未移位 Fock 重解一次轨道，以免把移位
  带进最终本征值。

更系统的做法是分阶段：好初猜（§18.2）→ 远区混合（EDIIS/ADIIS，§18.3）
→ 近区 CDIIS 或二阶方法（§18.3–18.4）。

---

## 18.2 初猜

SCF 是局部迭代，初猜决定落入哪个收敛域以及离解多远。按物理内容由少到
多排列：

| 初猜 | 要点 |
|---|---|
| Hcore | 对角化 \(T+V_{ne}\)，无电子屏蔽，最差但总是可用 |
| MINAO / SAD | 原子密度叠加；默认；密度合理但非某 Fock 的本征态 |
| 扩展 Hückel | 经验 \(F_{\mu\nu}\)（VSIP 与 Wolfsberg–Helmholz \(K=1.75\)），一次对角化得轨道 |
| SAP | 叠加原子有效势，\(F=T+V_{\mathrm{SAP}}\)，不依赖密度 |

Hückel 表覆盖 H–Ar；表外元素拒绝计算。SAP 覆盖 \(Z\le 118\)（Lehtola
等的有效势数据）。开壳层时，自旋对称的初猜算符不能单独产生自旋极化；
破缺仍依赖后续混合与占据。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2SVP;
rhf.init_guess = cpu::scf::InitGuess::Hueckel; // 或 SAP、MINAO、Hcore
rhf.run(frame);
```

水分子 def2-SVP、RHF、\(10^{-9}\) 容差下，Hückel 与 SAP 的迭代次数均低于
Hcore（约 14 对 16）。不同程序的 Hückel 参数化不同，应比较相对次序，
而非绝对迭代数。

---

## 18.3 混合：CDIIS、EDIIS 与 ADIIS

**CDIIS**（Pulay）以对易子误差 \(e_i=X^\top(F_i P_i S-S P_i F_i)X\) 做
无约束最小二乘，外推 Fock。近自洽时误差与密度近似线性，收敛接近二次。
远离自洽时线性化失效，系数可很大，一步把密度推出物理区域。

**EDIIS**（Kudin–Scuseria–Cancès）在单纯形 \(\sum c_i=1\)、\(c_i\ge 0\) 上
极小化能量的二次模型；新密度是历史密度的凸组合，不易过冲。**ADIIS**
（Hu–Yang）用类似约束的增广 Roothaan–Hall 能量模型，不显式需要能量
通道。二者适合**远区**；接近收敛时最优解退化为只保留最新一步，在误差
约 \(10^{-4}\) 附近失速。正确用法是升级切换：远区 ADIIS/EDIIS，
\(\max|e|\) 低于阈值后**永久**切回 CDIIS。

```cpp
cpu::scf::UHF uhf(5, 5);
uhf.mixer = cpu::scf::MixerKind::ADIIS_CDIIS;
uhf.mixer_switch_tol = 1e-2;
uhf.run(frame);
```

例：键长放大 2.6 倍的水分子、def2-SVP、UHF、Hcore 初猜。自旋对称解
不稳定；真极小为碎对称双自由基（\(\langle S^2\rangle\approx 1.97\)）。
纯 CDIIS 在 60 步内振荡；ADIIS→CDIIS 约 38 步收敛到
\(-75.718966\,\mathrm{Ha}\)（与 PySCF 一致到 \(10^{-6}\,\mathrm{Ha}\)）。

EDIIS 使用的能量通道在实现上滞后密度一步，远区启发式可接受；ADIIS
不依赖该近似。近区切换由 `mixer_switch_tol` 控制（默认 \(10^{-2}\)）。

---

## 18.4 二阶 SCF：轨道 Newton 法

前两节改善的是不动点路径。本质改动是把问题改写为**在轨道流形上极小化
能量**。保持 \(C^\top S C=I\) 的更新可写成

\[
C(\kappa)=C\,e^{\kappa},
\qquad
\kappa^\top=-\kappa,
\]

独立参数仅为占据–虚轨道块 \(\kappa_{ai}\)。能量对 \(\kappa\) 的二次展开
为 \(E_0+g^\top\kappa+\frac12\kappa^\top H\kappa+\cdots\)。闭壳层实轨道下
\(g_{ai}=4F_{ai}\)；Hessian 与稳定性分析中的 \({}^1(A+B)\) 相同。
Newton 步 \(\kappa=-H^{-1}g\) 在极小附近二次收敛。

直接构造 \(H\) 代价过高，且远区 \(H\) 未必正定。**CIAH**（共迭代增广
Hessian；Sun）求解增广本征问题

\[
\begin{pmatrix}0&g^\top\\ g&H\end{pmatrix}
\begin{pmatrix}1\\ \kappa\end{pmatrix}
=\lambda
\begin{pmatrix}1\\ \kappa\end{pmatrix},
\]

取最低根。Davidson 迭代只需要矩阵–向量积 \(H\cdot v\)，而该积正是
稳定性分析的 \(\sigma\)-向量。信赖域限制 \(\|\kappa\|_2\)；能量上升则
缩小信赖域并回退。指数映射用缩放–平方 Taylor 求值，保持正交归一。

Newton 步比一次 CDIIS 外推贵（含占据–虚轨道积分变换）。健康体系上
CDIIS 已够用。病理体系上，推荐与混合器同一逻辑：**先 DIIS，后 SOSCF**。

```cpp
cpu::scf::UHF uhf(5, 5);
uhf.soscf = true;
uhf.soscf_start = 0.25;   // max|DIIS 误差| 低于此值后切换
uhf.soscf_max_micro = 12;
uhf.soscf_trust = 0.4;
uhf.run(frame);
```

切换不早于第 1 次迭代：第 0 步先做 Roothaan 对角化，使密度型初猜与
系数一致。收敛后用未移位 Fock 再对角化一次，得到正则轨道能（SOSCF
过程中的本征值只是半正则的）。

对上述拉长水分子，约第 4 步切入 Newton，约 9 个宏迭代后到达碎对称
极小，总能与 PySCF `mf.newton()` 差小于 \(10^{-11}\,\mathrm{Ha}\)。
拉长 N₂（\(R=2.5\,\mathrm{Å}\)，6-31G，RHF）上 CDIIS 本身可收敛，
SOSCF 的 Newton 相约 4 步，与 PySCF 宏迭代数同量级。

当前 SOSCF 仅支持 HF（RHF/UHF）。Kohn–Sham 需要 XC 二阶核；在 RKS/UKS
上打开 `soscf` 会被拒绝。

---

## 数值对照摘要

| 体系 | 说明 |
|---|---|
| H₂O def2-SVP RHF | Hückel/SAP 迭代少于 Hcore |
| 2.6×H₂O def2-SVP UHF | ADIIS→CDIIS 收敛；纯 CDIIS 失速 |
| 同上 + SOSCF | 与 PySCF newton 能量一致到 \(10^{-8}\,\mathrm{Ha}\) 以上 |
| N₂ \(R=2.5\,\mathrm{Å}\) 6-31G RHF + SOSCF | 与 PySCF newton 能量一致 |

### 进一步阅读

- Jensen, *Introduction to Computational Chemistry*，§ SCF techniques。
- Pulay, *Chem. Phys. Lett.* **73**, 393 (1980)（DIIS）。
- Kudin, Scuseria, Cancès, *J. Chem. Phys.* **116**, 8255 (2002)（EDIIS）。
- Hu & Yang, *J. Chem. Phys.* **132**, 054109 (2010)（ADIIS）。
- Sun, arXiv:1610.08423（CIAH）。
- Helgaker, Jørgensen, Olsen, *Molecular Electronic-Structure Theory*，第 10 章。
