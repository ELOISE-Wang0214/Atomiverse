# Chapter 18 — SCF 技巧:初猜、混合与二阶方法

> **状态:半章**(初猜与 DIIS 家族由 hf-capability-gaps-07-guess 填充;
> 18.4 二阶方法归 09-soscf 环)。本章成文后,读者读完应相当于读完 Jensen
> 《Introduction to Computational Chemistry》§3.8 前半("SCF techniques"):
> 以**方法为主干**,把算法与技巧织进收敛叙事,而不是罗列功能。

## 18.1 SCF 为什么会发散

Roothaan-Hall 迭代 \(F(P_k)\,C = S C \varepsilon \to P_{k+1}\) 是一个不动点
迭代,而不动点迭代天生不保证收敛。发散的典型病理只有一种,却以两副面孔出现:

- **电荷来回倒(charge sloshing)**:两个近简并的电子构型 A、B 之间,
  \(P_A\) 生成的 Fock 算符偏爱构型 B,\(P_B\) 生成的又偏爱 A。密度在两个
  构型之间做周期 2 的振荡,能量曲线呈锯齿。小 HOMO-LUMO 能隙的体系
  (拉伸的共价键、双自由基、过渡金属)几乎必然表现出这一行为。
- **过冲(overshoot)**:每步走得太远,新的密度比旧的离自洽点更远,
  振幅逐步放大直至数值溢出。裸 Hcore 初猜的第一步常是重灾区——没有任何
  电子屏蔽的核哈密顿量会把全部电子塞进最深的几个轨道。

最古老的两副镇痛剂在第 2、3 章已经用过,这里回顾其机理与成员位置:

- **密度阻尼(damping)**:\(P \leftarrow (1-\alpha)P_{\text{new}} +
  \alpha P_{\text{old}}\),成员 `damping_factor`(PySCF 约定,α 是**保留旧
  密度**的比例;`UHF` 默认 0.7,DIIS 激活后自动关闭)。它压制振荡,但把
  收敛速度也一并压扁——纯阻尼是线性收敛。
- **能级移位(level shift)**:把虚轨道整体抬高 λ(成员 `level_shift`
  [Hartree];Saunders & Hillier 1973)。人为拉大能隙可以强行止住电荷
  来回倒,代价同样是变慢:λ 越大越稳、越慢,且移位下的"收敛"要在去掉
  移位后重新验证(引擎在收敛后自动用未移位 Fock 重解一次轨道)。

现代的标准配方是**三段式**:好初猜把起点放进收敛域附近(18.2),
远域混合器(EDIIS/ADIIS)负责把远离自洽的迭代拖进近域(18.3),
近域再交给 CDIIS 收尾。三段各司其职,任何一段用错位置都会失速。

## 18.2 初猜家族

SCF 是局部迭代,初猜决定落进哪个收敛盆、以及离盆底多远。四种初猜按
"物理含量"递增排列(Lehtola, *JCTC* **15**, 1593 (2019), DOI
[10.1021/acs.jctc.8b01089](https://doi.org/10.1021/acs.jctc.8b01089) 对
它们做了系统评测):

- **Hcore**(`InitGuess::Hcore`):对角化裸核哈密顿量 \(T + V_{ne}\)。
  没有电子-电子屏蔽,轨道被过度压向核,是所有初猜里最差的——它存在的
  意义是永远可用(不依赖任何表格)。
- **MINAO / SAD**(`InitGuess::MINAO`,默认):原子密度叠加——把预先算好
  的原子块密度投影到分子基组,得到一个不错的分子密度,但它不是任何
  Fock 算符的本征态。PySCF 的默认路线。
- **扩展 Hückel**(`InitGuess::Hueckel`):组装一个经验 Fock 矩阵
  \(F_{\mu\mu} = h_{\mu\mu}\)(VSIP,价态电离势静态表,Hartree),
  \(F_{\mu\nu} = \tfrac12 K (h_{\mu\mu}+h_{\nu\nu}) S_{\mu\nu}\),
  \(K = 1.75\)(Wolfsberg-Helmholz),对角化一次即得初猜轨道。与 MINAO
  相反:密度粗糙,但轨道**有能级结构**——占据模式从一开始就合理。
  表内元素 H–Ar;表外元素直接抛异常(无回退值)。
- **SAP**(`InitGuess::SAP`):原子势叠加——有效单电子势取
  \(V_{\text{SAP}}(\mathbf r) = \sum_A V_A(|\mathbf r - \mathbf R_A|)\),
  \(V_A\) 是 Lehtola 表列的交换泛函原子径向有效势(本引擎采用其
  高斯拟合 erf 展开形式,数据转录自 PySCF 的 `sap_grasp_large`;
  Lehtola-Visscher-Engel, DOI
  [10.1063/5.0004046](https://doi.org/10.1063/5.0004046)),
  \(F = T + V_{\text{SAP}}\) 在 DFT 求积网格上装配,对角化一次,
  **完全不需要密度**。全周期表 Z ≤ 118。

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2SVP;
rhf.init_guess = cpu::scf::InitGuess::Hueckel;   // 或 InitGuess::SAP
rhf.run(frame);
```

**开壳陷阱(占据平均化)**:Hückel 与 SAP 都是自旋对称的算符,单次对角化
给不出自旋极化。开壳层运行沿用现有规则——同一组初猜轨道分别填 α/β,
并做 HOMO-LUMO 混合以便打破对称。寻找碎对称解时,初猜只负责"起点合理",
对称破缺仍靠混合角与后续迭代(见 18.3 的困难例)。

H2O/def2-SVP RHF(tol=1e-9)的实测迭代数(PySCF 2.14.0 同几何对照;
PySCF 的 `huckel` 是 Lehtola 的免参数即时原子计算变体,本引擎是同文献的
静态 VSIP 表变体,故只应比较**相对 Hcore 的次序**,不比较绝对值):

| 初猜 | Atomiverse 迭代数 | PySCF 迭代数 |
|---|---:|---:|
| Hcore | 16 | 10 |
| MINAO | 16 | 8 |
| Hückel | **14** | 8 |
| SAP | **14** | 8 |

两条断言 `iters(Hückel) < iters(Hcore)`、`iters(SAP) < iters(Hcore)` 由
`regressions/scf/h2o_scf_convergence` 长期钉住。

## 18.3 DIIS 家族:CDIIS、EDIIS、ADIIS

**CDIIS 回顾**(Pulay 1980/1982):以对易子
\(e_i = X^\top(F_iP_iS - SP_iF_i)X\) 为误差向量,在历史子空间里做
无约束最小二乘,再以同一组系数外推 Fock:\(\tilde F = \sum_i c_i F_i\),
\(\sum_i c_i = 1\) 但 **c 可以为负、可以远大于 1**。自洽点附近 \(e\) 与
密度误差近似线性,这一外推近乎二次收敛——CDIIS 是无可争议的近域冠军。
但在远离自洽处,线性化失效,外推系数动辄 ±10,一步把密度甩出物理区域。

**EDIIS**(Kudin-Scuseria-Cancès 2002, DOI
[10.1063/1.1470195](https://doi.org/10.1063/1.1470195))换了目标函数:
在**单纯形约束**(\(\sum c_i = 1\), \(c_i \ge 0\))下最小化能量的二次模型

\[
\min_{c}\; \sum_i c_i E_i \;-\; \tfrac12 \sum_{ij} c_i c_j
\langle D_i - D_j \,|\, F_i - F_j \rangle .
\]

系数非负意味着新密度永远是历史密度的**凸组合**——绝不过冲,单调把
能量往下压。**ADIIS**(Hu-Yang 2010, DOI
[10.1063/1.3304922](https://doi.org/10.1063/1.3304922))用增广
Roothaan-Hall 能量替代显式能量,同一单纯形约束:

\[
\min_{c}\; 2\sum_i c_i \langle D_i - D_n | F_n \rangle
+ \sum_{ij} c_i c_j \langle D_i - D_n | F_j - F_n \rangle ,
\]

只需要 (D, F) 历史、不需要能量通道,是 PySCF 推荐的默认远域混合器。

**适用域陷阱——远域专用。** 凸组合既是它们的美德也是它们的天花板:
接近自洽时最优解退化到"全部权重压在最新一步",EDIIS/ADIIS 在
DIIS 误差 ~1e-4 附近失速,永远磨不出 1e-9。**用到底是错误用法。**
标准配方是升级切换:远域用 ADIIS/EDIIS,一旦 max|DIIS 误差| 降到阈值
以下就**永久**切回 CDIIS——

```cpp
cpu::scf::UHF uhf(5, 5);
uhf.mixer = cpu::scf::MixerKind::ADIIS_CDIIS;  // 默认 CDIIS,行为零变化
uhf.mixer_switch_tol = 1e-2;                   // max|DIIS err| 低于即切换
uhf.run(frame);
```

困难例:2.6 倍拉伸的 H2O/def2-SVP UHF(Hcore 初猜,无阻尼)。这个几何
的自旋对称解是不稳定驻点(PySCF 稳定性分析确认其存在内禀不稳定性),
真正的 UHF 极小是碎对称双自由基解(⟨S²⟩ ≈ 1.97)。纯 CDIIS 在 60 步预算
内始终绕着不稳定对称解振荡;ADIIS→CDIIS 则在 38 步收敛到 PySCF 钉定的
碎对称能量 −75.71896632 Ha(≤1e-6 Ha):

```molplot preset="molplot" theme="auto" width="100%" aspect="16:9"
mark: {type: line, point: true, tooltip: true}
data:
  values:
    - {iter: 2, series: CDIIS, E: -62.204}
    - {iter: 4, series: CDIIS, E: -80.771}
    - {iter: 6, series: CDIIS, E: -73.833}
    - {iter: 8, series: CDIIS, E: -76.205}
    - {iter: 10, series: CDIIS, E: -75.478}
    - {iter: 14, series: CDIIS, E: -75.357}
    - {iter: 18, series: CDIIS, E: -75.358}
    - {iter: 22, series: CDIIS, E: -75.358}
    - {iter: 26, series: CDIIS, E: -75.357}
    - {iter: 30, series: CDIIS, E: -75.358}
    - {iter: 34, series: CDIIS, E: -75.357}
    - {iter: 38, series: CDIIS, E: -75.357}
    - {iter: 42, series: CDIIS, E: -75.357}
    - {iter: 2, series: ADIIS-CDIIS, E: -62.204}
    - {iter: 4, series: ADIIS-CDIIS, E: -71.732}
    - {iter: 6, series: ADIIS-CDIIS, E: -75.354}
    - {iter: 8, series: ADIIS-CDIIS, E: -75.213}
    - {iter: 10, series: ADIIS-CDIIS, E: -75.365}
    - {iter: 14, series: ADIIS-CDIIS, E: -75.704}
    - {iter: 18, series: ADIIS-CDIIS, E: -75.718}
    - {iter: 22, series: ADIIS-CDIIS, E: -75.719}
    - {iter: 26, series: ADIIS-CDIIS, E: -75.720}
    - {iter: 30, series: ADIIS-CDIIS, E: -75.719}
    - {iter: 34, series: ADIIS-CDIIS, E: -75.71897}
    - {iter: 38, series: ADIIS-CDIIS, E: -75.71897}
    - {iter: 42, series: ADIIS-CDIIS, E: -75.71897}
encoding:
  x: {field: iter, type: quantitative, title: "SCF iteration"}
  y: {field: E, type: quantitative, title: "E (Hartree)", scale: {zero: false}}
  color: {field: series, type: nominal, title: null}
```

CDIIS(上方平台)被不稳定对称解捕获;ADIIS 的凸组合先把能量单调压过
对称平台,切换回 CDIIS 后二次收敛锁定碎对称极小(下方平台)。

工程注记:EDIIS 需要能量通道 \(E_i\),本引擎经 `MixStage` 的能量引用
传入,而能量在 density 阶段才结算,故 \(E_i\) 滞后其 (F, D) 对一次迭代
——这是远域启发式可接受的近似,但也是上面困难例中 EDIIS 表现不如
ADIIS 的原因之一;ADIIS 完全不用能量,无此近似。二者失速的另一面
(近域必须切回 CDIIS)由 `mixer_switch_tol` 负责,阈值默认 1e-2。

## 18.4 二阶方法:SOSCF / CIAH(骨架)

<!-- 填充:09-soscf — 轨道旋转上的 Newton-Raphson、增广 Hessian、信赖域、
     为何二次收敛、何时划算;DIIS vs SOSCF 收敛曲线 -->

## 章末:与 PySCF 对照(累积表)

| 体系/设置 | 量 | Atomiverse | PySCF(钉定) | 容差 | 环 |
|---|---|---|---|---|---|
| H2O/def2-SVP RHF, Hcore | 迭代数 | 16 | 10 | 次序断言 | 07 |
| H2O/def2-SVP RHF, MINAO | 迭代数 | 16 | 8 | 次序断言 | 07 |
| H2O/def2-SVP RHF, Hückel | 迭代数 | 14 | 8 | 次序断言 | 07 |
| H2O/def2-SVP RHF, SAP | 迭代数 | 14 | 8 | 次序断言 | 07 |
| H2O/def2-SVP RHF | E (Ha) | −75.96098399 | −75.96098398708 | 1e-7 | 07 |
| 2.6×H2O/def2-SVP UHF, ADIIS→CDIIS | E (Ha) | −75.71896632 | −75.71896632261 | 1e-6 | 07 |
| 2.6×H2O/def2-SVP UHF, 纯 CDIIS | 60 步内收敛? | 否(失速) | —(对称解不稳定) | — | 07 |
| _待 09-soscf 填充_ | 难例迭代数/终能量 | — | — | — | 09 |

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
