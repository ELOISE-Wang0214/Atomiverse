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

## 20.2 布居分析与 ESP 电荷(骨架)

<!-- 填充:04-charges — Mulliken 回顾、Löwdin、MK/CHELPG 拟合与病态陷阱 -->

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
| _待 13 填充_ | H2O 频率 | — | — | 1 cm⁻¹ | 13 |
| _待 14 填充_ | CH2O 激发能 | — | — | 1e-4 eV | 14 |

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
