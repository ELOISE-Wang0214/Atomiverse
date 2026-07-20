# Chapter 19 — 大体系 HF:从 Schwarz 筛选到线性标度库仑

> **状态:骨架章**(hf-capability-gaps 链 01 环创建)。各节由标注的链环填充;
> 本章对应 Jensen《Introduction to Computational Chemistry》§3.8 后半
> ("the SCF bottleneck" / linear scaling):方法为主干,把 direct SCF、密度
> 拟合与多极方法编成一条"库仑瓶颈如何被逐级击破"的叙事。

## 19.1 SCF 的瓶颈:四指标 ERI(骨架)

<!-- 填充:12-df 导言 — O(N⁴) 张量、direct SCF + Schwarz 筛选(现有实现,
     也正是 PySCF 的分子路线;诚实注明 PySCF 无 FMM) -->

## 19.2 密度拟合 RI(骨架)

<!-- 填充:12-df — 4-指标→3-指标压缩、Coulomb 度规二阶误差、JKFIT 取舍、
     μHa 误差预算表、精度+速度对照表 -->

## 19.3 长程库仑 I:Ewald 与 PME(骨架)

<!-- 填充:19-longrange — 周期库仑的条件收敛、实/倒空间劈分、PME O(N log N);
     NaCl Madelung −1.7476 动手;longrange/ 新居所导览 -->

## 19.4 长程库仑 II:快速多极方法(骨架)

<!-- 填充:20-fmm — 多极压缩远场、M2M/M2L/L2L、与 PME 的开边界对偶;
     误差-阶数收敛表(标量表格);molplot 标度曲线 -->

## 19.5 长程库仑 III:CFMM 加速 HF(骨架)

<!-- 填充:21-cfmm — 高斯分布延伸分类、CFMM-J 接入、硬门(逐元素复现
     direct J)、交叉点规模表;direct→RI→CFMM 叙事收束 -->

## 章末:与 PySCF 对照(累积表占位)

| 体系/设置 | 量 | Atomiverse | PySCF(钉定) | 容差 | 环 |
|---|---|---|---|---|---|
| _待 12-df 填充_ | RI-JK 能量 / 墙钟 | — | — | — | 12 |
| _待 21-cfmm 填充_ | CFMM-J 精度 / 交叉点 | — | — | — | 21 |

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
