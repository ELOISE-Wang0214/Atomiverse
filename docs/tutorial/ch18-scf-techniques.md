# Chapter 18 — SCF 技巧:初猜、混合与二阶方法

> **状态:骨架章**(hf-capability-gaps 链 01 环创建)。各节由标注的链环填充;
> 本章成文后,读者读完应相当于读完 Jensen《Introduction to Computational
> Chemistry》§3.8 前半("SCF techniques"):以**方法为主干**,把算法与技巧织进
> 收敛叙事,而不是罗列功能。

## 18.1 SCF 为什么会发散(骨架)

<!-- 填充:07-guess — 振荡/电荷来回倒的病理学,damping 与 level shift 回顾
     (现有实现),为什么"好初猜 + 远域混合器 + 近域 DIIS"是标准三段式 -->

## 18.2 初猜家族(骨架)

<!-- 填充:07-guess — Hcore 为何差、MINAO/SAD 思想、扩展 Hückel(VSIP)、
     SAP(Lehtola 原子势叠加);H2O/def2-SVP 各初猜迭代数表 -->

## 18.3 DIIS 家族:CDIIS、EDIIS、ADIIS(骨架)

<!-- 填充:07-guess — Pulay 对易子 DIIS 回顾、EDIIS/ADIIS 能量插值、
     适用域陷阱(远域专用,近收敛必须切回 CDIIS)、切换配方;
     molplot 收敛曲线 -->

## 18.4 二阶方法:SOSCF / CIAH(骨架)

<!-- 填充:09-soscf — 轨道旋转上的 Newton-Raphson、增广 Hessian、信赖域、
     为何二次收敛、何时划算;DIIS vs SOSCF 收敛曲线 -->

## 章末:与 PySCF 对照(累积表占位)

| 体系/设置 | 量 | Atomiverse | PySCF(钉定) | 容差 | 环 |
|---|---|---|---|---|---|
| _待 07-guess 填充_ | 迭代数(各初猜) | — | — | — | 07 |
| _待 09-soscf 填充_ | 难例迭代数/终能量 | — | — | — | 09 |

对照方法学见[附录:PySCF 对照方法学](appendix-pyscf-oracle.md)。
