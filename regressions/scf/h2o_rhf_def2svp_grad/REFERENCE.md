# Reference — H2O RHF/def2-SVP analytic gradient, spherical harmonics (5d)

Oracle protocol: docs/tutorial/appendix-pyscf-oracle.md (external run, pinned).

| 字段 | 值 |
|---|---|
| 几何 (Å) | O(0, 0, 0.1173); H(0, ±0.7572, −0.4692) |
| 基组 | def2-SVP, **spherical** (`cart=False`; the set's defining convention) |
| 收敛 | `conv_tol = 1e-12`, RHF |
| PySCF 版本 | 2.14.0 |
| 命令 | `gto.M(atom=..., unit='Angstrom', basis='def2-svp', cart=False)`; `mf=scf.RHF(mol); mf.conv_tol=1e-12; mf.kernel()`; `g = mf.nuc_grad_method().kernel()` |
| **E_RHF** | **−75.96098398708 Ha** (identical to the energy sibling — same run) |
| 对照容差 | 元素级 1e-6 Ha/Bohr;平移不变性 1e-8 Ha/Bohr |

## Pinned analytic gradient dE/dR [Ha/Bohr]

`mf.nuc_grad_method().kernel()` (PySCF has **no** `kick_start`; the entry point
is `.kernel()`). Molecule lies in the yz-plane, so all x components are
numerically zero (~1e-16).

| 原子 | ∂E/∂x | ∂E/∂y | ∂E/∂z |
|---|---|---|---|
| O  | 0            | 0             | +1.804450607022e-02 |
| H  | 0            | +1.107031890619e-02 | −9.022253035110e-03 |
| H  | 0            | −1.107031890619e-02 | −9.022253035110e-03 |

Σ over atoms = (−3e-30, +4.2e-15, −4.2e-15) Ha/Bohr — translationally invariant
to ~1e-15 (the regression gate is a looser 1e-8, covering the code's assembly).

`max|∂E/∂R| = 1.804450607022e-02 Ha/Bohr`.

## 为什么 1e-6 是判别阈值

PySCF 与 Atomiverse 都给**解析**梯度,同一 c2s 约定几何无关。正确的 c2s
sandwich(1e 路径)+ 密度反变换(2e 路径)复现 PySCF 到机器精度量级;错误的
约定/缺失的 Pulay 或算符导数项会产生 ≥ 1e-3 Ha/Bohr 的偏差,远高于本阈值。
d 壳层(l=2,6→5)在此几何上被激活,是约定错误最先暴露之处。

文献:Weigend & Ahlrichs, *PCCP* **7**, 3297 (2005), DOI 10.1039/B508541A
(def2-SVP 定义,球谐约定);Schlegel & Frisch, *Int. J. Quantum Chem.* **54**,
83 (1995), DOI 10.1002/qua.560540202(c2s 闭式系数);Pulay(能量加权密度 /
Pulay 力);Helgaker, Jørgensen, Olsen, *Molecular Electronic-Structure
Theory* §10(分子梯度)。
