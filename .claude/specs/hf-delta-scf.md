---
title: ΔSCF vertical ionization — worked example, tutorial section, PySCF regression
status: approved
created: 2026-07-02
---

# ΔSCF vertical ionization — worked example, tutorial section, PySCF regression

## Summary

本 spec 交付 ΔSCF 垂直电离能的完整教学工作流：新增算例 `examples/12_delta_scf.cpp`，对 He、H₂、H₂O、NH₃、CH₄（STO-3G 与 6-31G）在**中性几何**下分别求解中性 RHF 与阳离子 UHF，计算 I_ΔSCF = E_UHF(N−1) − E_RHF(N)，并与 Koopmans 定理的 −ε_HOMO 并排打印对比表（含轨道弛豫量）；在 `docs/tutorial/ch03-open-shell.md`（UHF 章）新增 Jensen 风格的 ΔSCF 小节，讲清为何 ΔSCF ≤ Koopmans（变分两行论证）、自旋记账、自旋污染与激发解两个陷阱，并交叉引用 ch02 的 Koopmans 小节；所有数字由引擎产生，并以 PySCF 2.13.1 钉死值（1e-7 Ha / ⟨S²⟩ 1e-4）做回归测试——引擎机制（`build_integrals` 的 `charge` 参数 + `uhf_solve`/`scf_solve`）已存在，本 spec 不改引擎，只补齐目前完全缺失的阳离子算例/测试覆盖。默认 ctest 门（~16 s）只保留 2+1 个锚点（单元级）；完整 5×2 PySCF 平价扫描属 e2e，放 `benchmarks/qc/rhf/validation/`（`val_scf_delta_scf_sweep`），验收时跑。

**依赖（链上游）：**
- **hf-koopmans**（前一 spec）：ch02 的 Koopmans 小节与 `examples/11_koopmans.cpp` 编号；且 **6-31G 基组表当前无 He**（已核实 `pople_631g.h` 第 39 行 `Z=2 unsupported`）——He/6-31G 条目由 hf-koopmans 负责补表。本 spec 的 He⁺/6-31G 扫描项以此为前置；若 hf-koopmans 未落地，E2E 扫描暂时跳过 He/6-31G 并留 TODO。
- **hf-robust-01**（若获批）：其 multiplicity resolver 可从 (charge, multiplicity) 推导 (n_α, n_β)。本 spec 采用手工自旋记账（本身就是教学点），并在算例注释中注明：resolver 落地后替换为一行调用。

## Domain basis

- **定义**：I_ΔSCF = E_UHF(阳离子, 中性几何) − E_RHF(中性)。阳离子为二重态：M = 2，n_α = n_β + 1（PySCF 侧 `charge=1, spin=1`）。两个态必须使用同一几何、同一基组、同一收敛阈值。
- **变分两行论证（ΔSCF ≤ Koopmans）**：Koopmans 的 −ε_HOMO 对应"冻结轨道"的 (N−1) 电子行列式能量；E_UHF(阳离子) 是对同一 Hilbert 空间变分优化的结果，故 E_relaxed ≤ E_frozen，即 I_ΔSCF ≤ −ε_HOMO。等号当且仅当无可弛豫自由度——He⁺、H₂⁺ 在 STO-3G 中是"基组内精确"的单电子问题，无处可弛豫，故两者**严格相等**（教学点）。
- **陷阱**：(i) 自旋污染 ⟨S²⟩ > 0.75，经既有 `SCFResult::spin_contamination`（= ⟨S²⟩ − S(S+1)）监测；(ii) SCF 可能收敛到激发/空穴移位解——用外部 `SCFInitialGuess` 以中性收敛密度播种阳离子（见 Design），或校验 aufbau 占据；(iii) 垂直 vs 绝热：6-31G 的 H₂O (10.98 eV)、CH₄ (13.48 eV) 只与垂直实验值 (12.6, 13.6 eV) 做定性比较。
- **几何（与 hf-koopmans 完全一致）**：H₂ 键长 0.7414 Å；H₂O：O (0, 0, 0.1173)、H (0, ±0.7572, −0.4692)；NH₃、CH₄ 标准几何；He 单原子。
- **钉死参考值**（PySCF 2.13.1，`cart=True`，`conv_tol=1e-12`；E⁺ [Ha] / I_ΔSCF [eV] / ⟨S²⟩）：

  | 体系 | STO-3G | 6-31G |
  |---|---|---|
  | He⁺ | −1.9317484501 / 23.838 / 0.750000 | −1.9936177758 / 23.444 / 0.750000 |
  | H₂⁺ | −0.5387095799 / 15.728 / 0.750000 | −0.5569099329 / 15.506 / 0.750000 |
  | H₂O⁺ | −74.6558987136 / 8.357 / 0.755202 | −75.5805492591 / 10.978 / 0.755280 |
  | NH₃⁺ | −55.1725944750 / 7.660 / 0.757355 | −55.8247003045 / 9.152 / 0.757187 |
  | CH₄⁺ | −39.2386120436 / 13.285 / 0.760065 | −39.6851175561 / 13.480 / 0.761532 |

- **断言**：对全部 10 对 (体系, 基组) 有 I_ΔSCF ≤ −ε_HOMO；He⁺/H₂⁺ STO-3G 严格相等（1e-9 Ha）。
- **参考文献**：T. Koopmans, *Physica* **1**, 104 (1934), DOI 10.1016/S0031-8914(34)90011-2；F. Jensen, *Introduction to Computational Chemistry*, 3rd ed., Wiley (2017), ISBN 978-1-118-82599-0（§3.4 轨道弛豫）；A. Szabo & N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996), §3.3；PySCF：Q. Sun *et al.*, *WIREs Comput. Mol. Sci.* **8**, e1340 (2018), DOI 10.1002/wcms.1340。
- **PySCF 生成命令**：`scf.UHF(gto.M(atom=geom, charge=1, spin=1, basis=b, cart=True)).run().e_tot`；`mf.spin_square()`。

## Design

不触碰引擎源码；全部交付物为算例、教程、测试三层。

- **阳离子求解路径（spec 决定并写入教程）**：算例与测试统一走 `cpu::scf::build_integrals(frame, basis, /*charge=*/+1)` + `scf_solve(ints_cation, uhf_config(ints_cation, n_alpha, n_beta), &guess)`，其中 `SCFInitialGuess` 以中性 RHF 收敛总密度播种：`P_alpha = P_beta = P_neutral / 2`（两自旋通道对半）。理由：密度热启动把 SCF 拉向基态离子解，规避激发/空穴移位解陷阱——这正是教程要演示的实践。教程同时说明裸 `uhf_solve`（Hcore 初猜、damping 0.7、DIIS 起于第 5 轮）对本组小体系同样收敛，但热启动是应当养成的习惯。
- **自旋记账**：中性 N 电子 → 阳离子 N−1 电子，n_β = (N−1)/2，n_α = n_β + 1（M = 2）。手工计算并在算例中以注释讲解；hf-robust-01 落地后换 resolver。
- **Koopmans 侧**：中性 `rhf_solve` 的 `RHFResult::eigenvalues` 取 ε_HOMO（索引 `n_occ−1`）。
- **单位**：全程 Hartree 计算，表格输出 eV，转换只用 `core/units.h` 的 `kHartree2eV`（禁止内联魔数）。
- **参考数据单一来源**：新头文件 `tests/cpu/scf/delta_scf_reference.h` 集中存放 5 个几何构造函数与 10 组钉死值（E⁺ [Ha]、I [eV]、⟨S²⟩），被两个测试文件共享；`tests/cpu/scf/pyscf_reference.py` 追加阳离子 UHF 生成块留证。
- **测试分层**（遵循 `tests/CMakeLists.txt` 尾部 `ATV_E2E_TESTS` 禁用机制）：
  - 默认门：`test_cpu_scf_delta_scf` — 锚点 He⁺/STO-3G、H₂O⁺/STO-3G，外加 H₂⁺/STO-3G（2 个基函数、毫秒级，为相等性断言所需）。
  - E2E 门：`test_cpu_scf_delta_scf_sweep` — 完整 5×2 扫描 + 全部不等式断言，注册后加入 `ATV_E2E_TESTS` 列表。
- **教程（内容用英文）**：在 `ch03-open-shell.md`（该章已讨论阳离子与自旋污染）追加 "ΔSCF: ionization energies beyond Koopmans" 小节：变分两行论证、自旋记账、`SCFInitialGuess` 热启动、⟨S²⟩ 监测、He/H₂ 相等教学点、垂直-绝热定性对比，交叉引用 ch02 的 Koopmans 小节与 `examples/12_delta_scf.cpp`；章内引用的数字必须与 `delta_scf_reference.h` 钉死值一致（docs 回归）。
- **算例注册**：`examples/CMakeLists.txt` 在 `if(ATV_BUILD_RHF)` 块内 `atv_add_cpp_example(12_delta_scf)`（11 号已被 hf-koopmans 占用）；`examples/README.md` 增加一行。
- 新文件公共符号按 Doxygen 规范（`@file`/`@brief`，能量参数注明 [Hartree]/[eV]）。

## Files to create or modify

- `examples/12_delta_scf.cpp` (new)
- `examples/CMakeLists.txt`
- `examples/README.md`
- `docs/tutorial/ch03-open-shell.md`
- `tests/cpu/scf/delta_scf_reference.h` (new)
- `tests/cpu/scf/test_delta_scf.cpp` (new)
- `tests/cpu/scf/test_delta_scf_sweep.cpp` (new)
- `tests/CMakeLists.txt`
- `tests/cpu/scf/pyscf_reference.py`

## Tasks

- [ ] Write shared reference header with the 5 geometries and 10 pinned PySCF triples (E⁺ [Ha], I [eV], ⟨S²⟩) (`tests/cpu/scf/delta_scf_reference.h`)
- [ ] Write failing anchor tests — He⁺/H₂⁺/H₂O⁺ STO-3G energy+⟨S²⟩ regression, ΔSCF ≤ −ε_HOMO, He/H₂ equality to 1e-9 Ha (`tests/cpu/scf/test_delta_scf.cpp`)
- [ ] Write failing E2E sweep tests — all 10 pairs, same tolerances + inequality for every pair (`tests/cpu/scf/test_delta_scf_sweep.cpp`)
- [ ] Extend PySCF provenance script with the cation UHF generation block (`tests/cpu/scf/pyscf_reference.py`)
- [ ] Register `test_cpu_scf_delta_scf` and `test_cpu_scf_delta_scf_sweep` in `tests/CMakeLists.txt`; add the sweep target to the `ATV_E2E_TESTS` disable list
- [ ] Implement `examples/12_delta_scf.cpp` — neutral RHF + `SCFInitialGuess`-seeded cation UHF per system/basis, print Koopmans-vs-ΔSCF table (eV via `kHartree2eV`) with relaxation column and ⟨S²⟩
- [ ] Register the example in `examples/CMakeLists.txt` (inside `ATV_BUILD_RHF`) and add its row to `examples/README.md`
- [ ] Write the ΔSCF tutorial section in `docs/tutorial/ch03-open-shell.md` (variational argument, spin bookkeeping, warm-start, contamination + excited-solution caveats, He/H₂ teaching point, vertical-vs-experiment framing; cross-ref ch02 Koopmans; quote only pinned numbers)
- [ ] Verify default ctest gate stays ~16 s (the sweep is not a ctest target) and `val_scf_delta_scf_sweep` passes under `-DATV_BUILD_BENCHMARKS=ON`
- [ ] Run full check + test suite

## Testing strategy

- **正常路径**：锚点测试对 He⁺/H₂⁺/H₂O⁺（STO-3G）断言 `scf_solve` 收敛、能量与钉死值差 ≤ 1e-7 Ha、`spin_contamination + 0.75` 与钉死 ⟨S²⟩ 差 ≤ 1e-4；E2E 扫描对全部 10 对重复同一断言。
- **领域验证**（`science.required`）：(a) 对每一对断言 I_ΔSCF ≤ −ε_HOMO（变分下界）；(b) He⁺/H₂⁺ STO-3G 断言 |I_ΔSCF − (−ε_HOMO)| ≤ 1e-9 Ha（单电子基组内精确 ⇒ 零弛豫）；(c) ⟨S²⟩ 回归本身即自旋污染的定量验证。
- **边界**：热启动与 Hcore 冷启动收敛到同一能量（1e-8 Ha，锚点体系，防激发解回归）；`converged == true` 在所有断言之前显式检查。
- **门控**：sweep 不是 ctest 目标，位于 `benchmarks/`，验收时运行。
- **docs 回归**：ch03 中引用的每个数值必须能在 `delta_scf_reference.h` 中逐字找到（评审时人工核对，acceptance ac-007）。

## Out of scope

- 引擎改动：multiplicity resolver（hf-robust-01）、aufbau 占据校验器、MOM/ΔSCF 激发态定向收敛——均为独立 spec。
- He/6-31G 基组表补行（hf-koopmans 的交付物；此处只消费）。
- 绝热电离能（需阳离子几何优化）、电子亲和能、UKS/DFT 的 ΔSCF——教程仅一句展望。
- CUDA 路径（`cuda::scf` 无 UHF 阳离子教学需求）与基准测试（benchmarks/）。
- 除交叉引用外不改动 `docs/tutorial/ch02-rhf.md`。
