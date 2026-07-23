---
slug: qmmm-link-atoms
status: approved
created: 2026-07-23
track: A-multiscale
grilled: false
---

# qmmm-link-atoms

## Summary

为静电/机械 QM/MM 增加 **link-atom（H 封端）边界**：跨 QM–MM 共价切割时自动放置封端 H、修正 QM 区拓扑/电荷，并在公开 Potential 路径上给出可回归的能量与力。将 textbook `QmmmElectrostatic` / `QmmmMechanical` 从「无共价边界」推进到最小可用边界方案。

## Domain basis

标准加性 QM/MM 与 H 封端：在切断的 QM 重原子–MM 重原子键上，沿键向量放置 H，键长 \(r_{\mathrm{QH}} = g\, r_{\mathrm{Q-MM}}\)（\(g\) 为比例因子，典型 ~0.7 量级，以实现/文献钉死）。力回传到真实原子需链式法则。参考：Senn & Thiel, *Angew. Chem. Int. Ed.* **48**, 1198 (2009), DOI 10.1002/anie.200802019；Field, Bash, Karplus, *J. Comput. Chem.* **11**, 700 (1990)。单位：Å、eV、e。

## Design

- **类型**：扩展现有 `QmmmElectrostatic` / 或新 `QmmmLinkAtom` Potential 组合现有 QM+MM；优先 **扩展** `QmmmElectrostatic` 接受 `link_bonds`（QM tag, MM tag）列表，避免第三套平行 API。
- **QM 引擎**：v1 保持 RHF/STO-3G 或 DFTB2 与现组件一致；link H 仅存在于 QM 私有几何，不写入用户 Frame 永久拓扑（或明确写回策略并测试）。
- **禁止** 公共 `QmmmConfig` bag / factory。
- **错误**：非法切割（QM–QM、缺失原子）throw。

### Reuse decision

- `reuse` `QmmmElectrostatic` / `QmmmMechanical`（`cpu/multiscale/`）。
- `reuse` `RHF` / `DFTB2` 作为 QM 力源。
- `reuse` MM `PairLJCut` 力学部分。
- `pattern` 现回归 `reg_qmmm_es_h2` / `reg_qmmm_h2_ar`。
- `new` link 几何与力投影 helper（private）。

## Files to create or modify

- `include/atomiverse/cpu/multiscale/qmmm_electrostatic.h`（及/或 mechanical）
- `src/cpu/multiscale/` 对应实现
- `include/atomiverse/cpu/multiscale/link_atom.h` (new, private-ish helper ok as internal)
- `tests/cpu/multiscale/test_link_atom.cpp` (new)
- `regressions/multiscale/qmmm_link_ch3oh/reg_qmmm_link_ch3oh.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for link H placement and force projection (tests/cpu/multiscale/test_link_atom.cpp)
- [ ] Implement link-atom geometry + force back-transform helper
- [ ] Extend QmmmElectrostatic (or sibling) public API for link bonds
- [ ] Document Doxygen (Å, e, force chain rule)
- [ ] Add regression regressions/multiscale/qmmm_link_ch3oh (hard-coded energy or force gate)
- [ ] Verify one C–C or C–O cut textbook fragment
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：给定 Q–MM 坐标，H 位置解析公式。
- **Unit**：平移 QM+MM+link 体系合力 ~0。
- **Regression**：最小有机片段（如甲醇甲基切割玩具）；硬编码；≤ 60 s。

## Out of scope

- 极化 QM/MM（`polarizable-qmmm`）
- ONIOM 减性方案
- 周期 QM/MM + PME 长程
- CUDA

---

**Librarian note：** `cpu/multiscale/`；confidence high。
