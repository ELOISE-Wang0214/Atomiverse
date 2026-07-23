---
slug: polarizable-qmmm
status: approved
created: 2026-07-23
track: A-multiscale
grilled: false
---

# polarizable-qmmm

## Summary

交付 **极化嵌入 QM/MM** textbook 路径：MM 位点携带可诱导偶极（或电荷-偶极），与 QM 密度/点电荷自洽极化，能量与力写回 Frame。在已有静电嵌入（固定 `point_charges`）之上增加 MM 极化自由度，不替代 link-atom spec。

## Domain basis

诱导偶极模型（Applequist / Thole 阻尼常见）：

\[
\boldsymbol{\mu}_i = \boldsymbol{\alpha}_i \mathbf{E}_i^{\mathrm{tot}},\quad
\mathbf{E}_i^{\mathrm{tot}} = \mathbf{E}_i^{\mathrm{QM}} + \sum_{j\neq i}\mathbf{T}_{ij}\boldsymbol{\mu}_j
\]

自洽迭代至 \(\max|\Delta\mu|<\varepsilon\)。能量含 \(-\frac12\sum_i\boldsymbol{\mu}_i\cdot\mathbf{E}_i^{\mathrm{perm}}\) 等标准项。参考：Thompson & Schenter, *J. Phys. Chem.* **99**, 6374 (1995)；Thole, *Chem. Phys.* **59**, 341 (1981)。单位：e、Å、eV。

## Design

- **类型**：`cpu::QmmmPolarizable` : `Potential`，或扩展 `QmmmElectrostatic` 的极化模式——优先 **新类** 保持机械/固定电荷路径简单。
- **QM**：v1 用 RHF + 外场/点偶极嵌入（reuse 现有 embedding 钩子若存在）。
- **MM**：各向同性 \(\alpha\) 成员参数 [Å³]；Thole 阻尼因子成员。
- **禁止** 公共 Config bag / factory。

### Reuse decision

- `reuse` 固定点电荷嵌入与 `QmmmElectrostatic` 分区逻辑。
- `reuse` `RHF` 外场/点电荷力路径（`reg_embed_h2o_grad` 家族）。
- `pattern` multiscale Potential + `last_qm_energy` 访问器。
- `new` 诱导偶极 SCC 循环（private）。

## Files to create or modify

- `include/atomiverse/cpu/multiscale/qmmm_polarizable.h` (new)
- `src/cpu/multiscale/qmmm_polarizable.cpp` (new)
- `tests/cpu/multiscale/test_qmmm_polarizable.cpp` (new)
- `regressions/multiscale/qmmm_pol_h2o/reg_qmmm_pol_h2o.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for induced-dipole SCC residual (tests/cpu/multiscale/test_qmmm_polarizable.cpp)
- [ ] Implement QmmmPolarizable Potential (header + src)
- [ ] Document polarizability units and damping members (Doxygen)
- [ ] Add regression regressions/multiscale/qmmm_pol_h2o (hard-coded energy shift vs non-polarizable)
- [ ] Verify polarization lowers energy vs fixed-charge embedding on toy system
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：单 MM 位点 + 均匀场 → \(\mu=\alpha E\)。
- **Unit**：SCC 收敛标志 / 迭代上限 throw。
- **Regression**：H2O QM + 偶极 MM 位点玩具；硬编码 \(\Delta E\)；≤ 60 s。

## Out of scope

- Drude 振子动力学 / 双恒温器
- AMOEBA 完整多极
- Link atoms（并行 spec）
- CUDA

---

**Librarian note：** `cpu/multiscale/`；confidence medium。
