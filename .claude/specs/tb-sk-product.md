---
slug: tb-sk-product
status: approved
created: 2026-07-23
track: A-tb
grilled: false
---

# tb-sk-product

## Summary

把 SCC-DFTB2 从 **H-only 玩具表** 推进到 **静态转录的最小多元素 SK 子集**（至少 H/C/N/O 的 mio 风格参数，in-tree `.h` 表、运行时不读文件），并保证 **解析力** 可驱动 Evaluate/MD/Minimize。GFN 全表与 DFTB3 明确不在本 spec（可 follow-up）。

## Domain basis

SCC-DFTB2 总能量（Elstner et al.）：

\[
E = \sum_i n_i\varepsilon_i + E_{\mathrm{rep}} + E_{\mathrm{2nd}}[\Delta q]
\]

SK 积分表 \(H_{\mu\nu}(R), S_{\mu\nu}(R)\) 与斥力对势；SCC 用 \(\gamma_{AB}\) 与 Mulliken 电荷。参考：Elstner et al., *Phys. Rev. B* **58**, 7260 (1998), DOI 10.1103/PhysRevB.58.7260；mio 参数集文档（DFTB.org）。单位：公开 eV、Å、eV/Å。

## Design

- **数据**：`cpu/semiempirical/basis_data/` 或 `dftb_tables.h` 静态数组（对齐 GPW `gpw_tables.h` 模式）；`GthPotential` 式 **lookup 构造** 仅当需要按元素取表——SK 对表用内部函数 OK，**禁止** `builtin_mio()` 公共 factory 返回平行世界模型。
- **类型**：扩展 `DFTB2` 成员/表覆盖；非 H 元素不再 throw（在支持集合内）。
- **力**：`with_forces=true` 路径对多元素几何正确；FD 对照。
- **错误**：不支持的元素对 throw，无静默 0。

### Reuse decision

- `reuse` `DFTB2` 类与 `reg_h2_dftb2`。
- `reuse` `Semiempirical` kind + Evaluate/MD 绑定模式（`AM1Evaluator` 先例）。
- `pattern` GPW 静态表 + OOP lookup（`gth_potential.h`）。
- `new` multi-element SK 表数据（转录，非运行时下载）。

## Files to create or modify

- `include/atomiverse/cpu/semiempirical/dftb2.h`
- `include/atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h` (new)
- `src/cpu/semiempirical/dftb2.cpp`（及表相关）
- `tests/cpu/semiempirical/test_dftb2_sk.cpp` (new)
- `regressions/semiempirical/h2o_dftb2_mio/reg_h2o_dftb2_mio.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for C/H SK lookup and H2O energy finite (tests/cpu/semiempirical/test_dftb2_sk.cpp)
- [ ] Add static mio-like SK tables for H/C/N/O (in-tree headers)
- [ ] Implement multi-element DFTB2 energy + analytic forces using tables
- [ ] Document supported elements and units (Doxygen)
- [ ] Add regression regressions/semiempirical/h2o_dftb2_mio (hard-coded energy/forces)
- [ ] Verify H2O DFTB2 vs cited reference window
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：不支持元素 throw；H2 能量仍匹配旧门。
- **Unit**：力 FD 一致性（单分量）。
- **Regression**：H2O；硬编码；≤ 60 s。

## Out of scope

- DFTB3 / spin DFTB / 周期 DFTB
- 完整 GFN2-xTB 产品表（GFN0 toy 保持独立）
- CUDA
- 运行时解析 .skf 文件

---

**Librarian note：** `cpu/semiempirical/`；GPW table pattern；confidence high。
