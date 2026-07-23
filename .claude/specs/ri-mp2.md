---
slug: ri-mp2
status: approved
created: 2026-07-23
track: A-post-hf
grilled: false
---

# ri-mp2

## Summary

在 `RHF::mp2` 上交付 **RI/密度拟合 MP2 相关能**生产路径：复用已有 `internal::DensityFitting`（SCF RI-J/JK）的三中心积分与度量，避免完整 O(N⁴) AO ERI 物化，使中等基组 MP2 可跑；回归钉住相对 canonical MP2 / 文献或 PySCF DF-MP2 硬编码值。

## Domain basis

闭壳层 canonical MP2：

\[
E^{(2)} = \sum_{ijab}
\frac{(ia|jb)[2(ia|jb)-(ib|ja)]}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}
\]

RI 近似 \((ia|jb)\approx\sum_P B_{ia}^P B_{jb}^P\)，\(B\) 来自 Coulomb-metric 拟合（与 SCF `DensityFitting` 同构）。参考：Feyereisen, Fitzgerald, Komornicki, *Chem. Phys. Lett.* **208**, 359 (1993)；Weigend et al. 相关 DF 文献。单位：Hartree 内部；与 `MP2Result` 一致。

## Design

- **API**：`RHF::mp2(frame)` 在 `ri != RI::None` 时走 RI-MP2；或显式 `RHF::ri_mp2(frame)` 成员——优先 **同一 `mp2()`** 读 `rhf.ri` 成员（与 SCF 一致），禁止 factory。
- **Reuse** `internal::DensityFitting` 构建 B 张量；MP2 只消费 ov 块，不公开 IntegralPack。
- **错误**：未收敛 RHF、线性相关 aux 基 → throw（已有 lindep 路径）。
- **梯度**：Out of scope（后续 `mp2-gradients`）。

### Reuse decision

- `reuse` `RHF::mp2` / `MP2Result` / `ao2mo_ovov`（canonical 对照）。
- `reuse` `internal::DensityFitting` + `RI` enum + `reg_h2o_ri_jk` 辅助基。
- `reuse` `reg_h2o_mp2_sto3g` 作 canonical 对照门。
- `pattern` 方法成员 post-SCF（非 Component）。
- `generalize` DensityFitting 若仅暴露 J/K 接口 → 增加 ov 投影/B 张量访问供 MP2。

## Files to create or modify

- `include/atomiverse/cpu/scf/mp2.h` / `mean_field_method.h` / `rhf.h` 文档
- `src/cpu/scf/` mp2 + density_fitting 实现
- `include/atomiverse/cpu/scf/internal/density_fitting.h`（若需导出 B）
- `tests/cpu/scf/test_ri_mp2.cpp` (new)
- `regressions/scf/h2o_ri_mp2/reg_h2o_ri_mp2.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for RI-MP2 vs canonical on tiny system (tests/cpu/scf/test_ri_mp2.cpp)
- [ ] Generalize DensityFitting to supply RI-MP2 B/ov intermediates if needed
- [ ] Implement RI-MP2 energy path on RHF::mp2 when ri fitting enabled
- [ ] Document Doxygen (Hartree, RI aux basis)
- [ ] Add regression regressions/scf/h2o_ri_mp2 (hard-coded e_corr)
- [ ] Verify |E_RI − E_canonical| within documented tol on STO-3G/def2-SVP water
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：H2/STO-3G RI vs canonical < 1e-5 Ha（或实现标定）。
- **Regression**：H2O；硬编码 e_corr；注释 oracle；≤ 60 s。
- 不在 tests 跑大基组。

## Out of scope

- MP2 解析梯度 / 密度
- UMP2 / frozen core
- Laplace / local MP2
- CUDA

---

**Librarian note：** `cpu/scf/`；DensityFitting reuse critical；confidence high。
