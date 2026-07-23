---
slug: multi-cv-metad
status: approved
created: 2026-07-23
track: A-sampling
grilled: false
---

# multi-cv-metad

## Summary

将 well-tempered metadynamics 从 **1D**（`Metadynamics1D`）推广到 **多集体变量**（至少 2D）：在 \(\mathbf{s}\in\mathbb{R}^{n}\) 上堆积高斯丘，公开 Bias 可绑定多个 `CollectiveVariable`，回归在 2D 双阱或距离-角度玩具上验证丘计数与偏置力维度。

## Domain basis

Well-tempered metadynamics（Barducci et al.）：

\[
V(\mathbf{s},t)=\sum_{k:t_k\le t} h_k
\exp\!\Big(-\sum_{\alpha=1}^{n}\frac{(s_\alpha-s_{k,\alpha})^2}{2\sigma_\alpha^2}\Big),
\quad
h_k=h_0\exp\!\Big(-\frac{V(\mathbf{s}_k,t_k)}{(γ-1)k_BT}\Big)
\]

参考：Laio & Parrinello, *PNAS* **99**, 12562 (2002); Barducci, Bussi, Parrinello, *Phys. Rev. Lett.* **100**, 020603 (2008), DOI 10.1103/PhysRevLett.100.020603。单位：eV、CV 原生单位、K。

## Design

- **类型**：`cpu::bias::Metadynamics`（多维）或 `MetadynamicsND`；构造：`vector<unique_ptr<CV>>` + 每维 \(\sigma\) + \(h_0,T,γ,pace\)。
- **Reuse** `Metadynamics1D` 实现可 generalize 为 n=1 特例，或保留 1D 类为薄包装。
- **力**：对每个原子累加 \(\sum_\alpha (\partial V/\partial s_\alpha)(\partial s_\alpha/\partial x)\)。
- **禁止** factory；CV 仍为用户实现的 `CollectiveVariable` 协议。

### Reuse decision

- `reuse` `CollectiveVariable` / `DistanceCV` / `AngleCV` / `DihedralCV`。
- `generalize` `Metadynamics1D` → 共享丘存储与 WT 高度逻辑。
- `pattern` `HarmonicBias` 持有 CV 的所有权模型。
- `reuse` `reg_metad_double_well` 作 1D 回归不回归。

## Files to create or modify

- `include/atomiverse/cpu/bias/metadynamics.h` (new) 或扩展 `metadynamics_1d.h`
- `src/cpu/bias/metadynamics.cpp` (new)
- `include/atomiverse/cpu/bias/metadynamics_1d.h`（若改为包装）
- `tests/cpu/bias/test_metadynamics_nd.cpp` (new)
- `regressions/sampling/metad_2d/reg_metad_2d.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for 2D hill deposit and force size (tests/cpu/bias/test_metadynamics_nd.cpp)
- [ ] Generalize hill/WT logic from Metadynamics1D to N-dimensional Bias
- [ ] Implement public multi-CV metadynamics Bias API
- [ ] Document Doxygen (eV, sigma per CV, gamma)
- [ ] Add regression regressions/sampling/metad_2d (hard-coded n_hills / bias range)
- [ ] Verify 2D run deposits hills and finite forces
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：n=1 与 `Metadynamics1D` 数值一致（若 generalize）。
- **Unit**：2D 丘中心维度与 \(\nabla V\) 长度。
- **Regression**：双 CV 玩具；≤ 60 s；硬编码。

## Out of scope

- OPES / multiple-walkers / bias-exchange
- Multi-window WHAM（并行 spec）
- PLUMED 桥
- CUDA

---

**Librarian note：** `cpu/bias/`；confidence high。
