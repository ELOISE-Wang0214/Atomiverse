---
slug: multi-window-wham
status: approved
created: 2026-07-23
track: A-sampling
grilled: false
---

# multi-window-wham

## Summary

交付 **多窗口 WHAM** 自由能重构：引擎侧从多窗 umbrella（`HarmonicBias` + CV）流式写出 CV 样本；**WHAM 数学在 molrs**（与 RDF 归属一致），C++ 仅 FFI/薄封装。回归用解析双阱或硬编码合成直方图验证 \(F(s)\) 重建。

## Domain basis

WHAM（Kumar et al. 1992；Souaille & Roux 2001）：自洽

\[
P(s) \propto \sum_i n_i(s) \Big/ \sum_j N_j e^{-\beta(w_j(s)-f_j)},\quad
e^{-\beta f_i} = \int P(s) e^{-\beta w_i(s)}\,ds
\]

\(w_i=\frac12 k(s-c_i)^2\)。参考：Kumar et al., *J. Comput. Chem.* **13**, 1011 (1992), DOI 10.1002/jcc.540130812。单位：eV 或 kT（文档固定一种并钉住）。

## Design

- **Engine**：umbrella 已有；本 spec 增加 **多样本导出**（文件或 molrs accumulator API）——Observer 流式写 CV，不在 C++ 重实现 WHAM 迭代。
- **molrs**：`wham` 模块（Rust）自洽求解 \(f_i\) 与 \(F(s)=-\beta^{-1}\ln P(s)\)；经现有 molrs CXX 桥暴露。
- **C++ 门面**：`cpu::compute::wham` 或 `cpu::Wham` 薄封装（PIMPL 如 trajectory_metrics），禁止在 `cpu/bias` 内嵌完整 WHAM 数值核。
- **架构**：遵守 “analysis math in molrs”。

### Reuse decision

- `reuse` `HarmonicBias` + CV + `reg_umbrella_*`。
- `reuse` molrs 流式 accumulator 模式（RDF/MSD）。
- `pattern` `cpu/observe/trajectory_metrics.h` PIMPL。
- `new` molrs WHAM solver + C++ FFI 包装。

## Files to create or modify

- `molrs/` WHAM 源（Rust，路径按 molrs 布局）
- molrs CXX 桥 / `molrs.hpp` 生成面（若需新 API）
- `include/atomiverse/cpu/observe/wham.h` (new) 或 `compute` 命名空间
- `src/cpu/observe/wham.cpp` (new)
- `tests/cpu/observe/test_wham.cpp` (new) — 单函数包装；合成数据
- `regressions/sampling/wham_double_well/reg_wham_double_well.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for WHAM free-energy rebuild on synthetic histograms (tests/cpu/observe/test_wham.cpp)
- [ ] Implement WHAM self-consistent solver in molrs
- [ ] Expose thin cpu:: WHAM API via molrs FFI (no C++ reimplementation of core math)
- [ ] Document units (kT or eV) and window input format (Doxygen)
- [ ] Add regression regressions/sampling/wham_double_well (hard-coded F(s) points)
- [ ] Verify multi-window umbrella samples reconstruct known double-well PMF shape
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：两窗合成高斯 → \(f_i\) 已知解析对照。
- **Regression**：双阱 + 多 `HarmonicBias` 窗短采样 + WHAM；硬编码若干 \(F(s)\)；≤ 60 s。
- 不在 tests 跑长 MD 收敛。

## Out of scope

- MBAR / umbrella integration 全套
- OPES
- 多维 WHAM（可先 1D CV 多窗）
- CUDA

---

**Librarian note：** math in molrs；engine streams；confidence medium (cross-repo molrs)。
