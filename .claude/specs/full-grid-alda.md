---
slug: full-grid-alda
status: done
created: 2026-07-23
track: A-response
grilled: false
---

# full-grid-alda

## Summary

把 `RKS::tda` 的 ALDA 交换核从 **Dirac 对角近似** 升到 **完整格子上的 MO 乘积积分** \(f^{\mathrm{xc}}\)，使 LSDA TDA 激发能相对对角近似有可测位移，并在公开 API 上通过 `ExcitedStatesResult::alda == true` 与回归门钉住。仍限于闭壳层 LSDA；杂化 TDDFT 与解析力不在本 spec。

## Domain basis

Tamm–Dancoff TDDFT 的 A 矩阵（闭壳层单重态，原子单位）：

\[
A_{ia,jb} = (\varepsilon_a-\varepsilon_i)\delta_{ij}\delta_{ab} + 2(ia|jb) + (ia|f^{\mathrm{xc}}|jb)
\]

ALDA 核（LSDA，\(n_\uparrow=n_\downarrow=n/2\)）在实空间对角：

\[
(ia|f^{\mathrm{xc}}|jb) = \int \phi_i(\mathbf{r})\phi_a(\mathbf{r})\,
f^{\mathrm{xc}}_{\mathrm{ALDA}}(n(\mathbf{r}))\,
\phi_j(\mathbf{r})\phi_b(\mathbf{r})\,d\mathbf{r}
\]

\(f^{\mathrm{xc}} = \delta^2 E_{\mathrm{xc}} / \delta n^2\) 由现有 Padé LDA / Dirac 交换参数化的二阶导数给出。参考：Casida 1995；Dreuw & Head-Gordon, *Chem. Rev.* **105**, 4009 (2005), DOI 10.1021/cr0505627。单位：内部 Hartree；公开激发能 eV（与 `tdhf.h` 一致）。

## Design

- **实体**：`cpu::scf::RKS::tda`、`internal` TDA/ALDA 路径、`ExcitedStatesResult`（已有 `alda` 标志）。
- **行为**：LSDA 参考下，对每个 \((ia,jb)\) 通道在 KS 格子上累加 MO 乘积 × \(f^{\mathrm{xc}}\)，取代当前仅对角 \(ia=jb\) 的 Dirac 位移。
- **生命周期**：一次 `run()` 后 `tda(n)` 只读收敛密度 / 轨道 / 格子；不分配在 `launch` 热路径上的新全局状态。
- **错误**：非 LSDA 参考保持现有行为（`alda=false` 或显式 throw，与 `test_alda_diff` 一致）；无静默回退到对角近似。

### Reuse decision

- `reuse` `RKS::tda` / `ExcitedStatesResult`（`include/atomiverse/cpu/scf/rks.h`, `tdhf.h`）— 扩展实现，不改用户入口签名。
- `reuse` 现有 KS 格子与 LDA 参数化（`cpu/scf` grid + `cpu/xc`）— 同一套 `f_xc` 导数源。
- `reuse` 对角 ALDA 回归 `reg_h2_rks_alda_tda` / `reg_h2_alda_vs_cis` — 收紧期望：完整 ALDA 与纯 CIS@RKS 仍可区分，且相对对角近似 \(\Delta\omega\) 非零（钉硬编码）。
- `pattern` `RHF::tda` A 矩阵组装 — 耦合项插入点对齐 HF 的 \((ia|jb)\) 块布局。
- `new` 无 — 不新开 Driver / 不新 Phase。

## Files to create or modify

- `include/atomiverse/cpu/scf/rks.h`（文档：完整 ALDA）
- `src/cpu/scf/rks.cpp`（或现有 TDA 实现 TU）— MO 乘积 ALDA 积分
- `include/atomiverse/cpu/scf/internal/` 下现有 ALDA/TDA helper（若已抽出）或同文件内联至第二调用点再抽
- `tests/cpu/scf/test_alda_diff.cpp` — 扩展/替换对角-only 断言
- `tests/cpu/scf/test_full_grid_alda.cpp` (new) — 单函数：格子积分对称性 / 对角极限
- `regressions/scf/h2_full_grid_alda/reg_h2_full_grid_alda.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for full-grid ALDA kernel (tests/cpu/scf/test_full_grid_alda.cpp → TestFullGridAlda)
- [ ] Write failing extension of tests/cpu/scf/test_alda_diff.cpp (grid vs diagonal shift)
- [ ] Implement full MO-product ALDA contribution on RKS LSDA TDA path in src/cpu/scf/
- [ ] Document RKS::tda ALDA semantics (Doxygen, units eV) in include/atomiverse/cpu/scf/rks.h
- [ ] Add regression example regressions/scf/h2_full_grid_alda/reg_h2_full_grid_alda.cpp (public API; hard-coded ω [eV])
- [ ] Verify against H2 LSDA textbook window and non-zero Δ vs diagonal-only path
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：`test_full_grid_alda` — 单函数核：对称性 \(K_{ia,jb}=K_{jb,ia}\)、单点密度下对角块与解析 Dirac 极限一致（容差 1e-8 Hartree）。
- **Unit**：`test_alda_diff` — LSDA 完整 ALDA 与对角近似激发能差 > 阈值阈值。
- **Regression**：`regressions/scf/h2_full_grid_alda` — H2 / LSDA / 小基；硬编码 S1 [eV]（PySCF 或本引擎对角路径对照，注释工具/版本/日期）；`alda==true`；≤ 60 s。
- **非 e2e under tests/**。

## Out of scope

- 杂化 TDDFT（B3LYP/PBE0 的 HF 交换核在 TDA 中）
- 完整线性响应 TDDFT（含 B 矩阵）
- 解析 TDA / TDDFT 力与 NAC 向量（见 `analytic-tda-forces`, `true-nac-vectors`）
- CUDA / 开壳层 UKS TDA
- GGA kernel（PBE 二阶核）

---

**Librarian note：** blueprint 2026-07-23；placement `cpu/scf/`；confidence high。
