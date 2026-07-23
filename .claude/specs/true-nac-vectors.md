---
slug: true-nac-vectors
status: approved
created: 2026-07-23
track: A-response
grilled: false
---

# true-nac-vectors

## Summary

在闭壳层 RHF TDA 根之间提供 **真非绝热耦合向量** \(\mathbf{d}_{ij}=\langle\Psi_i|\nabla_R|\Psi_j\rangle\) [1/Å]（每原子 3 分量），替代/补全仅标量的 `tda_nac_strength`。为 fewest-switches surface hopping 解除“只有强度、没有方向”的阻塞；本 spec **不**实现 hopping 驱动本身。

## Domain basis

一阶导数耦合（NACV），对绝热态 \(i\neq j\)：

\[
\mathbf{d}_{ij}
= \frac{\langle\Psi_i|\nabla_R \hat{H}|\Psi_j\rangle}{E_j-E_i}
\quad (E_i\neq E_j)
\]

CIS/TDA 下可用过渡密度与基态/激发态力差形式构造（Chernyak–Mukamel；Send–Furche）。相位：固定规范使某一参考分量连续为正。参考：Send & Furche, *J. Chem. Phys.* **132**, 044107 (2010), DOI 10.1063/1.3292571。公开单位 1/Å（与现有 strength 文档一致）；内部可用 Bohr⁻¹。

## Design

- **API**：`RHF::tda_nac_vectors(frame, i, j) -> std::vector<double>` 长度 3N，或小型 `NacResult{vectors, energy_gap}` 结果结构（`tdhf.h` 风格，结果-only header）。
- **保留** `tda_nac_strength`：可定义为 \(\|\mathbf{d}_{ij}\|\) 或保持旧 FD 标量并文档化为 legacy；禁止静默改语义导致回归假绿。
- **依赖**：优先复用解析 TDA 力基础设施（若 `analytic-tda-forces` 已合入则 reuse；否则本 spec 内用 FD 力差形式实现 v1，并在 Tasks 标明）。
- **错误**：简并对、未收敛、\(i=j\) → throw。

### Reuse decision

- `reuse` `RHF::tda` / `tda_nac_strength` 落点与回归 `reg_h2_tda_nac`。
- `reuse` 或 `generalize` 力差外壳自 `tda_forces`（解析优先）。
- `pattern` `ExcitedStatesResult` / `MP2Result` 结果型。
- `new` `tda_nac_vectors` 成员 + 可选 `NacResult` — 无平行 System 类型。

## Files to create or modify

- `include/atomiverse/cpu/scf/rhf.h`
- `include/atomiverse/cpu/scf/tdhf.h`（若放 NacResult）
- `src/cpu/scf/` TDA/NAC 实现 TU
- `tests/cpu/scf/test_tda_nac.cpp` — 扩展向量 API
- `tests/cpu/scf/test_nac_vectors.cpp` (new)
- `regressions/scf/h2_nac_vectors/reg_h2_nac_vectors.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for NAC vector shape and antisymmetry (tests/cpu/scf/test_nac_vectors.cpp)
- [ ] Implement RHF::tda_nac_vectors (and NacResult if used) in cpu/scf
- [ ] Wire tda_nac_strength to ||d_ij|| or document legacy FD scalar without silent break
- [ ] Document Doxygen units [1/Å] and phase convention
- [ ] Add regression regressions/scf/h2_nac_vectors (hard-coded norm and/or component gate)
- [ ] Verify H2 two-root NAC norm consistent with literature or FD oracle
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：维度 3N；\(\mathbf{d}_{ji}\approx -\mathbf{d}_{ij}\)；能隙因子有限。
- **Unit**：\(\|\mathbf{d}\|\) 与 strength API 一致（若重定义为范数）。
- **Regression**：H2 两根；硬编码；≤ 60 s。
- 不实现表面跳跃动力学测试。

## Out of scope

- Fewest-switches / decoherence hopping Driver
- 开壳层 / TDDFT NAC
- 二阶导数耦合
- CUDA

---

**Librarian note：** `cpu/scf/`；depends logically on forces; confidence high。
