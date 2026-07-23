---
slug: nph-ensemble
status: approved
created: 2026-07-23
track: A-cell-md
grilled: false
---

# nph-ensemble

## Summary

提供 **NPH（等焓-等压）** 可达路径：恒温器与气压计解耦，使 `MD` 可在无原子热浴的情况下演化体积/胞，公开验证扩展焓 \(H = E + PV\)（及 barostat 动能项若适用）守恒窗口。解除 roadmap「NPH unreachable because thermostat fused into MTK」债务。

## Domain basis

NPH：固定 \(N,P,H\)。实现上采用 **无 Nosé–Hoover 原子链的 barostat 积分**（对角 MTK barostat 或 PR 胞动力学的 thermostat-free 变体）。守恒量为扩展焓（含气压计动能时为完整扩展哈密顿量）。参考：Martyna et al., *J. Chem. Phys.* **101**, 4177 (1994)。单位：eV、Å、fs。

## Design

- **类型**：`cpu::md::NphIntegrator` 或 `NoseHooverNPT` 的显式模式 `thermostat=false` **仅当** 不破坏现有 NPT 默认——优先 **独立 Integrator 类**，避免布尔上帝参数。
- **API**：成员 `P_target` [pressure unit 与 MTK 一致]、`tau_P`、可选胞质量；**无** \(T\) 目标。
- **依赖**：`parrinello-rahman-npt` 若已合入可 reuse 胞更新；否则 v1 仅 **各向同性** 体积自由度即可关闭 NPH 空洞（文档标明；各向异/NPH-PR 可 follow-up）。
- **Criterion**：不声明收敛；步数预算终止。

### Reuse decision

- `reuse` 压力张量 / virial 路径同 NPT。
- `reuse` `Box` 体积缩放原语。
- `pattern` `NoseHooverNPT` phase 绑定与 Trotter 切片（去掉原子 NHC）。
- `new` 独立 `NphIntegrator`（或等价命名）— 不扩 `*Config`。

## Files to create or modify

- `include/atomiverse/cpu/dynamics/nph_integrator.h` (new)
- `src/cpu/dynamics/nph_integrator.cpp` (new)
- `tests/cpu/dynamics/test_nph_integrator.cpp` (new)
- `regressions/md/nph_lj/reg_nph_lj.cpp` (new) + `REFERENCE.md` (new)

## Tasks

- [ ] Write failing unit tests for NPH enthalpy conservation helper (tests/cpu/dynamics/test_nph_integrator.cpp)
- [ ] Implement thermostat-free NPH Integrator (isotropic v1 minimum)
- [ ] Document ensemble and conserved H in Doxygen
- [ ] Add regression regressions/md/nph_lj (hard-coded relative drift of H)
- [ ] Verify short LJ NPH run keeps H within gate
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：单步体积变化后 \(H\) 定义自洽。
- **Regression**：LJ；相对漂移 \(\Delta H / |H|\) 硬编码上限；≤ 60 s。

## Out of scope

- 完整生物膜表面张力协议
- CUDA
- 与 thermostat 同时打开的「半 NPH」模糊模式

---

**Librarian note：** `cpu/dynamics/`；confidence medium-high（v1 isotropic enough）。
