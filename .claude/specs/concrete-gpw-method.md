---
slug: concrete-gpw-method
status: approved
created: 2026-07-23
track: C-architecture
grilled: false
---

# concrete-gpw-method

## Summary

将 GPW DFT 用户路径从 **`GpwConfig` + `GPWSolver` bag** 收敛为 **concrete method**（成员参数 + Frame `run` / `GPWEvaluator`），与 `RHF` / 本战役 `concrete-pw-method` 同构。不改 Γ-point 闭壳层物理与解析力；删除公共 Config bag 暴露。

## Domain basis

无新方程。架构债务：notes 与 CLAUDE.md 禁止公共 `*Config` / bag solver。GPW 现有 `GthPotential` / `GaussianBasis` lookup 构造保持（那是表行查找，不是 factory world-model）。

## Design

- **类型**：`cpu::gpw::GPW`（或 `GpwDft`）方法对象；成员吸收原 config（cutoff、相对网格、functional、`local_pp_on_grid` 等）。
- **Evaluator**：`GPWEvaluator` 绑定 method；力开关为成员。
- **Tombstone**：公共 `GPWSolver`/`GpwConfig` 删除或 `#error`。
- **Kind 表 / static data**：不变。

### Reuse decision

- `reuse` `GPWEvaluator`、内部 SCF/collocate、DFT-D3 表。
- `reuse` `GthPotential` / `GaussianBasis` lookup 构造。
- `pattern` `concrete-pw-method` + `RHF`。
- `new` 公共方法类；Config 内化。

## Files to create or modify

- `include/atomiverse/cpu/gpw/gpw_method.h` (new) 或改写 `gpw_solver.h`
- `include/atomiverse/cpu/gpw/gpw_solver.h` — tombstone
- `include/atomiverse/cpu/gpw/gpw_evaluator.h`
- `src/cpu/gpw/` 实现
- `regressions/` GPW 案例 API 更新
- `tests/cpu/gpw/test_gpw_method.cpp` (new)
- examples 若引用 GPWSolver

## Tasks

- [ ] Write failing unit tests for concrete GPW method API (tests/cpu/gpw/test_gpw_method.cpp)
- [ ] Implement concrete GPW method with member parameters
- [ ] Retarget GPWEvaluator; tombstone public GpwConfig/GPWSolver bag
- [ ] Update GPW regressions/examples to method API
- [ ] Document Doxygen public surface
- [ ] Verify existing GPW energy/force regression gates on new API
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：方法可 `run` toy 体系。
- **Regression**：既有 GPW 能量/力门改 API，金标不变。

## Out of scope

- k-points / spin / hybrid GPW
- Driver-phase 全拆
- PW/NDDO（并行 specs）
- CUDA

---

**Librarian note：** `cpu/gpw/`；confidence high。
