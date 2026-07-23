---
slug: concrete-nddo-method
status: approved
created: 2026-07-23
track: C-architecture
grilled: false
---

# concrete-nddo-method

## Summary

将半经验 NDDO 路径的剩余 **bag 形 `NDDOSolver`** 收敛为 **concrete methods**（`AM1` / `RM1` 已部分 evaluator 化）：公共用户故事为成员参数 + Frame，内部 workspace 私有；`NDDOSolver` 退出公共头或仅作 internal 基类。对齐 `RHF` 与 notes 债务表。

## Domain basis

无新参数化。AM1/RM1 方程与现实现一致；本 spec 纯 API/层边界。

## Design

- **公共**：`cpu::semiempirical::AM1` / `RM1` 方法型（或强化现有 `AM1Evaluator` 为唯一公共入口，solver 全 private）。
- **Internal**：`NDDOSolver` → `semiempirical/internal/`；公共 include 若仍见 `nddo_solver.h` 则 tombstone `#error`。
- **Evaluator**：`AM1Evaluator` 继续 `Semiempirical` kind；可持有 concrete method。
- **禁止** `NDDOSystem` 回流公共 API（已 tombstone 则保持）。

### Reuse decision

- `reuse` `AM1Evaluator` / `RM1Evaluator` / 参数表 / 力路径。
- `reuse` `reg_am1_water_*` 金标。
- `pattern` RHF concrete + PW/GPW concrete specs。
- `generalize` 若 AM1/RM1 仅差参数表 → 共享 internal NDDO engine（已存在则 reuse）。

## Files to create or modify

- `include/atomiverse/cpu/semiempirical/nddo_solver.h` — move/tombstone
- `include/atomiverse/cpu/semiempirical/am1_solver.h` / `am1_evaluator.h`
- `include/atomiverse/cpu/semiempirical/rm1_*.h`
- `src/cpu/semiempirical/` 实现路径调整
- `tests/cpu/semiempirical/` API 测试
- `regressions/semiempirical/` 若直接 include solver 则更新

## Tasks

- [ ] Write failing unit tests asserting public headers do not expose NDDOSolver bag (tests/cpu/semiempirical/test_nddo_api_surface.cpp)
- [ ] Move NDDOSolver to internal; concrete AM1/RM1 own parameters as members
- [ ] Retarget evaluators and regressions to non-bag public API
- [ ] Tombstone public nddo_solver.h if needed (#error)
- [ ] Document Doxygen AM1/RM1 public surface
- [ ] Verify am1_water_* regressions still pass
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：公共头 grep/隔离测试（可仿 scf isolation）。
- **Regression**：现有 AM1 能量/力/BOMD 门不变。

## Out of scope

- PM6/PM7 新参数
- DFTB/xTB（`tb-sk-product`）
- 解析 Hessian
- CUDA

---

**Librarian note：** `cpu/semiempirical/`；confidence high。
