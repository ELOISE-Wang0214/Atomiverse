---
slug: concrete-pw-method
status: approved
created: 2026-07-23
track: C-architecture
grilled: false
---

# concrete-pw-method

## Summary

将平面波 DFT 的用户路径从 **`PWConfig` + `PWSolver` bag** 收敛到与高斯侧一致的 **concrete method**：参数为成员、`run(frame)` / `PWEvaluator` 消费方法对象、Frame 唯一 I/O。`PWConfig` 降为 private 或删除公共暴露；不改变物理（LDA PW 能量/力），只改 API 形状与债务边界。

## Domain basis

无新物理方程。架构约束见 notes：`[2026-07-17] SCF is a Driver; methods are concrete` 与 CLAUDE.md「No public `*Config` bags」。平面波自洽仍为内部循环（Driver 形后续可再迁）；本 spec 只做 **方法对象 + 成员参数** 对齐 `RHF` 用户故事。

## Design

- **类型**：`cpu::pw::PlaneWave` 或 `cpu::pw::PW`（命名与目录 `planewave` 对齐，避免 `PWSolver` 暗示“万能求解器 bag”）。
- **成员**：`ecut_ry`, `k_mesh`, `tol`, `max_iter`, `mix_beta`, `pulay_depth`, `upf_path` 等原 `PWConfig` 字段。
- **API**：`build(frame)` / `run(frame)` → 结果结构；`PWEvaluator` 持有/引用 method 或自身即薄包装。
- **迁移**：旧 `PWSolver{PWConfig}` 公共符号删除或 `#error` tombstone（与 `IntegralPack` 先例）；更新 examples/regressions。
- **禁止** `make_pw` / `pw_config()` factory。

### Reuse decision

- `reuse` `internal::` Poisson/Davidson/… 步 helpers 与 workspace。
- `reuse` `PWEvaluator` 作为 Component 入口（改构造）。
- `pattern` `cpu::scf::RHF` 成员参数 + `run(frame)`。
- `new` 公共方法类名；Config 内化。

## Files to create or modify

- `include/atomiverse/cpu/planewave/pw_method.h` (new) 或重写 `pw_solver.h`
- `include/atomiverse/cpu/planewave/pw_solver.h` — tombstone 或 private
- `include/atomiverse/cpu/planewave/pw_evaluator.h`
- `src/cpu/planewave/` 对应 TU
- `examples/` / `regressions/pw/` 用户路径
- `tests/cpu/planewave/` API 测试更新
- `docs/` 若有 PW 教程片段（最小 docs 任务可并入 Doxygen-only）

## Tasks

- [ ] Write failing unit tests for concrete PW method member API (tests/cpu/planewave/test_pw_method.cpp)
- [ ] Implement cpu::pw concrete method with member parameters and run(frame)
- [ ] Retarget PWEvaluator to method API; tombstone public PWConfig/PWSolver bag
- [ ] Update regressions/examples to Frame + method members only
- [ ] Document Doxygen public method surface
- [ ] Add or retarget regression regressions/pw (existing Si case) to new API with same hard-coded energy gate
- [ ] Run full check + test suite

## Testing strategy

- **Unit**：方法默认成员可构造；`run` 在 toy UPF 下有限能量。
- **Regression**：现有 Si/平面波回归改 API，**金标能量不变**。
- 无新物理 oracle。

## Out of scope

- 将 PW SCF 重写为 `Driver<pw::Phase>` 全阶段拆分（可 follow-up）
- USPP/PAW / 自旋 / hybrid PW
- CUDA API 对齐
- GPW/NDDO（并行 concrete-* specs）

---

**Librarian note：** `cpu/planewave/`；RHF pattern；confidence high。
