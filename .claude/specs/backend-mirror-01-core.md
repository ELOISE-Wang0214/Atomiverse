---
spec: backend-mirror-01-core
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 01-core
grilled: true
---

# backend-mirror-01-core — 合同层迁入 core，删除 facade

## Summary

把现有顶层 `include/atomiverse/facade/` 中与设备无关的合同类型（Plan、Stage、Role、Verdict、IterationStatus、molrs 版本契约，以及无后端算法的 Driver 模板/实现辅助）**整体迁入** `include/atomiverse/core/`，删除 `facade/` 目录；全仓 include 与测试按新路径重写。本环 **不做向后兼容**：无 `#include` 转发头、无旧路径别名、无 deprecated 宏。

## Domain basis

无新物理。本环只动类型布局与构建图。

## Design

**目标布局（本环交付后）：**

```text
include/atomiverse/core/
  …既有 precision/units/schema/binding…
  plan.h                 # 自 facade/plan.h
  stage.h                # 自 facade/stage.h
  role.h                 # 自 facade/component_role.h（文件名单词/清晰）
  verdict.h              # 自 facade/loop.h 中 Verdict/Outcome/ConvergenceMetrics
  iteration_status.h     # 自 facade/
  molrs_contract.h       # 自 facade/
  driver_template.h      # 自 facade/driver.h（DeviceTag 模板调度原语；非后端实现）
  driver_template_impl.h # 自 facade/driver_impl.h
```

**删除：** `include/atomiverse/facade/` 整树。

**不在本环：** 目录英文重命名（`md`→`dynamics` 等）、Backend 拆解、顶层 `longrange/`/`potentials/` 吸收、`CPUDriver` 改名（见后续子 spec）。

**规则：**

- 所有 `#include "atomiverse/facade/..."` 改为 `atomiverse/core/...`；grep 门禁在本环末锁定零 `facade` 引用。
- `tests/facade/` 迁到 `tests/core/`（`stage_typing_scrub.sh`、`test_plan.cu` 等），内容按新路径改写，不保留 `tests/facade` 空壳。
- CMake / 文档中指向 facade 的路径同步更新（本环至少改到能通过配置与相关单测；叙事全文可在链末 docs 子 spec 收口，但 **不得** 留下可编译的 facade 路径）。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| `facade/plan.h` Plan/PlanBuilder | **reuse** — 迁 `core/plan.h`，语义不变 |
| `facade/stage.h` Phase enums + slot_of | **reuse** — 迁 `core/stage.h` |
| `facade/loop.h` Verdict/Outcome/metrics | **reuse** — 迁 `core/verdict.h`（或保留逻辑文件名但禁止再叫 facade） |
| `facade/component_role.h` | **reuse** — 迁 `core/role.h` |
| `facade/iteration_status.h`, `molrs_contract.h` | **reuse** — 迁 `core/` |
| `facade/driver.h` + `driver_impl.h` | **reuse** — 迁 `core/driver_template.h`（+ impl）；后端 `cpu::CPUDriver` 仍 include 新路径 |
| `core/binding.h` | **reuse** — 已在 core，不动 |

## Files to create or modify

- `include/atomiverse/core/plan.h` (new or moved)
- `include/atomiverse/core/stage.h` (new or moved)
- `include/atomiverse/core/role.h` (new or moved)
- `include/atomiverse/core/verdict.h` (new or moved)
- `include/atomiverse/core/iteration_status.h` (new or moved)
- `include/atomiverse/core/molrs_contract.h` (new or moved)
- `include/atomiverse/core/driver_template.h` (new or moved)
- `include/atomiverse/core/driver_template_impl.h` (new or moved)
- Delete tree: `include/atomiverse/facade/`
- All in-tree includes of `atomiverse/facade/*` (src/, include/, tests/, examples/, regressions/, python/, app/, tools/ as present)
- `tests/core/stage_typing_scrub.sh` (moved/rewritten from `tests/facade/`)
- `tests/core/test_plan.cu` or host equivalent (moved/rewritten)
- CMake test registration paths for former `tests/facade` targets

## Tasks

- [ ] Write failing unit/scrub expectations under `tests/core/` that assert public includes resolve via `atomiverse/core/{plan,stage,role,verdict}.h` and that `include/atomiverse/facade` is absent
- [ ] Move facade contract headers into `include/atomiverse/core/` with the target filenames; update internal includes among them
- [ ] Rewrite every repository `#include "atomiverse/facade/..."` to `atomiverse/core/...` (no shim headers)
- [ ] Delete `include/atomiverse/facade/` entirely
- [ ] Rewrite `tests/facade/*` into `tests/core/*` (stage_typing_scrub, plan tests); remove `tests/facade/`
- [ ] Wire CMake/ctest entries to `tests/core/` paths; drop facade targets
- [ ] Add regression example `regressions/layout/core_contract_includes/` — minimal public consumer that includes only `core/plan.h` + `core/stage.h` and links (hard-coded: compiles; documents new include root)
- [ ] Run full check + test suite

## Testing strategy

- **Unit/scrub:** `tests/core/stage_typing_scrub.sh` continues to fail the build on hand-written Phase `static_cast<int>` but greps new paths; plus a scrub that `rg facade` under include/src/tests is empty (except this spec text if any — code trees only).
- **Happy path:** existing plan/stage unit coverage still green after move.
- **No e2e under tests/.**
- **Regression:** `regressions/layout/core_contract_includes/` — out-of-library include of core contracts only.

## Out of scope

- Renaming `CPUDriver` / `CUDADriver` / Backend dissolve
- Top-level `longrange/` / `potentials/` absorption
- Domain directory renames (`md`→`dynamics`, …)
- CLAUDE.md / AGENTS.md 全文重写（链末 docs；但本环不得依赖 facade 路径仍存在）
- Torch backend redesign
- Electronic-structure API debt (`SCFSolver` bags)
