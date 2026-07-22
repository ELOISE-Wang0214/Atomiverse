---
spec: backend-mirror-07-docs
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 07-docs
depends_on: backend-mirror-06-gates
grilled: true
---

# backend-mirror-07-docs — 架构叙事与 harness 对齐

## Summary

重写 `CLAUDE.md` / `AGENTS.md` 的 Architecture 相关段落，以及 `.claude/notes/architecture.md`（或触发 `/mol:map` 刷新），使文档与已落地的布局一致：**pillar = frame 数据结构**；无 facade/shared/runtime/loop；完整词目录；`Driver`/`Criterion`；Backend 已拆；cpu/cuda 镜像；**明确无向后兼容**。删除过时的 CPUDriver/LoopHandler/kspace/三支柱旧说辞。

## Domain basis

无物理。

## Design

文档必须陈述：

1. 顶层仅 `core/` + `cpu/` + `cuda/`（+ `torch/` 实验说明）+ `molrs.hpp`  
2. Pillar：`frame/` only  
3. Infra 完整词：`communication`, `asynchronous`, `device`, `neighbor`, …  
4. Domain 完整词：`dynamics`, `montecarlo`, `minimize`, `planewave`, `observe`, `record`, …  
5. 合同在 `core/`；实现分后端  
6. 测试三层不变；架构 scrub 列表  
7. AGENTS 与 CLAUDE 不再互相矛盾（Criterion vs Handler）

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| CLAUDE managed architecture section | **generalize** rewrite to new truth |
| `.claude/notes/architecture.md` | **reuse** via map refresh after code lands |
| notes.md SCF driver notes | **reuse** — do not re-open SCF bag debate except path names if mentioned |

## Files to create or modify

- `CLAUDE.md` (Architecture, layout, key APIs)
- `AGENTS.md` (same; resolve Handler/Criterion drift)
- `.claude/notes/architecture.md` (refresh inventory)
- Optionally `.claude/notes/notes.md` one superseding layout note (only if required by /mol:note norms — prefer architecture sections)

## Tasks

- [ ] Rewrite CLAUDE.md Architecture + repository layout to the new tree and vocabulary
- [ ] Rewrite AGENTS.md Architecture to match CLAUDE (single Criterion story; Driver not CPUDriver)
- [ ] Refresh `.claude/notes/architecture.md` inventory to match disk after 01–06
- [ ] Grep docs/ and examples README for facade/CPUDriver/kspace/legacy short paths; fix narrative stragglers that claim old layout
- [ ] Add `tests/architecture/docs_vocab_scrub.sh` (or extend gates) asserting CLAUDE/AGENTS Architecture sections lack facade/CPUDriver/LoopHandler/kspace-as-canonical claims
- [ ] Run full check + test suite

## Testing strategy

- `type: docs` acceptance: banned legacy terms absent from CLAUDE/AGENTS architecture sections.
- Full suite still green.
- No science regressions required.

## Out of scope

- Code moves (01–06)
- Tutorial rewrite of every examples/N beyond broken includes (includes fixed in earlier phases)
- Torch redesign
