---
spec: backend-mirror-06-gates
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 06-gates
depends_on: backend-mirror-05-dissolve
grilled: true
---

# backend-mirror-06-gates — 架构门禁与测试标准固化

## Summary

用 **可执行的 architecture scrub / 单测** 把新布局锁死，并完成任何仍指向旧布局的测试清理。门禁覆盖：无 facade、无顶层 shared math、无 cpu/cuda 路径漂移、无 CPUDriver/LoopHandler、无 runtime/loop/shared 目录、目录名单单词完整词表、测试路径镜像源码。**测试按新标准重写**，不为旧布局留豁免名单（除明确的 CPU-only 完整词目录白名单）。

## Domain basis

无物理。门禁是工程契约。

## Design

在 `tests/architecture/`（及必要时 `tests/core/`）落地脚本，建议集合：

| 门禁 | 断言 |
|------|------|
| `no_facade.sh` | 无 `include/atomiverse/facade` |
| `no_top_level_shared_math.sh` | 无顶层 `longrange/`、`potentials/` 模块 |
| `mirror_paths.sh` | cpu/cuda 共同域相对路径同名；禁止 `kspace`；禁止一侧 `md` 一侧 `dynamics` |
| `full_word_dirs.sh` | 禁用目录名单：`md,mc,comm,async,vib,pw,observer,recorder,minimizer,kspace` |
| `driver_vocab.sh` | 禁用 `CPUDriver,CUDADriver,LoopHandler,add_handler`；要求 `add_criterion` |
| `no_runtime_loop_shared_dirs.sh` | 无 `runtime/`,`loop/`,`shared/` 目录 |
| `tests_mirror_layout.sh` | `tests/cpu/<dir>` 仅在 `include/atomiverse/cpu/<dir>` 存在时出现（反之对已测模块） |

**CPU-only 白名单**（允许 cuda 缺失）：`montecarlo,vibration,semiempirical,gpw,io,bias,path,replica,task,linalg,evaluate`（以 03 落地为准，门禁配置化）。

既有 `no_factory_scrub.sh`、`stage_typing_scrub` **reuse**，仅更新路径。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| `tests/architecture/no_factory_scrub.sh` | **reuse** |
| stage typing scrub | **reuse** under tests/core |
| New scrubs | **new** — lock discuss decisions |

## Files to create or modify

- `tests/architecture/*.sh` (new scrubs)
- CMake/ctest registration
- Any remaining tests still on old paths (final sweep)
- Optional: `tests/architecture/cpu_only_dirs.txt` whitelist

## Tasks

- [ ] Write failing architecture scrubs listed in Design (RED)
- [ ] Implement scrub scripts + ctest registration so they run in default test config
- [ ] Final sweep: fix any remaining old-path or old-symbol test leftovers (rewrite, no skip)
- [ ] Ensure unit tests under tests/ still single-function and path-mirrored; move stragglers
- [ ] Add regression `regressions/layout/architecture_gates/` — documents that scrubs are part of ctest (or invokes a dry list); hard-coded expected gate names
- [ ] Run full check + test suite

## Testing strategy

- Scrubs are the primary tests; each fails closed on violation.
- No e2e science under tests/.
- Regression lists gate names as literals.

## Out of scope

- Implementing missing CUDA domains (stubs optional, not required beyond whitelist)
- CLAUDE/AGENTS full rewrite (07)
- Science feature work
