---
spec: backend-mirror-03-paths
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 03-paths
depends_on: backend-mirror-02-absorb
grilled: true
---

# backend-mirror-03-paths — 完整英文目录重命名（infra + domain）

## Summary

在 `cpu/` 与 `cuda/`（及镜像的 `src/`、`tests/`）上，把缩写/旧目录名 **重命名为完整英文单词**，两边相对路径严格相同；**无兼容别名、无旧路径转发、测试整树按新路径重写**。本环只改路径与 include/namespace 对齐的必要符号，不改算法语义。

## Domain basis

无新物理。

## Design

### 强制重命名表（include + src + tests 同步）

| 旧目录 | 新目录（完整词） | 两端 |
|--------|------------------|------|
| `comm/` | `communication/` | cpu, cuda |
| `async/` | `asynchronous/` | cpu, cuda |
| `md/` | `dynamics/` | cpu, cuda |
| `mc/` | `montecarlo/` | cpu（cuda 若尚无实现：建同名 stub 头 **或** 仅在有代码时建——优先 **有实现才建目录**；镜像门禁白名单记录 CPU-only 完整词目录） |
| `minimizer/` | `minimize/` | cpu, cuda |
| `vib/` | `vibration/` | cpu |
| `pw/` | `planewave/` | cpu, cuda |
| `observer/` | `observe/` | cpu, cuda |
| `recorder/` | `record/` | cpu, cuda |
| `semiempirical/` | `semiempirical/` | 已是完整词，保持 |
| `longrange/` | `longrange/` | 保持（02 已统一） |
| `potentials/` | `potentials/` | 保持完整词 |
| `neighbor/` `frame/` `driver/` `component/` `device/` `replica/` `task/` `bias/` `path/` `linalg/` `scf/` `gpw/` `xc/` `io/` | 保持 | |

**CPU-only 完整词目录**（允许 cuda 暂缺，但不得在 cuda 侧用别名）：  
`montecarlo/`, `vibration/`, `semiempirical/`, `gpw/`, `io/`, `bias/`, `path/`, `replica/`, `task/`, `linalg/`（若 cuda 无则白名单）。

**命名空间：** 新代码与本环触及的符号应向 `ATV::cpu::communication`、`ATV::cpu::dynamics` 等与目录一致的后缀收敛。允许 **本环内** 为减少爆炸用 `namespace dynamics = md` 过渡？——**否：用户要求无向后兼容**。直接改命名空间；全仓替换。

**预设文件随域走：**

- `md/md.h` → `dynamics/md.h` 或 `dynamics/dynamics.h`（类名可暂留 `MD` 作为缩写类型名；目录必须是 `dynamics`）
- `minimizer/minimize.h` → `minimize/minimize.h`
- `mc/*.h` → `montecarlo/`
- Evaluate / PathDriver：若仍在 `driver/`，本环 **移到** `evaluate/` 与 `path/`（完整词目录），避免 driver 杂物间。

**类名 `CPUDriver` 不在本环改**（04）。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| All domain/infra implementations under old dirs | **reuse** content, **rename** paths only |
| Closest pattern `frame/` mirror | **pattern** — same relative path on cpu/cuda |
| Existing tests under `tests/cpu/md` etc. | **rewrite** to `tests/cpu/dynamics` etc.; delete old paths |

## Files to create or modify

- Rename trees under `include/atomiverse/cpu|cuda/`, `src/cpu|cuda/`, `tests/cpu|cuda/` per table
- Create `include/atomiverse/cpu/evaluate/`, `path/` as needed when moving presets out of driver
- Update every include, CMakeLists, python bindings, examples, regressions paths
- Architecture notes only if required for build; full docs in 07

## Tasks

- [ ] Write failing path-mirror scrub listing required renames (old names must not exist; new names must)
- [ ] Rename infra dirs `comm→communication`, `async→asynchronous` on cpu and cuda (include/src/tests); fix includes/namespaces
- [ ] Rename domain dirs per table (`md→dynamics`, `mc→montecarlo`, `minimizer→minimize`, `vib→vibration`, `pw→planewave`, `observer→observe`, `recorder→record`); fix includes/namespaces
- [ ] Move Evaluate/single_point and PathDriver out of `driver/` into `evaluate/` and `path/`; update consumers
- [ ] Rewrite tests tree paths to mirror new layout; delete old test directories
- [ ] Update examples, regressions, python bindings, cmake — no old-path shims
- [ ] Add regression `regressions/layout/path_words/` — includes `atomiverse/cpu/dynamics/` and `atomiverse/cpu/communication/` headers only under new paths
- [ ] Run full check + test suite

## Testing strategy

- Scrub forbids old directory names under include/atomiverse/cpu and cuda.
- All moved unit tests keep prior assertions; only paths/namespaces change.
- Regression locks new include roots.

## Out of scope

- CPUDriver/CUDADriver → Driver rename
- Backend dissolve
- Criterion vs LoopHandler unify (04)
- Top-level shared (done in 02)
- Full AGENTS/CLAUDE rewrite (07)
