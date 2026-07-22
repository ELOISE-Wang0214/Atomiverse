---
spec: backend-mirror-02-absorb
status: approved
created: 2026-07-22
chain: backend-mirror
phase: 02-absorb
depends_on: backend-mirror-01-core
grilled: true
---

# backend-mirror-02-absorb — 删除顶层 shared 算子，并入各后端

## Summary

删除顶层 `include/atomiverse/longrange/` 与 `include/atomiverse/potentials/`（以及对应 `src/longrange/`、跨后端 `tests/longrange/`），把其中数学 **吸收进各后端自己的实现树**：CPU 侧进入 `cpu/longrange/`；CUDA 侧 pair kernel 进入 `cuda/potentials/`，并把 `cuda/kspace/` **改名为** `cuda/longrange/`（同相对路径镜像）。**无 shared/、无顶层算子库、无转发头、不向后兼容。** 测试按新路径重写。

## Domain basis

无新物理公式；Ewald/PME/FMM/pair kernel 数值行为本环默认 **bit-level 或既有容差保持**（搬家不改算法）。若移动迫使 include 图变化，以现有 unit 金标为准。

## Design

**禁止再出现：**

```text
include/atomiverse/longrange/
include/atomiverse/potentials/
src/longrange/
tests/longrange/          # 跨后端共享测试树
```

**吸收映射：**

| 现状 | 目标 |
|------|------|
| `include/atomiverse/longrange/*` | `include/atomiverse/cpu/longrange/`（可 `internal/` 子目录）+ 若 CUDA 需要则 **复制/改写** 到 `cuda/longrange/`，禁止 CPU include CUDA 或反向 |
| `src/longrange/*` | `src/cpu/longrange/`（及 cuda 对等若有） |
| `tests/longrange/*` | `tests/cpu/longrange/*`（CUDA 有则 `tests/cuda/longrange/*`） |
| `include/atomiverse/potentials/pair_kernel.h` | `include/atomiverse/cuda/potentials/pair_kernel.cuh`（或 `.h` 若刻意 CUDA-free bridge；但 **不得** 留在顶层）。CPU pair 继续只用 `cpu/potentials/` 自有实现 |
| `include/atomiverse/cuda/kspace/` | **rename** → `include/atomiverse/cuda/longrange/`；`src/cuda/kspace/` → `src/cuda/longrange/`；`tests/cuda` 中 kspace 引用同步 |
| `include/atomiverse/cpu/longrange/` | 保留目录名 `longrange`（完整词），吸收顶层数学 |

**Torch：** 若 `torch/` 依赖顶层 `potentials/`，本环必须改为依赖 `cuda/potentials` 或 torch 自有拷贝，或暂时从构建中排除该依赖并在 Out of scope 说明——**不得** 为 torch 保留顶层 shared。

### Reuse decision

| Candidate | Decision |
|-----------|----------|
| top-level `longrange/*` operators | **generalize（absorb）** into `cpu/longrange/` (+ cuda copy only if needed) |
| `potentials/pair_kernel.h` | **generalize（absorb）** into `cuda/potentials/` |
| `cpu/longrange/{ewald,pme,pppm}` components | **reuse** — keep as backend-owned components; retarget includes |
| `cuda/kspace/*` | **reuse + rename** to `cuda/longrange/` |
| CFMM in `cpu/scf/internal` | **reuse** — update includes to cpu/longrange after absorb |

## Files to create or modify

- Move/absorb: `include/atomiverse/longrange/*` → `include/atomiverse/cpu/longrange/` (and cuda if required)
- Move: `src/longrange/*` → `src/cpu/longrange/`
- Move: `tests/longrange/*` → `tests/cpu/longrange/*`
- Move: `include/atomiverse/potentials/pair_kernel.h` → under `cuda/potentials/`
- Rename tree: `include/atomiverse/cuda/kspace/` → `longrange/`; `src/cuda/kspace/` → `longrange/`
- Update all consumers (cpu longrange components, scf cfmm, cuda pair_lj, cmake)
- Delete empty top-level longrange/potentials trees
- CMake source lists for longrange/pair

## Tasks

- [ ] Write failing scrub under `tests/architecture/` that forbids `include/atomiverse/longrange` and `include/atomiverse/potentials` as top-level module dirs and forbids `cuda/kspace` path
- [ ] Absorb top-level longrange headers/sources into `cpu/longrange/` (and cuda longrange only if a cuda TU needs them); update includes
- [ ] Move `pair_kernel` into `cuda/potentials/`; rewrite cuda pair consumers; ensure CPU does not depend on it
- [ ] Rename `cuda/kspace` → `cuda/longrange` (include + src + tests); zero `kspace` path leftovers in product trees
- [ ] Rewrite `tests/longrange/` into `tests/cpu/longrange/` (and cuda mirror if needed); delete cross-backend `tests/longrange/`
- [ ] Fix torch/bindings/cmake consumers of deleted top-level paths (no shims)
- [ ] Add regression `regressions/layout/backend_owned_longrange/` — CPU public header include path under `atomiverse/cpu/longrange/` only
- [ ] Run full check + test suite

## Testing strategy

- Scrub: no top-level shared math dirs; no `cuda/kspace`.
- Existing FMM/PME unit numerical tests move under `tests/cpu/longrange/` and stay within prior tolerances.
- Regression documents backend-owned include root.

## Out of scope

- Domain renames (`md`→`dynamics`, …)
- Driver/Backend class redesign
- New PME algorithms or FMM features
- Full CLAUDE architecture prose rewrite
