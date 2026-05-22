---
title: cuda/torch Frame cleanup, materializers, and three-Frame alignment
status: code-complete
created: 2026-05-22
follows: rust-frame-lowering-03-python-path
---

# cuda/torch Frame 清理、物化器与三 Frame 对齐

## Summary

`rust-frame-lowering` 链第 4 步(末)。把 `cuda::Frame` 里违反「Frame 只装系统模型」原则的 `SCFWorkspace`、`PairList` 拆出去,实现 cuda/torch 物化器,并让三个执行 Frame 结构对齐。

## Domain basis

无新物理、无方程。CLAUDE.md 规则:「热路径数据结构只装基本物理量——绝不装预计算索引、缓存或派生映射」。SCF 工作区(Fock/ERI/本征向量)与邻居表是算出来的、方法相关的 scratch,不属于 Frame。

## Design

**移除 cuda::Frame 的 scratch。** `cuda::Frame` 现内嵌 `SCFWorkspace`(含 n⁴ 的 `ERI`、`H`、`P`、`C`、本征求解 scratch、基组数据)和 `PairList`:

- `SCFWorkspace` → 由 CUDA SCF 组件持有(`cpu` 侧 SCF 工作区本就在 `SCFSolver`,对齐之)。
- `PairList` → 由 CUDA `NeighborList` 组件持有(`cpu` 侧本就如此——`cpu::Frame` 不含 pair list)。

**三 Frame 对齐。** `cuda::Frame`、`torch::Frame` 内嵌 `core/system_scalars.h` 的 `SystemScalars`(`-02` 引入);补 `tag`、`image`、`conserved`,与 `cpu::Frame` 对称。修订后三个执行 Frame 的字段集只差「设备表示」(xtensor / 裸设备指针 / `at::Tensor`),不差「装哪些物理量」。

**物化器。** 新建 `cuda/materialize.{cuh,cu}`、`torch/materialize.*`(各置 backend 根,与 `cpu/materialize` 平行)。`molrs.Frame`(主机)→ 设备 `Frame`:物化时显式 H2D;`write_back` 时显式 D2H。结构与 `cpu/materialize` 对称——按已知 block/column 名读写。

## Files to create or modify

- `include/atomiverse/cuda/frame.cuh` — 移除 `SCFWorkspace`/`PairList`;加 `tag`/`image`/`conserved`;内嵌 `SystemScalars`
- `include/atomiverse/torch/frame.cuh` — 同上对齐
- `include/atomiverse/cuda/scf/*` — 接管原 `SCFWorkspace`
- `include/atomiverse/cuda/md/neighbor_list*`(及相关 `.cu`)— 接管原 `PairList`
- `include/atomiverse/cuda/materialize.cuh` + `src/cuda/materialize.cu` — 新建
- `include/atomiverse/torch/materialize.*` — 新建
- `CLAUDE.md` / `.claude/notes/notes.md` — 更新 Frame / 数据驻留描述
- `tests/cuda/`、`tests/torch/` — 物化往返测试

## Tasks

- [x] Remove SCFWorkspace from cuda::Frame and relocate it to CUDA SCF components
- [x] Remove PairList from cuda::Frame and relocate it to CUDA NeighborList
- [x] Align cuda::Frame and torch::Frame with cpu::Frame (SystemScalars, tag/image/conserved)
- [x] Implement cuda::materialize and torch::materialize with explicit H2D / D2H
- [x] Update CLAUDE.md and notes for the finalized Frame model
- [x] Run full check + test suite
- [x] Hygiene: pure relocation refactor — no legacy/shims/dead code introduced
      (transfer.cu's now-unused copy_*_prec helpers deleted in the rewrite)

## Testing strategy

- **结构断言**:grep `include/atomiverse/cuda/frame.cuh` 无 `SCFWorkspace`/`ERI`/`PairList` 符号。
- **对齐**:`cuda::Frame`/`torch::Frame` 含 `tag`/`image`/`conserved` 且内嵌 `SystemScalars`。
- **物化往返**(CUDA toolkit 在场时):`molrs.Frame` → `cuda::Frame` → 写回逐元素一致;torch 同理。
- **回归**:完整 `ctest` 全绿;CUDA 测试在 toolkit 在场时通过。

## Out of scope

- 按需物化优化(链外)。
- CUDA driver taxonomy 迁移(`cuda-driver-taxonomy.md`)。
