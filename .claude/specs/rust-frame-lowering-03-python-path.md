---
title: run() takes molrs.Frame; remove _convert.py and the CpuFrame binding
status: approved
created: 2026-05-22
revised: 2026-05-22
follows: rust-frame-lowering-02-cpu-materialize
---

# run() 接收 molrs.Frame;移除 _convert.py 与 CpuFrame 绑定

## Summary

`rust-frame-lowering` 链第 3 步,完成 CPU「入口」。把 driver `run()` 改为接收 `molrs.Frame`:在 cxx 侧借用 Python 传来的同一个 `FrameRef` 句柄(`-01` 的 PyCapsule),内部经 `-02` 的 `cpu::materialize` / `write_back` 完成物化与写回。Atomiverse 停止绑定 Frame 类型本身——删 `_convert.py`、移除暴露给 Python 的 `CpuFrame`。无向后兼容:旧路径在本 spec 内一并删除。

## Domain basis

无新物理、无方程。AIMD 约束(`MD` + `SCFEvaluator` 共用同一个 `cpu::Frame`,见 `aimd.md`)在本 spec 必须保持可跑通。

## Design

**句柄穿透。** Python 用户持有 `molrs.Frame`(molrs-python,内部即 `molrs_ffi::FrameRef`)。调用 Atomiverse 的 `MD.run(frame)` 时,Python 绑定层取 `-01` 暴露的 PyCapsule(承载 `*mut molrs_ffi::FrameRef`),把 `FrameRef` 句柄经 cxx 交给 C++;C++ 据此重建一个指向同一 `Store` 的 `FrameRef`、借用其中同一个 `Frame` 跑完 `run()`——零拷贝,不 copy-out / copy-in。借用期严格限定在单次 `run()` 内。

**run() 签名。** 内部 `run()`(`facade/driver.h`)现接收 `cpu::Frame&`。面向用户的入口改为接收 `molrs::FrameRef`:`run()` 入口 `cpu::materialize` 成 `cpu::Frame`、执行、出口 `cpu::write_back`。`cpu::Frame` 降级为纯执行产物,不再是用户接触的类型——更新其 file/class 文档(现 docstring 称其为「primary data container users interact with」,需改写)。

**移除 Frame 绑定。** 删 `python/src/atomiverse/_convert.py`;移除绑定层暴露的 `CpuFrame`(`_cpp.CpuFrame` 及其 `set_*` / `get_*`)。`MD.run()` / `RHF()` 等 Python 调用点改为直接传 `molrs.Frame`。Atomiverse 仍绑定 driver / component API(`MD`、`set_pair`、`LJ`、`SCFConfig` …)——只是不再绑定 Frame 数据类型。

## Files to create or modify

- `include/atomiverse/facade/driver.h`(及 cpu driver / `driver_impl.h`)— 面向用户的 `run()` 入口接收 `molrs::FrameRef`
- Atomiverse cxx 侧 `FrameRef` 句柄重建(借用 `molrs-python` 的 PyCapsule;经 `molrs.hpp`)
- `python/src/atomiverse/_convert.py` — 删除
- Atomiverse Python 绑定层(`_binding` / bind 源)— 移除 `CpuFrame`;`run()` 走 `molrs.Frame`
- `MD.run()` / `RHF()` 等 Python 调用点
- `tests/` — Python 集成测试:从 `molrs.Frame` 起跑 MD

## Tasks

- [ ] Wire the molrs.Frame FrameRef PyCapsule handle through the cxx bridge into the CPU driver
- [ ] Change the user-facing run() entry to accept a molrs.Frame and materialize / write-back internally
- [ ] Remove the CpuFrame type from the Atomiverse Python binding surface
- [ ] Delete python/src/atomiverse/_convert.py and update MD.run() / RHF() call sites
- [ ] Add a Python integration test running MD from a molrs.Frame end to end
- [ ] Run full check + test suite

## Testing strategy

- **Python 集成**:构造 `molrs.Frame` → `MD` 跑若干步 → 结果写回同一 `molrs.Frame`,能量有限。
- **AIMD 不回归**:`MD` + `SCFEvaluator` 共用一个物化出来的 `cpu::Frame`,H2 STO-3G NVE 满足 `aimd` 能量守恒判据。
- **无 CpuFrame**:Python 侧不再能 import / 构造 `CpuFrame`;`_convert.py` 不存在。
- **回归**:完整 `ctest` 全绿。

## Out of scope

- cuda/torch Frame 清理与物化器(见 `-04`)。
- Spec/Plan/Driver 拆分(链外,单独 spec)。
