---
title: molrs CXX Frame ingress + Python handle exposure
status: approved
created: 2026-05-22
---

# molrs CXX Frame 入口 + Python 句柄暴露

## Summary

`rust-frame-lowering` 链第 1 步(共 4)。本链取代同名旧草案(typed `Slot<Domain,Tag,ValueSpec>` 需求系统 + Rust-first lowering 路线)——旧草案的需求-类型元编程整体删除,它解决的是一个不存在的问题。

收敛方向:`molrs.Frame` 是唯一规范主机端交换 Frame;执行 Frame 按 backend 划分(`cpu`/`cuda`/`torch` 各一,所有 driver 共用);每 `run()` 物化一次、不进 step 循环;Frame 只装系统模型 + 通用结果。

本 spec 的目标:让 `molrs.Frame` 能被 Atomiverse C++ 通过 CXX 读写。`molrs-cxxapi/src/bridge.rs` 现在只有输出(`AtvMolRec` staging + `xyz_write`/`molrec_write_zarr`),**没有 `Frame` 入口**。本 spec 在 bridge 新增 `Frame` 的零拷贝读接口 + 建/写接口,并让 `molrs-python` 暴露底层 `Frame` 句柄,供 `-03` 把同一个 `Frame` 同时交给 pyo3 和 cxx。

**本 spec 全部改动在 molrs 仓库(`/Users/roykid/work/molcrafts/molrs`),不改 Atomiverse。**

## Domain basis

无新物理、无方程。`molrs::Frame`(`molrs-core/src/frame.rs`)是 `HashMap<String,Block>` + `meta` + `simbox` 的通用 block 字典,每个 `Block` 含异构 typed column(f32/f64/i64/bool)。单位由上层约定(eV/Å/fs/amu/K),bridge 只做无语义的数据搬运。

## Design

**读接口(零拷贝)。** cxx 支持返回 `rust::Slice<const T>`。新增按 dtype 的零拷贝列读取:

- `frame_block_names(frame) -> Vec<String>`、`frame_has_block(frame, name) -> bool`
- `frame_block_columns(frame, block) -> Vec<String>`、`frame_block_nrows(frame, block) -> i64`
- `frame_column_f64(frame, block, col) -> &[f64]`、`frame_column_i64(...) -> &[i64]` 等
- `frame_simbox(frame) -> [f64; 9]`(按现有 `SimBox` 形状)

**建 / 写接口。** `frame_new() -> Box<Frame>`;`frame_set_column_f64(frame, block, col, &[f64])` 等按 dtype 的写入(create-or-update:目标 `Block`/column 不存在则创建,存在则就地更新);`frame_set_simbox(...)`。这套接口同时服务生产侧写回与测试侧构造。

**Python 句柄暴露。** `molrs-python` 的 `Frame` pyo3 类新增一个方法,返回承载底层 `*mut Frame` 的 PyCapsule。`-03` 的 Atomiverse cxx bridge 据此借用同一个 `Frame`。两侧 `Frame` 同源于 `molrs-core`——`molrs-python` 与 `molrs-cxxapi` 必须锁同一 `molrs-core` 版本。

## Files to create or modify

- `molrs-cxxapi/src/bridge.rs` — 新增 `Frame` 读 + 建/写的 `extern "Rust"` 声明
- `molrs-cxxapi/src/lib.rs` — 实现上述函数
- `molrs-python/src/frame.rs` — 暴露底层 `Frame` 句柄(PyCapsule)
- molrs 仓库相应 Rust 测试

## Tasks

- [ ] Add a zero-copy Frame column-read interface to the molrs-cxxapi bridge
- [ ] Add a Frame create / write interface (frame_new + create-or-update set_column) to the bridge
- [ ] Expose the underlying Frame handle as a PyCapsule from molrs-python
- [ ] Add Rust tests covering the bridge Frame read/write round-trip
- [ ] Run cargo test in the molrs workspace

## Testing strategy

- molrs 仓库 `cargo test`:`frame_new` 构造 `Frame`,经 `frame_set_column_*` 写入各 dtype 列,经 `frame_column_*` 读回,逐元素一致;`frame_simbox` 往返一致。
- PyCapsule:从 `molrs-python` 取句柄,确认指向的 `Frame` 与 Python 侧是同一实例(改一侧另一侧可见)。

## Out of scope

- Atomiverse 侧的物化器(见 `-02`)。
- `run()` 签名 / Python 路径改动(见 `-03`)。
- cuda/torch Frame(见 `-04`)。
