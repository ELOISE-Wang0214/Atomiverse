---
title: cpu::Frame materializer, charge field, shared SystemScalars
status: approved
created: 2026-05-22
follows: rust-frame-lowering-01-molrs-bridge
---

# cpu::Frame 物化器 + charge 字段 + 共享 SystemScalars

## Summary

`rust-frame-lowering` 链第 2 步。在 `-01` 把 molrs CXX bridge 的 `Frame` 入口做好之后,本 spec 在 Atomiverse 侧实现 `molrs.Frame` ↔ `cpu::Frame` 的物化与写回,补齐 `cpu::Frame` 的 `charge` 字段,并抽出三 backend 共享的 `SystemScalars`。

收敛设计:`molrs.Frame` 是唯一规范交换层;执行 Frame 按 backend 划分;物化每 `run()` 一次、不进 step 循环;Frame 只装系统模型 + 通用结果,算出来的方法相关量留在组件里。

## Domain basis

无新物理、无方程。单位 eV/Å/fs/amu/K(`core/units.h`)。

Frame 承载的物理量:每原子量(位置、速度、力、质量、类型、tag、image、电荷)、拓扑(bonds/angles/dihedrals)、`box`、通用系统标量(`energy`、`kinetic`、`conserved`、`stress`、`stress_kinetic`)。

唯一领域约束来自 AIMD(`aimd.md`):Born-Oppenheimer MD 是 `cpu::md::MD` driver + `cpu::scf::SCFEvaluator` 组件跑在**同一个 `cpu::Frame`** 上。这条约束决定执行 Frame 按 backend 而非 MD/SCF domain 划分。

## Design

**charge 字段。** `cpu::Frame::Atoms` 新增 `charge`(`xt::xtensor<double,1>`,单位 e)。`Capability::Charge` 已在 `facade/component_role.h`,但 `cpu::Frame` 至今无对应字段——补齐,使能力位与字段一一对应。

**共享 SystemScalars。** 新建 `core/system_scalars.h`,定义 `SystemScalars { energy, kinetic, conserved, stress[6], stress_kinetic[6] }`,保留现有 Doxygen 注释(虚功/动能应力的时间切片语义)。`cpu::Frame` 内嵌之,取代散落的标量成员。`cuda`/`torch` 在 `-04` 跟进。

**物化器。** 新建 `cpu/materialize.{h,cpp}`,置于 cpu backend 根(与 `transfer.h` 平级——属共享基础设施,符合 CLAUDE.md backend 分层规则)。

- `cpu::Frame cpu::materialize(const molrs::FrameRef&)`:按已知 block/column 名读 `molrs.Frame` 填充——`atoms.{x,y,z,vx,vy,vz,fx,fy,fz,mass,type,tag,image,charge}`、`bonds.{atomi,atomj,type}`、`angles.*`、`dihedrals.*`、`box`。缺失的可选字段(速度、电荷)留零;缺失的必需字段(位置、类型)抛异常并指明字段名(沿用 `_convert.py` 的「失败响亮」风格)。
- `void cpu::write_back(const cpu::Frame&, molrs::FrameRef&)`:把位置、速度、力、电荷、系统标量写回同一个 `molrs.Frame`。
- v1 物化固定全字段集;「只物化 Plan 用到的字段」是后续优化,链外。

`molrs::FrameRef` 是 Atomiverse 侧对 `-01` bridge 句柄的薄封装(经 `molrs.hpp`)。

## Files to create or modify

- `include/atomiverse/cpu/frame.h` — 加 `charge`;内嵌 `core/system_scalars.h` 的 `SystemScalars`
- `include/atomiverse/core/system_scalars.h` — 新建
- `include/atomiverse/cpu/materialize.h` + `src/cpu/materialize.cpp` — 新建
- `tests/cpu/test_materialize.cpp` — 新建
- `tests/CMakeLists.txt` — 注册新测试

## Tasks

- [ ] Add a charge field to cpu::Frame::Atoms with Doxygen units
- [ ] Extract core/system_scalars.h and embed SystemScalars in cpu::Frame
- [ ] Implement cpu::materialize (molrs.Frame -> cpu::Frame) over known block/column names
- [ ] Implement cpu::write_back (cpu::Frame -> molrs.Frame)
- [ ] Write round-trip materialization tests and register them in tests/CMakeLists.txt
- [ ] Run full check + test suite

## Testing strategy

- **往返**:`molrs.Frame` → `cpu::materialize` → `cpu::write_back` → 同一 `molrs.Frame`,位置/速度/力/电荷/拓扑/box/标量逐元素一致(dtype cast 容差内)。
- **缺字段**:缺 `atoms.x` 抛异常并指明字段名;缺速度则物化出的 `cpu::Frame` 速度为零、不抛。
- **直接起跑**:用物化出来的 `cpu::Frame` 经 `cpu::md::MD` 跑数步,能量有限。
- **回归**:完整 `ctest` 全绿。

## Out of scope

- `run()` 签名 / Python 路径(见 `-03`)。
- cuda/torch Frame 与物化器(见 `-04`)。
- 按需物化优化(链外)。
