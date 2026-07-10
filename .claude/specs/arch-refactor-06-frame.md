---
title: arch-refactor-06-frame — 删除 cpu::Frame 死拓扑 shim + 钦定 write_back 契约 + Frame 文档修正
status: approved
date: 2026-07-10
chain: arch-refactor (06/09)
depends_on: [arch-refactor-05-context]
---

# arch-refactor-06-frame — 删除 cpu::Frame 死拓扑 shim + 钦定 write_back 契约 + Frame 文档修正

## Summary

本篇是 arch-refactor 链（06/09）的 Frame 收尾：（1）删除 `cpu::Frame` 中嵌套的
`Bonds/Angles/Dihedrals` flat-xtensor 死存储 shim（include/atomiverse/cpu/frame/frame.h:82-104），
将仍靠它编译的 legacy 测试 / benchmark / 调试工具 fixture 全部迁移到
`frame.atoms.add_bond / add_angle / add_dihedral` per-atom API，并修正同文件自相矛盾的头注释
（:26-28 声称子结构 "are deleted"，:8 又说保留为编译 shim）；（2）钦定并以 golden round-trip 测试
固化 `cpu::write_back` 的字段契约——**image 与 tag 必须回写**（当前静默丢失是缺陷），
type/mass/topology/box 钦定为不回写并文档化，scalars 保持只写通道，per-frame metadata 继续走
`write_frame_xyz_meta` 既有路径；`cuda::write_back`（= `to_cpu` + `cpu::write_back`）自动继承该契约；
（3）修正三处与实现不符的注释/文档：molrs.hpp 的 bridge_f32/f64 configure 选择说法（实际单一
always-f64 bridge）、CUDA 双拓扑现实（单卡 flat list + 多卡 per-atom SoA）、以及散落各处引用
`frame.bonds` 的陈旧文档。材料化/回写是冷路径（每 `run()` 一次），但 frame.h 被热路径 include，
删 shim 后须以 PE stage-5 全量基准确认无扰动（链级门禁）。本篇逻辑独立，排在 05 之后仅为避免
头文件合并冲突；与 `.claude/specs/precision-migration.md` 正交——只文档化现状精度 tier，不改 tier。

## Domain basis

本篇为基础设施重构，不引入新物理，但触及两个有物理语义的字段约定，
均沿用 LAMMPS 的 Atom 布局与 image-flag 约定
（Thompson et al., *Comput. Phys. Commun.* **271** (2022) 108171,
DOI: 10.1016/j.cpc.2021.108171）：

- **image flags**（无量纲整数计数，`image(i,d)` = 原子 i 沿 d 方向穿越盒面的净次数）：
  非包裹坐标重建 `x_unwrap = x_wrapped + image · h`（h 为 cell 矩阵，单位 Å）。
  MD 中随步演化，write_back 丢失它即丢失扩散/MSD 类分析所需的信息——这是本篇修复的缺陷。
- **tag**（全局原子恒等 ID）：多卡域分解迁移后本地行序可变，tag 是跨设备/跨 rank 的唯一
  身份锚；write_back 不回写 tag 则 molrs 侧列与 cpu 侧行的对应关系在迁移后失真。
- per-atom 拓扑布局（`num_bond[i]` / `bond_atom(i,m)`=partner 全局 tag / `bond_type(i,m)`，
  angle 由中心原子持有、dihedral 由第 2 原子持有，`nspecial/special` 为 1-2/1-3/1-4 排除表）
  即 LAMMPS `Atom` 等价物，见 include/atomiverse/cpu/frame/atom.h:198-246（运行时自适应
  stride :53-62，lazy 分配 :258-267）。

单位约定不变：Å、fs、amu、eV、eV/Å（core/units.h）。

## Design

### 1. 删除死拓扑 shim（frame.h）

`cpu::Frame` 的嵌套 `struct Bonds/Angles/Dihedrals`（frame.h:82-104）是 "Transitional
dead-storage shims"（:68-80）：无任何 force kernel / halo / migration / NL 路径读取
（bonded 力全部迭代 `frame.atoms.num_bond[i]` 等 per-atom 布局），仅剩 fixture 侧写入。
全仓核实的存活消费者（比注释声称的 "~30" 少，多数早已迁移）：

- tests/cpu/test_md_no_alloc.cpp:107-113 — `frame.bonds` 填充（死存储，BondHarmonic 实际读
  per-atom，故该 no-alloc 测试目前**并未真正执行键力路径**；迁移后 bonded 路径被真实纳入
  no-alloc 门禁，断言本身不变）。
- benchmarks/md/lj/validation/bench_nve_limit_lammps_align.cpp:213-217, 273-277, 324-328 —
  三个 case 填充死存储 → 零键力 → `NveLimitLammpsAlign_VelocityVerletConservesEnergy` 的
  非零初始能量守卫必然失败（benchmark 默认 OFF 所以潜伏）。**迁移即修复**，恢复其设计意图。
- tools/debug_amorphous/lammps_loader.h:107-210（写入）+ run_atomiverse.cpp:245-249
  （`frame.bonds.size()` 日志读取，改为 per-atom 计数求和）。该工具由
  `ATV_BUILD_DEBUG_AMORPHOUS`（CMakeLists.txt:757-759）门控，须一并迁移以保持可编译。

迁移范式对照既有样板 tests/cpu/md/test_pe_single_chain.cpp:99-111：
`frame.atoms.map_clear()` → `map_set(0, nlocal)`（add_* 经 tag→local map 解析属主行）→
`add_bond/add_angle/add_dihedral(全局 tag…, type)` → `build_special_list()`。

同时重写 frame.h 头注释：删除 :4-11 的 "Transitional state" 段与 :68-80 的 shim 说明，
保留并收敛 :26-31 的 "子结构已删除、数据在 per-atom 布局" 陈述，消除自相矛盾。

### 2. write_back 契约（钦定，写入 materialize.h 文档与 golden test）

现状（src/cpu/frame/materialize.cpp:214-255）：只写 x/y/z/vx..vz/fx..fz/charge（:222-232）
+ scalars/stress/stress_kinetic 块（:240-254）；**静默漏 tag/image**；materialize 侧
（:77-212）从不读 scalars 块。钦定契约如下：

| 字段 | materialize（molrs→cpu） | write_back（cpu→molrs） | 依据 |
|---|---|---|---|
| x/y/z | 必需，缺失抛异常 | 写 | 现状保持 |
| vx/vy/vz, fx/fy/fz | 可选，缺省零 | 写 | 现状保持 |
| charge | 可选，缺省零 | 写 | 现状保持 |
| image（ix/iy/iz 三列） | 可选，缺省零（i32→int8） | **新增：写**（int8→i32 加宽，与读取格式对称） | MD 中演化，丢失是缺陷 |
| tag | 可选，缺省恒等（i32→TagInt） | **新增：写**（TagInt→i32 受检窄化，溢出抛 `std::runtime_error`，不静默截断） | 多卡迁移后行序可变 |
| type | 必需 | 不写，文档化 | per-run 不可变输入 |
| mass | 可选，缺省零 | 不写，文档化 | per-run 不可变输入 |
| topology（bonds/angles/dihedrals 块） | 可选 → per-atom scatter | 不写，文档化 | per-run 静态 |
| box | 可选 | 不写，文档化 | per-run 静态（钦定） |
| scalars / stress / stress_kinetic | **不读**（只写通道，文档化不对称性） | 写 | 现状保持 |
| per-frame metadata（step/energy） | — | 不并入；继续走 `write_frame_xyz_meta`（src/cpu/io/frame_writer.cpp:34） | 既有路径 |

**CUDA 继承**：`cuda::write_back`（src/cuda/frame/materialize.cu:21-26）= `to_cpu` 到瞬态
host frame + `cpu::write_back`，自动继承契约。已核对 `cuda::to_cpu`
（src/cuda/frame/transfer.cu:252-297）回传 x/v/f/tag/charge/image/scalars/box——覆盖契约全部
"写" 字段；瞬态 host frame 的 type/mass 为零值，但契约不回写二者，故无垃圾流出 molrs 侧。
此一致性论证写入 include/atomiverse/cuda/frame/materialize.cuh 注释。

**Golden test**（tests/cpu/test_materialize.cpp 扩展）：materialize → 修改契约内字段
（x/v/f/charge/image/tag 置换/scalars）**且**修改契约外字段（type/mass）→ write_back →
re-read molrs 列：契约内字段呈现修改值（image 经 ix/iy/iz 往返、tag 呈现置换），
契约外的 type/mass/topology/box 列保持 materialize 前原值；另加 tag > INT32_MAX
（bigbig 默认 TagInt=int64）时 write_back 抛异常的边界 case。

### 3. 注释 / 文档修正

- **include/atomiverse/molrs.hpp:5-6**：声称 "precision-selected at CMake configure time
  (bridge_f32 or bridge_f64)"——实际 CMakeLists.txt:658-659 注明 molrs-core 硬编码 f64、
  单一 always-double bridge、FFI 边界零拷贝。改写注释与之一致。
- **CUDA 双拓扑现实**：include/atomiverse/cuda/frame/frame.cuh:8-13 声称字段集与
  cpu::Frame/torch::Frame "完全一致"——CPU shim 删除后 CPU 不再有 flat 拓扑块而 CUDA 有。
  在 frame.cuh 头注释写清双拓扑：单卡路径用 flat `BondList/AngleList/DihedralList`
  （:34-58, :148-151），多卡域分解用 per-atom SoA（:96-129）；CUDA kernel 尚未迁移到
  per-atom 拓扑，故 `cpu::to_cuda` 每次 H2D 重新展平（transfer.cu:70-78 已有局部注释，
  提升到文件头注释）。此注释为 09-harness 修正 CLAUDE.md "三 Frame 字段集相同" 表述提供依据
  （本篇不改 CLAUDE.md）。
- **陈旧 `frame.bonds` 引用**（实现均已读 per-atom，纯注释/文档漂移）：
  include/atomiverse/cpu/comm/halo_exchange.h:25,81,131、src/cpu/comm/halo_exchange.cpp:206、
  src/cpu/comm/domain.cpp:370、include/atomiverse/cpu/potentials/pair_lj_cut.h:26、
  docs/developer/cpu-internals.md:301、docs/tutorial/ch10-classical-md.md:247。
  统一改述为 per-atom 布局（`frame.atoms.num_bond / bond_atom / …`）。

### 生命周期 / 所有权

materialize/write_back 是冷路径（每 `run()` 一次，性能不敏感）；但 frame.h 被热路径
include——删 shim 仅移除无人读取的成员、无运行时行为变化，按链级门禁以 PE stage-5 全量基准
（`polyethylene_tension`，对照 benchmarks/md/pe/performance/runs/ 最新 baseline）确认无扰动。
本篇不新增公共符号；`cpu::write_back` 签名不变，仅扩大写出字段集。

## Files to create or modify

- include/atomiverse/cpu/frame/frame.h
- tests/cpu/test_md_no_alloc.cpp
- benchmarks/md/lj/validation/bench_nve_limit_lammps_align.cpp
- tools/debug_amorphous/lammps_loader.h
- tools/debug_amorphous/run_atomiverse.cpp
- include/atomiverse/cpu/frame/materialize.h
- src/cpu/frame/materialize.cpp
- tests/cpu/test_materialize.cpp
- include/atomiverse/cpu/comm/halo_exchange.h
- src/cpu/comm/halo_exchange.cpp
- src/cpu/comm/domain.cpp
- include/atomiverse/cpu/potentials/pair_lj_cut.h
- docs/developer/cpu-internals.md
- docs/tutorial/ch10-classical-md.md
- include/atomiverse/molrs.hpp
- include/atomiverse/cuda/frame/frame.cuh
- src/cuda/frame/transfer.cu
- include/atomiverse/cuda/frame/materialize.cuh

## Tasks

- [ ] Write failing golden round-trip tests for the write_back contract in
      tests/cpu/test_materialize.cpp: materialize → mutate x/v/f/charge/image/tag(permuted)/
      scalars **and** type/mass → write_back → re-read; assert contract fields carry the
      mutations (image via ix/iy/iz, tag via permutation) while type/mass/topology/box
      molrs columns stay untouched; add a tag > INT32_MAX overflow case expecting
      std::runtime_error.
- [ ] Implement image (int8→i32 widen, ix/iy/iz columns) and tag (TagInt→i32 checked
      narrowing, throw on overflow) write-back in src/cpu/frame/materialize.cpp, and
      document the full field contract table (incl. type/mass/topology/box designated
      non-write-back, scalars write-only asymmetry, metadata via write_frame_xyz_meta) in
      include/atomiverse/cpu/frame/materialize.h per doxygen style with units.
- [ ] Migrate all shim-compiling fixtures to the per-atom topology API
      (map_clear → map_set → add_bond/add_angle/add_dihedral → build_special_list, per the
      tests/cpu/md/test_pe_single_chain.cpp pattern): tests/cpu/test_md_no_alloc.cpp,
      benchmarks/md/lj/validation/bench_nve_limit_lammps_align.cpp (fixes the latent
      zero-bond-force cases), tools/debug_amorphous/lammps_loader.h, and replace the
      frame.bonds.size() logging in tools/debug_amorphous/run_atomiverse.cpp with per-atom
      count sums; keep every existing assertion unchanged.
- [ ] Delete the nested Bonds/Angles/Dihedrals sub-structs from
      include/atomiverse/cpu/frame/frame.h and rewrite its header comment to remove the
      "compile shim" / "are deleted" contradiction.
- [ ] Correct stale frame.bonds/angles/dihedrals doc references to the per-atom layout in
      include/atomiverse/cpu/comm/halo_exchange.h, src/cpu/comm/halo_exchange.cpp,
      src/cpu/comm/domain.cpp, include/atomiverse/cpu/potentials/pair_lj_cut.h,
      docs/developer/cpu-internals.md, and docs/tutorial/ch10-classical-md.md.
- [ ] Fix the include/atomiverse/molrs.hpp header comment: single always-f64 bridge
      (molrs-core hardcodes f64, zero-copy FFI), no configure-time bridge_f32/bridge_f64
      selection, matching CMakeLists.txt:658-659.
- [ ] Document the CUDA dual-topology reality (single-GPU flat BondList/AngleList/
      DihedralList + multi-GPU per-atom SoA; kernels not yet ported, per-transfer
      re-flatten in cpu::to_cuda) in include/atomiverse/cuda/frame/frame.cuh and
      src/cuda/frame/transfer.cu header comments, and note the write_back contract
      inheritance (to_cpu field coverage vs contract) in
      include/atomiverse/cuda/frame/materialize.cuh.
- [ ] Verify against the PE stage-5 performance benchmark: run polyethylene_tension after
      the frame.h change and confirm stage-5 throughput is within the historical noise
      band of the latest baseline under benchmarks/md/pe/performance/runs/.
- [ ] Run full check + test suite (cmake --build build && ctest --test-dir build
      --output-on-failure).

## Testing strategy

- **Happy path**：全量 `ctest`（默认 ~16 s 门禁）绿色；`test_cpu_materialize` 既有 4 个
  case（round-trip / fail-loud / optional-zero / runnable）在契约扩展后不回归。
- **Golden round-trip（新）**：materialize → 修改 → write_back → re-read，逐字段按契约表断言
  （写字段呈现修改；type/mass/topology/box 列保持原值；scalars 只写不读的不对称性以
  "materialize 不消费 scalars 块" 断言固化）。
- **边界**：tag 超出 i32 范围抛异常且异常信息点名字段；空 topology / 零原子外的最小 frame
  write_back 不抛。
- **迁移行为**：`test_cpu_md_no_alloc` 断言不变（bonded 路径迁移后真实进入 no-alloc 门禁）；
  `bench_nve_limit_lammps_align` 全部 case 通过——其中能量守恒 case（域验证，
  science.required）迁移后须真实产生键力并满足既有 1e-4 守恒界，即恢复 LAMMPS 对齐设计意图。
- **CPU↔CUDA 回归**：GPU 节点重跑 `test_cuda_materialize`（tests/cuda/test_materialize.cu，
  已覆盖 tag/image/charge 对齐与 D2H write_back）与 `tests/cuda/test_frame.cu`、
  `tests/cuda/md/test_multi_gpu.cu`（flat 拓扑传输 + 多卡 per-atom 路径），确认契约继承无破坏。
- **域验证 / 性能门禁**：PE stage-5（`polyethylene_tension`）吞吐对照
  benchmarks/md/pe/performance/runs/ 最新 baseline，无扰动（噪声带内，±5%）——frame.h 为
  热路径 include，此为链级不变量要求的全量基准确认。

## Out of scope

- CUDA kernel 迁移到 per-atom 拓扑（`cpu::to_cuda` 每次 H2D 重展平的现状保持；仅文档化）。
- 精度 tier 变更——与 `.claude/specs/precision-migration.md` 正交，本篇只修注释不动 tier。
- box / topology 的回写支持（钦定 per-run 静态；NPT/Deform 的盒演化经由 cuda::to_cpu 回到
  执行 frame，但不进入 molrs.Frame 回写契约）。
- `torch::Frame` 的拓扑表示与其 `bonds` 成员（torch 后端自有实存储，非本 shim）。
- image 的 int8 存储位宽（±127 穿越上限）扩展。
- CLAUDE.md "三 Frame 字段集相同" 表述修正与 `/mol:map` 刷新——09-harness 的职责，本篇仅在
  frame.cuh 中提供事实依据。
- `write_frame_xyz_meta` 元数据路径的任何改动。
