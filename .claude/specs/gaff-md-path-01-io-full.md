---
title: LAMMPS atom_style full 读入 — 电荷列与 GAFF 拓扑
slug: gaff-md-path-01-io-full
status: approved
chain: gaff-md-path
link: 01/07
created: 2026-07-21
---

# LAMMPS atom_style full 读入 — 电荷列与 GAFF 拓扑

## Summary

扩展 `cpu::load_lammps_data`，在既有 `atom_style=molecular` 路径之外正确识别并解析 `atom_style=full` 的 Atoms 行（`id mol-id type charge x y z`），把部分电荷写入 `kAtom_charge`，并在拓扑装载后调用 `build_special_list()`。完成后可直接读入 molpack-paper 的 `peo25.data`（GAFF2/AM1-BCC 单链）得到完整 Frame，而不再把 charge 误解析为坐标。

## Domain basis

- LAMMPS data 格式：`atom_style full` 行字段为 `atom-ID molecule-ID type q x y z`（可选镜像标志省略）；`molecular` 为 `atom-ID molecule-ID type x y z`。两者靠第 4 个浮点 token 之后是否仍有三个坐标浮点区分，**禁止**用 molrs 的 float-heuristic 自动猜 style（历史上会把 mol_id 当 type）。
- 电荷 \(q_i\) 单位为元电荷 \(e\)；与 `CoulLong` / `PairCoulCut` 一致，不在本环做单位换算。
- 拓扑：Bonds / Angles / Dihedrals 与 molecular 路径相同；`Frame::build_special_list()` 从 per-atom 键表生成 1–2 / 1–3 / 1–4 special 列表，供后续 special_bonds 缩放使用。
- 目标夹具：molpack-paper `figures/relaxation/chain_gaff/out/peo25.data`（177 atoms, 5 types, 有电荷）。本环测试可使用仓库内精简 fixture，不强制运行时依赖 paper 路径。

## Design

- **`load_lammps_data` 双路径**：在 Atoms 节解析前窥视首行 token 数 / 字段布局：
  - 6 字段且第 4–6 为坐标 → molecular（现有行为，0-index type，无 charge 写入或 charge=0）。
  - 7 字段 `id mol type q x y z` → full：写 `kAtom_charge`，type 仍 0-index。
  - 其它布局 → `std::runtime_error`，错误信息点明期望 full 或 molecular。
- **不要**经 molrs 自动猜 style；保持 C++ 自包含 parser（现文件注释已说明原因）。
- 拓扑装载结束后若存在 bonds，调用 `frame.build_special_list()`；无 bonds 时跳过。
- 头文件 `lammps_data.h` 更新 Doxygen：文档说明支持的 atom_style 与返回列。
- **tag / mol-id**：现有 molecular 路径用 row 当 tag；full 路径保持同一约定。molecule-id 若 Frame 尚无列则本环不新增 schema 列（Out of scope）；仅解析后丢弃 mol_id（与 molecular 现状一致）。

## Files to create or modify

- `include/atomiverse/cpu/io/lammps_data.h` — 文档更新
- `src/cpu/io/lammps_data.cpp` — full 路径 + special_list
- `tests/cpu/io/test_lammps_data_full.cpp` (new) — full 与 molecular 夹具
- `tests/cpu/io/fixtures/mini_full.data` (new) — 小 full 数据
- `tests/cpu/io/fixtures/mini_molecular.data` (new) — 回归 molecular
- `tests/CMakeLists.txt` — 注册测试

## Tasks

- [ ] Write failing tests in `tests/cpu/io/test_lammps_data_full.cpp`：full 夹具读出 nlocal、type、charge、坐标、键表；molecular 夹具仍绿；malformed 行抛异常
- [ ] Add fixtures `tests/cpu/io/fixtures/mini_full.data` 与 `mini_molecular.data`（≥2 原子、1 bond，full 含非零 q）
- [ ] Implement full Atoms 解析与 `kAtom_charge` 写入 in `src/cpu/io/lammps_data.cpp`；保留 molecular 路径
- [ ] Call `build_special_list()` after topology sections when bonds present
- [ ] Update Doxygen on `include/atomiverse/cpu/io/lammps_data.h`（@brief 写明 full/molecular）
- [ ] Register test in `tests/CMakeLists.txt`
- [ ] Run full check + test suite

## Testing strategy

- **快乐路径**：`mini_full.data` 读入后 charge 与文件逐原子相等（1e-12 相对）、type 0-index、bond/angle 计数正确；`nspecial` 在有键时非全零。
- **回归**：`mini_molecular.data` 与现有 PE 风格 6 列行为不变；既有依赖 `load_lammps_data` 的 amorphous/PE 路径不破。
- **边界**：缺 Masses、Atoms 列数错误、空文件 → 抛 `runtime_error`。
- **领域**：本环不比对库仑能量；只保证 \(q_i\) 与拓扑进入 Frame 的契约正确，供 02–05 使用。

## Out of scope

- `atom_style charge` / `atomic` / 倾斜盒 / image flags / Velocities 节。
- molecule-id 写入独立 Frame 列。
- 解析 LAMMPS `include` 力场文件或 `pair_coeff`（→ 04）。
- special_bonds 缩放实现（→ 02/03）。
- CUDA 侧 I/O。
