---
title: Multi-GPU bonded — per-atom topology on device + tag→local map + re-materialization
status: in-progress
created: 2026-06-06
chain: multigpu-bonded
link: 2
---

# Multi-GPU bonded — per-atom topology on device + tag→local map + re-materialization

## Summary
链条 `multigpu-bonded` 第 2 环(承接已落地的 `-01-reverse-halo`)。让键合拓扑在多卡域分解下"活着":
当前 `cuda::Frame` 只有**扁平** tuple 列表(`BondList.atomi/atomj` 等,局域索引,单卡 tag==row),没有逐原子拓扑;
`migrate` 不搬运键合数据;ghost 只带 x/y/z/tag/type。于是分区/迁移后,跨卡的键对无法解析。本环加三块基础设施:
(1) **逐原子键合 SoA** 上 GPU(`num_bond`/`bond_atom`=全局 tag/`bond_type`,angle/dihedral/special 同),分区时按
**owner 过滤**从 CPU frame materialize,partner 存**全局 tag**;(2) `migrate` 搬运这些逐原子数组(拓扑随 owner 走);
(3) 每个 reneighbor 在 halo 之后,用**设备 tag→local 映射**(owned+ghost 排序 tag + 二分)把逐原子拓扑**重新展开**成
扁平 tuple 列表(`bonds.atomi/atomj/...`),partner tag 解析到本卡 local(owned 或 ghost)行;解析失败的 tuple **丢弃并计数**
(halo 过薄诊断)。本环**不**算键力(留给 -03/-04)——只保证拓扑正确就位、可解析。

## Domain basis
无新物理。可证伪基础:重新展开出的扁平 tuple 集合,在 2/4 卡上(对每个 owner 卡)必须与单卡 materialize 的 tuple 集合
**逐元素等价**(同一组 (owner_tag, partner_tag, type),partner 解析到正确的 local 行,位置一致);键合 owner 规则
(bond=min-tag、angle=center、dihedral=2nd atom)与单卡/CPU 一致,保证每条 tuple 恰由一个 owner 卡持有(不重不漏)。

## Design
- **逐原子键合 SoA**(加到 `cuda::Frame::Atoms`,设备指针,容量 = `capacity_`):
  `int* num_bond`、`std::int64_t* bond_atom`(`capacity*kMaxBond`,全局 tag)、`int* bond_type`;
  `int* num_angle`、`std::int64_t* angle_atom`(`capacity*2*kMaxAngle`)、`int* angle_type`;
  `int* num_dihedral`、`std::int64_t* dihedral_atom`(`capacity*3*kMaxDihedral`)、`int* dihedral_type`;
  `int* nspecial`(`capacity*3`)、`std::int64_t* special`(`capacity*kMaxSpecial`)。caps 复用 cpu/atom.h
  (`kMaxBond=6/Angle=6/Dihedral=8/Special=16`)。**单卡不受影响**:这些数组仅多卡路径分配+使用(单卡走现有扁平 materialize)。
- **materialize at partition**(`multi_gpu.cu::partition` + helper):每张卡对其 owned 原子,从全局 CPU frame 拷入
  owner 过滤的逐原子拓扑(partner 存全局 tag,不预解析)。`migrate` 之后由迁移搬运维持。
- **migrate carries topology**:扩展 `MigBuf` + `gather_state`/scatter,把上述逐原子数组随原子搬到新 owner 卡。
- **device tag→local map**(新 `tag_map.cu`):每个 reneighbor 在 `build_halo` 之后,对 `[0, n_local+n_ghost)` 的
  `tag` 排序得 `(tag, local_idx)`,kernel 用二分 `resolve(tag)→local`(越界/缺失返回 -1)。
- **re-materialize**(新 `bonded_materialize.cu`):kernel 遍历 owned 原子,对每条 bond/angle/dihedral,把 partner 全局 tag
  经 tag→local 解析为 local 行,写入扁平 `bonds.atomi/atomj/type`(等);任一 partner 解析为 -1 → 丢弃该 tuple 并
  `atomicAdd` drop 计数器(teardown 时按 -01 的 DropTally 风格告警)。`n_bonds_/n_angles_/n_dihedrals_` 更新为本卡解析成功数。
- 调用点:`MultiGpu::reneighbor()` 末尾(halo 重建后)调用 build-tag-map + re-materialize;`partition` 末尾首次建立。

## Files to create or modify
- `include/atomiverse/cuda/frame.cuh` — `Atoms` 增逐原子键合 SoA 指针。
- `src/cuda/frame.cu`(或 materialize/destroy 所在)— 分配/释放新数组(按 capacity)。
- `src/cuda/transfer.cu` — 单卡 materialize 保持现状;新增"逐原子拓扑拷贝"helper 供多卡 partition 用。
- `src/cuda/md/multi_gpu.cu` — partition materialize 逐原子拓扑;`migrate` 搬运;`reneighbor` 调 tag-map + re-materialize。
- `include/atomiverse/cuda/md/multi_gpu.cuh` — 声明。
- `src/cuda/md/tag_map.cu` (new) — 设备 tag→local(排序+二分)。
- `src/cuda/md/bonded_materialize.cu` (new) — 重新展开扁平 tuple + drop 计数。
- `tests/cuda/md/test_cuda_multi_gpu_bonded.cu` (new) — tag map 解析 + 2-GPU 重展开等价单测(GH200)。

## Tasks
- [x] **P2.1** Add per-atom bonded SoA to `cuda::Frame::Atoms` + allocate/free by capacity (single-GPU path unaffected) — GH200 build green, LJ 13/13
- [x] **P2.2** Materialize owner-filtered per-atom topology (partner=global tag) on each subdomain at `partition` (`materialize_per_atom_bonded`, gated on `has_bonded`) — GH200 build green, LJ 13/13 (fill correctness validated once observable via re-materialize, P2.5)
- [ ] **P2.3** Extend `migrate` (`MigBuf`/`gather_state`) to carry the per-atom bonded arrays to the new owner
- [ ] **P2.4** Device tag→local map (sorted tags + binary-search resolve over owned+ghost), rebuilt each reneighbor
- [ ] **P2.5** Re-materialize flat tuple lists from per-atom topology with tag→local resolve + drop tally; wire into `reneighbor`/`partition`
- [ ] **P2.6** RED→GREEN test: tag-map resolves a cross-device tuple; 2-GPU re-materialized tuple set matches single-GPU

## Testing
- `test_cuda_multi_gpu_bonded` (new, GH200): (a) tag→local resolves owned + ghost tags, returns -1 for absent; (b) for a
  small gen_pe system on 2 GPUs, the union of per-GPU re-materialized tuples (mapped back to global tags) equals the
  single-GPU tuple set, each owned by exactly one GPU (no dup/loss); (c) `migrate` then re-materialize preserves topology.
- `test_cuda_multi_gpu` (LJ, 13/13) stays green (per-atom SoA only allocated on the multi-GPU bonded path).

## Out of scope
- Bonded force evaluation + reverse_comm wiring (→ `-03-bond`, `-04-angle-dihedral`).
- special-bonds exclusion (→ `-04`).
- Optimization of the map/re-materialize (→ `-05`).
