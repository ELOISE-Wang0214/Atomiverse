---
title: Multi-GPU bonded — Angle + Dihedral cross-GPU + special-bonds exclusion
status: approved
created: 2026-06-06
chain: multigpu-bonded
link: 4
---

# Multi-GPU bonded — Angle + Dihedral cross-GPU + special-bonds exclusion

## Summary
链条第 4 环。把 -03 的模式推广到 `AngleHarmonic`(center=owner)与 `DihedralMultiHarmonic`(2nd atom=owner):
per-GPU 重展开 angle/dihedral tuple(成员解析到 owned/ghost)→ 核 atomicAdd 所有成员(含 ghost)→ `reverse_comm`。
并在多卡 `PairLJCut` 接入 special-bonds 跨域排除:排除 partner 全局 tag 经 tag→local 解析(含 ghost),
`setup_force_engine` 调用 `set_special_bonds`。验证完整 PE(LJ+键+排除)多卡 == 单卡 == CPU。

## Domain basis
angle/dihedral 局域多体力,各由唯一 owner 卡持有。special-bonds 排除必须解析跨域 partner(ghost),否则
近邻成键对会被 LJ 误算(单卡 MSU-PE 崩溃的同源问题)。完整 PE 能量分解(E_pair/E_bond/E_angle/E_dihedral)守恒一致。

## Design
- angle/dihedral force:复用单卡 `harmonic_angle_k`/`multi_harmonic_dihedral_k`,输入为 -02 重展开的 per-GPU tuple;
  reverse_comm 已覆盖 ghost 力回传(三/四成员同理)。
- special-bonds:多卡 PairLJCut 的排除表用 per-atom `nspecial/special`(tag)解析到 local(含 ghost)构建;每 reneighbor 重建。
- 完整 PE 多卡 vs 单卡 vs CPU 对照测试。

## Files to create or modify
- `src/cuda/md/multi_gpu.cu`、`.cuh`(angle/dihedral force + special-bonds setup)。
- `src/cuda/potentials/pair_lj_cut.cu`(若排除构建需 owned+ghost 范围参数化)。
- `tests/cuda/md/test_cuda_multi_gpu_bonded.cu`(完整 PE 等价)。

## Tasks
- [ ] **P4.1** Angle cross-GPU (re-materialized tuples + reverse_comm)
- [ ] **P4.2** Dihedral cross-GPU
- [ ] **P4.3** special-bonds cross-domain exclusion in multi-GPU PairLJCut (resolve excluded tags incl. ghosts)
- [ ] **P4.4** RED→GREEN: full PE (LJ+bonded+exclusions) 2/4-GPU == single == CPU

## Testing
- Full-PE energy decomposition + forces parity 2/4-GPU vs single vs CPU; special-bonds-excluded pairs verified; ΣF≈0.

## Out of scope
- Optimization (→ -05); the benchmark campaign (→ benchmark-matrix).
