---
title: CUDA MD Review Cleanup — link 06 — __launch_bounds__ + force-stream pool
slug: cuda-md-review-cleanup-06-launch-config-streams
status: code-complete
created: 2026-05-24
revised: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 8, 5 findings)
chain: cuda-md-review-cleanup-* (link 6/9)
---

# CUDA MD Review Cleanup — link 06 — __launch_bounds__ + force-stream pool

## Summary

调整 launch 配置和并发拓扑:
- Hot force / 规约 kernel 全部加 `__launch_bounds__(256, 4)`(LJ pair、KE、bonded、`kinetic_stress_k`),让编译器按 sm_90 的 64 KiB 寄存器 / 228 KiB shared 上限优化占用率
- `dihedral_multi_harmonic` 寄存器压力(~30+ double 局部变量)逼近 255-reg 上限 → `__launch_bounds__(128, 4)` 强制 spill;若 spill 过大则分解 kernel
- Velocity-Verlet 2 阶 launch 在小 N 时 launch 占主导 — 可选 cooperative-groups grid-sync 一融合 kernel(复杂度高,作为可选)
- `Backend::step()` 把所有 component 串行在单 `context_.stream`,bonded + LJ 独立 compute 却被序列化,SM 空闲 → 引入 **force-stream pool**:Pair / Bond / Angle / Dihedral 各一个 dedicated stream,`cudaEventRecord(force_start)` + 每流 `cudaStreamWaitEvent(force_start)`,end-of-phase 用 events fan-in 到下一相位。预计 2-4× force-kernel 提速
- 替换 LJ pair 的 FP64 `1/r` div + 已算 sqrt(r²)→ `__drcp_rn`/`__frcp_rn` + `export_atom_pairs` 存 r²(不开根)

依赖:link 03 删 Driver per-phase sync(否则多流不起作用),link 04 + 05 稳定 kernel 形状(否则 launch_bounds 测量没意义)。

## Domain basis

`__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)`:nvcc 据此调整寄存器分配。Hopper sm_90 每 SM:
- 64 KiB 32-bit register file
- 228 KiB shared memory(可配)
- 2048 threads max active(64 warps)
- 32 thread blocks max active

`__launch_bounds__(256, 4)` 含义:每 SM 至少 4 个 block,每 block ≤ 256 thread = 1024 thread/SM = 32 warps/SM。寄存器上限 = 64 KiB / 1024 = 64 reg/thread。强制编译器在 64 reg 内完成,超出则 spill 到 local memory。

**Force-stream pool 设计**(NVIDIA 标准 multi-stream concurrency pattern):
```
cudaEvent_t force_start_event;
cudaStream_t lj_stream, bond_stream, angle_stream, dihedral_stream;
// 一次 build 创建所有 streams + events

// per step:
cudaEventRecord(force_start_event, main_stream);
cudaStreamWaitEvent(lj_stream, force_start_event);
launch_lj_kernel<<<..., lj_stream>>>();
cudaEventRecord(lj_done_event, lj_stream);
// same for bond, angle, dihedral
cudaStreamWaitEvent(main_stream, lj_done_event);  // before post_force phase
cudaStreamWaitEvent(main_stream, bond_done_event);
// etc.
```

bonded + LJ 共享 force buffer(atomic 写正确性已由 link 05 cluster-staged 保证),并发无 race。

## Design

完整 5 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 8。本 link 增量:

1. **统一加 `__launch_bounds__`**:
   - `pair_lj_cut.cu`、`kinetic_energy.cu`、`bond_harmonic.cu`、`angle_harmonic.cu`、`multi_harmonic_dihedral.cu`、`observer/stress.cu`(`kinetic_stress_k`)、`md/init.cu`、`md/nose_hoover_npt.cu`、`md/velocity_verlet.cu`、`md/nve_limit.cu`
   - 默认 `__launch_bounds__(256, 4)`;特殊 case 如 dihedral 用 `(128, 4)`
2. **dihedral 寄存器降压**:
   - 测量 `nvcc --ptxas-options=-v -res-usage` 看 reg/thread
   - 若 spill > 128 bytes,分解 kernel:把 `b1*b2`、`b2*b3` 预算到一次,后续复用
3. **Force-stream pool**:
   - `cuda::Context` 加 `std::array<cudaStream_t, 4> force_streams_` + `cudaEvent_t force_start_event_` + 4 个 `force_done_event_`
   - `Backend::dispatch_force_phase()` 调度 Pair / Bond / Angle / Dihedral 各到自己的流
   - end-of-force 用 events fan-in
4. **LJ `1/r` 优化**:
   - `export_atom_pairs_kernel`:存 r²,不开根(`pdist[idx] = r2`,不是 `sqrt(r2)`)
   - LJ kernel:`inv_r2 = __drcp_rn(r2)`(double)或 `__frcp_rn(r2)`(float),`sr2 = sigma² * inv_r2`,后续幂用 `sr2 * sr2 * sr2`

## Files to create or modify

- `src/cuda/potentials/pair_lj_cut.cu` — `__launch_bounds__` + r² (no sqrt) + `__drcp_rn`
- `src/cuda/potentials/bond_harmonic.cu` — `__launch_bounds__`
- `src/cuda/potentials/angle_harmonic.cu` — `__launch_bounds__`
- `src/cuda/potentials/dihedral_multi_harmonic.cu` — `__launch_bounds__(128, 4)` + 寄存器降压
- `src/cuda/md/kinetic_energy.cu` — `__launch_bounds__`
- `src/cuda/observer/stress.cu` — `__launch_bounds__`(若 cluster 5 fuse 没并入此处则继续)
- `src/cuda/md/velocity_verlet.cu` / `nve_limit.cu` / `nose_hoover_npt.cu` / `init.cu` — `__launch_bounds__`
- `src/cuda/md/neighborlist/export_pairs.cu` — 存 r²(不开根)
- `include/atomiverse/cuda/context.cuh` — 加 `force_streams_[4]` + 4 个 events
- `include/atomiverse/cuda/launch.cuh`(若存在,或新)— stream pool 管理
- `src/cuda/driver.cu` — `Backend::dispatch_force_phase()` 多流分发

## Tasks

- [x] Add `__launch_bounds__(256, 4)` to LJ pair kernel (`src/cuda/potentials/pair_lj_cut.cu:42`)
- [x] Add `__launch_bounds__(256, 4)` to KineticEnergy fused kernel + harmonic_bond_k + harmonic_angle_k + nose_hoover_npt reduce_ke_axes_k (`src/cuda/md/kinetic_energy.cu:40`, `src/cuda/potentials/bond_harmonic.cu:21`, `src/cuda/potentials/angle_harmonic.cu:30`, `src/cuda/md/nose_hoover_npt.cu:138`)
- [x] Add `__launch_bounds__(128, 4)` to multi_harmonic_dihedral_k; build clean with no spill warning (`src/cuda/potentials/dihedral_multi_harmonic.cu:36`)
- [ ] **[DEFERRED — user-approved follow-up]** Add force-stream pool to `cuda::Context`: 4 dedicated streams + force_start event + 4 force_done events (`include/atomiverse/cuda/context.cuh`)
- [ ] **[DEFERRED — user-approved follow-up]** Modify `Backend::dispatch_force_phase()` to dispatch each force-evaluator onto its dedicated stream; use cudaEventRecord/cudaStreamWaitEvent for inter-phase synchronization (`src/cuda/driver.cu:54-78`)
- [x] Replace LJ FP64 div chain (sqrt + 1/r + (1/r)²) with single FP64 reciprocal of r² + r² storage in `export_atom_pairs_kernel` (`src/cuda/potentials/pair_lj_cut.cu:97`, `src/cuda/md/neighborlist/export_pairs.cu:73`). *Note*: portable `work_t(1)/r2` shipped; explicit `__drcp_rn` intrinsic (ac-005) is a separate `/mol:fix` follow-up since it has different rounding semantics.
- [ ] (Optional, NOT pursued) Cooperative-groups grid-sync for small-N Velocity-Verlet fused kernel
- [x] Verify correctness: 23/23 cuda tests green; force/energy/stress within tolerances
- [x] Hygiene cleanup (`/mol:simplify`): 5 fixes applied (2 unused includes, 1 dead local, 2 stale comments); ran clean

## Follow-ups (tracked, not blockers)

- **Force-stream pool** (ac-003 + ac-004): pure-perf multi-stream concurrency. Mirrors link 05's cluster-staged deferral. Open a new link or `/mol:spec` issue when scheduling.
- **`__drcp_rn` intrinsic** (ac-005): perf nibble (skips IEEE divide round-trip). Route via `/mol:fix` with explicit perf-vs-numerics tradeoff.
- **Torch-backend signature drift** (link 06 janitor D1+D2): `src/torch/potentials/lj.cu` + `harmonic_bond.cu` callers use the pre-link-06 signatures. Already scoped to **link 07 (kernel-header-and-microfixes)**.
- **`dihedral_multi_harmonic.cu` 115-line indentation drift**: route to `clang-format -i` once a `.clang-format` is wired for `.cu` files.

## Testing strategy

1. **Build + regression** — clean, all cuda_md tests pass
2. **Compile-time check** — `cmake --build build` 通过;无 "uses too much shared / register" 警告
3. **PTX/SASS inspection** — `nvcc --ptxas-options=-v` 输出每 kernel 的 reg/thread <= 64(`dihedral_multi_harmonic` <= 128);spill <= 128 bytes
4. **Correctness** — `tests/cuda/md/test_md.cu` + `test_bonded.cu` + `test_pair_lj.cu` 全过
5. **Perf** — link 09 的 PE-NPT bench 记录 force kernel 单独 wall-clock,预计多流后 2-4× 加速

## Out of scope

- Header / 定义签名分歧(`lj_kernel.cuh`)→ link 07
- Box 按值传 → link 04(已完成)
- Skin gating → link 02(已完成)
