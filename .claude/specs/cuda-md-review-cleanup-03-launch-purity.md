---
title: CUDA MD Review Cleanup — link 03 — launch() GPU-purity (zero D2H, zero sync)
slug: cuda-md-review-cleanup-03-launch-purity
status: done
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 2, 8 findings)
chain: cuda-md-review-cleanup-* (link 3/9)
depends_on: cuda-md-review-cleanup-01-component-capability, cuda-md-review-cleanup-02-nlist-skin-gating
---

# CUDA MD Review Cleanup — link 03 — launch() GPU-purity (zero D2H, zero sync)

## Summary

清除 5 个 CUDA MD 组件 hot-path 上的全部 `cudaStreamSynchronize` + 同步 D2H,把 `launch()` 真正变成 "GPU-pure"(`CLAUDE.md ## Architecture § CUDA execution model` 契约 + `component.cuh:122-129`)。具体:
- `cuda::PairLJCut::compute()` 不再 D2H+sync 来定 grid 尺寸(改用 `pair_capacity` 上界 + kernel 早退)
- `cuda::NoseHooverNPT` 一步 5 个 sync + 12 个小 D2H 全部消除,改用单 8-double pinned buffer + `cudaMemcpyAsync` + `cudaEventQuery` 轮询
- `Backend::step()` 每相位双流阻塞 sync 全部删除,只在 `run()` 结束时 sync
- `cuda::AsyncBridge` 的 `defer_stream` sync 改用 `cudaLaunchHostFunc` 触发 host 回调
- `cuda::ZeroMomentum` 把 COM 计算搬到 device(或在 NOTES 中正式认定为 "integrity component" 例外)
- NHC chain 在 host 跑 `std::exp` + 交错 kernel launch → 融合为单 kernel
- `extended_hamiltonian_` / `end_of_step` 的 D2H/H2D 全部转 device-only 或批量到一次 pinned + event

复用 `msu-pe-cuda-components-07` 引入的 `cuda::stream_synchronize` 计数器作为验证仪表,把 `tests/cuda/async/test_hot_path_no_sync.cu` 扩展到驱动完整 NPT 循环(NoseHooverNPT + NeighborList + LJ + HarmonicBond + HarmonicAngle + MultiHarmonicDihedral + VirialStress)并断言每步 sync 计数 == 0。

## Domain basis

`CLAUDE.md ## Architecture § CUDA execution model` 把 CUDA 组件分为 compute(device-only launch)、CPU bridge(D2H + AsyncBridge callback)、host-dependent compute(documented exception)三类。第一类必须 `launch()` 内零 sync 零 D2H;第二类必须走 `AsyncBridge` 的 defer stream 异步路径;第三类必须在 `.claude/NOTES.md` 显式分类。当前实现里:
- LJ pair(应为 compute)走了同步 D2H + sync — 违例
- NPT(应为 compute)5 次 sync + 12 次 D2H — 违例,且无 NOTES 例外条目
- Backend(infra)每相位 sync 双流 — 违例(没有 component 需要,纯保守)
- AsyncBridge `defer_stream` sync 阻塞主 stream — 违例(应用 stream callback)
- ZeroMomentum 同步 D2H — 可争议的 integrity exception,需正式认定

NHC chain 在 host 跑是合理的(算 `std::exp` 串行链),但 **交错** 多个 kernel launch(mtchain=3 → 3 个 kernel)增加了 ~30-50 µs/step,可融合。

性能模型:Hopper GH200 `cudaStreamSynchronize` 在 stream 已空时仍 ~5-10 µs/call;NPT 一步 5 次 + Driver 12 次 = 17 次 = ~100 µs/step 纯 host 等待。每步 1 ms 量级的小系统损失 ~10%。

## Design

完整 8 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 2。本 link 增量设计:

1. **LJ pair grid sizing**(`pair_lj_cut.cu:298-301`):
   - 删除 D2H + sync;改用 `grid_size = (pair_capacity + block - 1) / block`
   - kernel 已有 `if (idx >= *n_pairs) return;`(`pair_lj_cut.cu:63`),正确处理过启动
2. **NPT 8-double 规约 buffer**(`nose_hoover_npt.cu`):
   - 新 device kernel `compute_npt_reductions_k`:同时归约 `{ke_total, ke_x, ke_y, ke_z, stress_xx, stress_yy, stress_zz, energy}` 到一个 8-double 输出
   - pinned host buffer 大小 8 doubles;`cudaMemcpyAsync(defer_stream)` 拷出;`cudaEventQuery` 轮询(host 不阻塞)
   - host 端 NHC `update_eps_dot_` 在轮询到完成后才读取,否则用上一步的值(1-step 滞后对 Nose-Hoover 反馈环安全)
3. **box 对角线 device 版本计数器**(`nose_hoover_npt.cu:505-512`):
   - `f.box` 加一个 `device_version_counter`(单 int);`propagate_x_and_box_` 写盒子后递增
   - host 缓存上次读到的版本;不变则不 D2H
4. **Driver per-phase sync 删除**(`driver.cu:71-72`):
   - 删除 6 phase × 2 stream 的 unconditional sync
   - 只在 `run()` 结束时 sync 一次
   - CPU-bridge 组件用自己的 `cudaEventQuery` 决定 host 可读时机
5. **AsyncBridge stream callback**(`async_bridge.cu:39`):
   - 改用 `cudaLaunchHostFunc` 注册 `on_data_ready` 回调,GPU 完成后异步触发,不阻塞主 stream
6. **ZeroMomentum on-device COM**(`init.cu:233`):
   - 现状:device 算 `{Σmv, Σm}` → D2H → host 算 `v_cm = P/M` → 写回 device
   - 修复:device 算 `{Σmv, Σm}` → device 单线程 kernel 写 `v_cm` 到 resident scratch → `remove_com` 直接读 scratch
7. **NHC chain 融合**(`nose_hoover_npt.cu:300-334, 336-382`):
   - `nhc_atoms_half_` / `nhc_baro_half_`:host 算 `eta_t_dot[mtchain]` 数组 → 单 `scale_vel_chain_k` kernel 接整个数组
8. **extended_hamiltonian on-device**(`nose_hoover_npt.cu:540-545, 598-602`):
   - device 直接算 H,不 D2H;`H` 存在 resident scalar,需要时再读一次

## Files to create or modify

- `src/cuda/potentials/pair_lj_cut.cu` — 删除 `compute()` 内 D2H+sync(行 298-301)
- `src/cuda/md/nose_hoover_npt.cu` — 新 device 规约 kernel;pinned buffer + event 轮询;box 版本计数器;NHC chain 融合;extended_hamiltonian on-device
- `include/atomiverse/cuda/md/nose_hoover_npt.cuh` — `State` 加 pinned buffer、event、版本计数器
- `src/cuda/driver.cu` — `Backend::step()` 删除 per-phase sync(行 71-72)
- `src/cuda/async_bridge.cu` — 改用 `cudaLaunchHostFunc`(行 39)
- `src/cuda/md/init.cu` — `ZeroMomentum::launch()` on-device COM(行 233)
- `include/atomiverse/cuda/md/init.cuh` — `State` 加 device scratch for v_cm
- `tests/cuda/async/test_hot_path_no_sync.cu` — 扩展 NPT 循环
- (新)`tests/cuda/test_pair_lj_no_d2h.cu` — 静态 grep gate
- (新)`tests/cuda/test_driver_no_per_phase_sync.cu` — runtime 测试

## Tasks

- [ ] Remove D2H+sync in `cuda::PairLJCut::compute()`; launch with `grid_size(pair_capacity)`; verify kernel early-exit (`src/cuda/potentials/pair_lj_cut.cu:298-301`)
- [ ] Add `compute_npt_reductions_k` device kernel that writes `{ke_total, ke_x, ke_y, ke_z, stress_xx, stress_yy, stress_zz, energy}` to one 8-double buffer (`src/cuda/md/nose_hoover_npt.cu`)
- [ ] Add pinned host buffer + `cudaEvent_t` + `cudaEventQuery` polling for the 8-double readback; NHC `update_eps_dot_` reads pinned values (1-step lag) (`src/cuda/md/nose_hoover_npt.cu:248,266,399-404`)
- [ ] Add device version counter on `cuda::Box`; `propagate_x_and_box_` increments after write; host caches last-read version and skips D2H if unchanged (`src/cuda/md/nose_hoover_npt.cu:505-512`)
- [ ] Fuse NHC chain: host computes `eta_t_dot[mtchain]` array → single `scale_vel_chain_k` kernel ingests array (`src/cuda/md/nose_hoover_npt.cu:300-334, 336-382`)
- [ ] Move `extended_hamiltonian_` H computation to device-only; remove 2 D2H + 2 H2D (`src/cuda/md/nose_hoover_npt.cu:540-545, 598-602`)
- [ ] Delete per-phase `cudaStreamSynchronize` on both main + defer streams in `Backend::step()`; sync only at `run()` end (`src/cuda/driver.cu:71-72`)
- [ ] Replace `cudaStreamSynchronize(defer_stream)` before `on_data_ready` with `cudaLaunchHostFunc` registering the host callback (`src/cuda/async_bridge.cu:39`)
- [ ] `ZeroMomentum::launch()` on-device COM removal: device computes `{Σmv, Σm}` → single-thread kernel writes `v_cm` to resident scratch → `remove_com` reads scratch directly (no D2H) (`src/cuda/md/init.cu:233`). Alternative: add formal "integrity component" classification row to `.claude/NOTES.md` and keep the host barrier.
- [ ] Extend `tests/cuda/async/test_hot_path_no_sync.cu` to drive a 100-step NPT MD loop with `NoseHooverNPT` + `NeighborList` + `PairLJCut` + `HarmonicBond` + `HarmonicAngle` + `MultiHarmonicDihedral` + `VirialStress`; assert `n_stream_syncs_during_md_step_ == 0`
- [ ] Add `tests/cuda/test_pair_lj_no_d2h.cu` as a static grep gate: assert no `cudaStreamSynchronize` or D2H-direction `cudaMemcpy*` in `src/cuda/potentials/pair_lj_cut.cu` post-fix
- [ ] Add `tests/cuda/test_driver_no_per_phase_sync.cu`: runs 50-step `cuda::md::MD`, asserts `n_stream_syncs_during_md_step_ == 0` after each step

## Testing strategy

1. **Build** — clean
2. **Zero-sync runtime gate** — extended `test_hot_path_no_sync.cu` 跑 NPT 100 步,sync 计数 == 0
3. **Per-component static gate** — `test_pair_lj_no_d2h.cu` 静态 grep
4. **Driver gate** — `test_driver_no_per_phase_sync.cu` 50 步 NVE,sync 计数 == 0
5. **Correctness regression** — `tests/cuda/md/test_cuda_nose_hoover_npt.cu` 全过,T/P 漂移在原误差带内
6. **Perf softgate** — GH200 上跑 PE-NPT(link 09 加),记录 baseline-vs-after,预计 ≥ 1.20× 加速

## Out of scope

- 共享 kernel header 签名分歧 → link 07(cluster 9)
- Warp-shuffle 规约迁移 → link 05(cluster 5)
- `__launch_bounds__` 标注 → link 06(cluster 8)
- 力流池(per-component stream)→ link 06(cluster 8)

## Dependencies

- 必须 **在 link 01 后**:NPT 测试需要 `NoseHooverNPT` 通过 capability(`reads()==Stress`)正确触发 VirialStress 自动注入
- 必须 **在 link 02 后**:不门控邻居列表重建,每步 radix sort 是不可消除的 sync barrier,零-sync 断言无法通过
