---
title: CUDA MD Review Cleanup — link 09 — PE-NPT benchmark + umbrella regression gate
slug: cuda-md-review-cleanup-09-bench-and-regression
status: draft
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 14 + umbrella ac-015, 1 finding + gate)
chain: cuda-md-review-cleanup-* (link 9/9)
depends_on: cuda-md-review-cleanup-01..08 (all)
---

# CUDA MD Review Cleanup — link 09 — PE-NPT benchmark + umbrella regression gate

## Summary

最终关门 link:
1. **添加 PE-realistic benchmark 用例**(`benchmarks/bench_md_step.cu`,cluster 14 唯一 finding):现有 bench 只测 N=100/1000/10000、cutoff=3.0、density=0.8,缺稠密 PE 场景 —— 添加 `bench_md_step_pe_npt`:N=10000-50000、cutoff 10.5 Å + skin 0.3、完整 bonds+angles+dihedrals 拓扑(从 MSU PE deck 抽样)、`cuda::NoseHooverNPT`、T=300 K、P=0 atm
2. **添加 4 个守护测试**(覆盖 link 01-08 的 critical 修复):
   - `tests/cuda/test_launch_purity.cu` — 100 步 NPT MD,断言 `n_stream_syncs_during_md_step_ == 0`
   - `tests/cuda/test_pair_lj_no_d2h.cu` — 静态 grep gate
   - `tests/cuda/test_driver_no_per_phase_sync.cu` — 50 步 MD,sync 计数 == 0
   - `tests/cuda/md/test_neighborlist_skin_gating.cu` — 30 步固定原子,radix sort 调用 ≤ 2(由 link 02 引入,本 link 验证已存在)
3. **Umbrella regression gate** — 重跑 `/mol:review --axis=perf,arch` 在 inventory 同样 scope 上,断言残留 `0 🚨 / ≤3 🔴 / ≤8 🟡`。残留必须在本 link 的 `acceptance.md` 中显式枚举作为已知 backlog(否则 fail)
4. **GH200 perf record** — baseline-vs-after 在 PDC GH200 节点(用户已配 env)上跑 PE-NPT bench;记录每 cluster wall-clock delta;断言:
   - cluster 2(launch-purity)≥ 1.20×
   - cluster 3(skin-gating)≥ 1.30×
   - cluster 5+6+7 合并 ≥ 1.10×
   - 其他 cluster ≥ 0.95×(不退化)

通过本 link 后,inventory `cuda-md-review-cleanup.md` 可改为 `status: done` 或删除。

## Domain basis

`/mol:review` 是 multi-axis static review skill,默认 axis 全开;此处显式 `--axis=perf,arch` 与原 review 范围一致。

`bench_md_step.cu` 现有结构:Google Benchmark 框架,`BENCHMARK(...)->Args({100})`/`Args({1000})`/`Args({10000})`。添加新 case 不应破坏现有 API。

MSU PE deck topology:`benchmarks/md/polyethylene_tension.cpp` 已经构造完整 PE 拓扑;本 link 提取其拓扑 builder 复用到 `bench_md_step_pe_npt`。

GH200 perf 测量约定:在 PDC `aarch64` node 直接运行(per user MEMORY `[[pdc-slurm-accounts]]`,不走 SLURM);CUDA backend 已在 `[[atomiverse-cuda-build]]` 配好。

## Design

完整 finding 见 inventory cluster 14;增量:

1. **PE-NPT bench**:
   - 在 `benchmarks/bench_md_step.cu` 添加 `bench_md_step_pe_npt(benchmark::State& state)`:
     - 复用 `polyethylene_tension.cpp` 的 topology builder
     - `cuda::md::MD` + `cuda::NoseHooverNPT{T=300, tau_T=100, P=0, tau_P=1000, dt=0.5}` + `set_pair(cuda::PairLJCut{...})` + `set_bond/set_angle/set_dihedral`
     - state loop 跑 1000 步预热 + 测 10 块 100 步,记录 mean step time
   - Args:`->Args({10000})->Args({30000})->Args({50000})`
2. **4 个 guard test**:见 Tasks
3. **Umbrella gate**:
   - 提供一个 `scripts/rerun_review_check.sh`:`/mol:review` 在 cuda-md scope 跑一遍,parse 输出,断言总数 + 严重度分布
   - 残留枚举:在 acceptance 中列出预期 ≤8 🟡 已知 backlog,带 file:line + 不修原因
4. **GH200 perf script**:
   - `scripts/bench_md_step_perf.sh`:在 GH200 上跑 baseline(git checkout cuda-md-review-cleanup-pre)和 after,记录 mean step time
   - 脚本输出 JSON:`{cluster_2_speedup: ..., cluster_3_speedup: ..., ...}`
   - acceptance ac-perf-* 断言各项 ≥ 阈值

## Files to create or modify

- `benchmarks/bench_md_step.cu` — 添加 `bench_md_step_pe_npt`
- `tests/cuda/test_launch_purity.cu`(新)
- `tests/cuda/test_pair_lj_no_d2h.cu`(新)
- `tests/cuda/test_driver_no_per_phase_sync.cu`(新)
- `tests/cuda/md/test_neighborlist_skin_gating.cu`(由 link 02 创建,本 link 仅验证存在)
- `tests/CMakeLists.txt` — 注册新 test
- `benchmarks/CMakeLists.txt` — 注册新 bench(已注册 `bench_md_step` 即可)
- `scripts/rerun_review_check.sh`(新,可选)
- `scripts/bench_md_step_perf.sh`(新)
- 最终更新 `cuda-md-review-cleanup.md` inventory:把 status 从 `reference` 改为 `done`

## Tasks

- [ ] Add `bench_md_step_pe_npt(benchmark::State& state)` to `benchmarks/bench_md_step.cu`: N=10000/30000/50000, cutoff 10.5 + skin 0.3, full PE bonds+angles+dihedrals topology, `cuda::NoseHooverNPT` thermostat+barostat
- [ ] Wire PE topology builder reuse from `benchmarks/md/polyethylene_tension.cpp`
- [ ] Add `tests/cuda/test_launch_purity.cu`: 100-step `cuda::md::MD` NPT loop with full force suite (`NoseHooverNPT` + `NeighborList` + `PairLJCut` + `HarmonicBond` + `HarmonicAngle` + `MultiHarmonicDihedral` + `VirialStress`); assert `cuda::detail::sync_count() == 0` after every step
- [ ] Add `tests/cuda/test_pair_lj_no_d2h.cu`: static grep gate — at test runtime, read `src/cuda/potentials/pair_lj_cut.cu` and assert no `cudaStreamSynchronize` / `cudaMemcpy.*HostFromDevice` / `cudaMemcpy.*DeviceToHost` literal occurs (or move to CI lint)
- [ ] Add `tests/cuda/test_driver_no_per_phase_sync.cu`: 50-step `cuda::md::MD` NVE loop with `NeighborList` + `PairLJCut`; assert `cuda::detail::sync_count() == 0` per step
- [ ] Verify `tests/cuda/md/test_neighborlist_skin_gating.cu` (created in link 02) is registered and passes
- [ ] Register all new tests in `tests/CMakeLists.txt`
- [ ] Add `scripts/rerun_review_check.sh`: invoke `/mol:review --axis=perf,arch` on cuda-md scope; parse output to extract 🚨 / 🔴 / 🟡 / 🟢 counts; exit non-zero if 🚨 > 0 or 🔴 > 3 or 🟡 > 8
- [ ] Enumerate expected ≤8 🟡 residual backlog in `acceptance.md` (each with file:line + reason for deferral)
- [ ] Add `scripts/bench_md_step_perf.sh`: run `bench_md_step_pe_npt` on GH200 baseline vs after; emit JSON with per-cluster speedup ratios
- [ ] Run the perf script on PDC GH200; collect speedup numbers
- [ ] Assert speedups: cluster 2 ≥ 1.20×, cluster 3 ≥ 1.30×, clusters 5+6+7 combined ≥ 1.10×, others ≥ 0.95×
- [ ] Update `cuda-md-review-cleanup.md` inventory: change `status: reference` to `status: done`; add a "Completion" section with link summaries + perf delta JSON

## Testing strategy

1. **New bench compiles & runs** — `cmake --build build --target bench_md_step` clean; `./build/benchmarks/bench_md_step --benchmark_filter=pe_npt` 跑通,记录基线
2. **Guard tests pass** — 4 个 guard test 全部 `ctest --test-dir build -R '^test_(launch_purity|pair_lj_no_d2h|driver_no_per_phase_sync|neighborlist_skin_gating)$'` 退 0
3. **Umbrella gate** — `scripts/rerun_review_check.sh` 退 0(0 🚨 / ≤3 🔴 / ≤8 🟡 残留)
4. **GH200 perf gate** — `scripts/bench_md_step_perf.sh` 输出的 JSON 通过 ac-perf-* 阈值断言

## Out of scope

- 任何新功能开发 — 本 link 仅做 benchmark + 守护 + 验证
- 修任何残留 finding — 残留必须显式枚举为 backlog,要么后续 spec 处理,要么留在 inventory 已完成版本中
- 多 GPU bench — CUDA 端单 GPU only(per cluster 13)

## Notes

GH200 perf 测量需要用户在 PDC 上手动运行 `scripts/bench_md_step_perf.sh`(需要 SLURM/直接节点访问);CI 不会自动跑。本 link 的 perf 类 acceptance criteria 标 `type: runtime`,但实际执行是 operator-gated。
