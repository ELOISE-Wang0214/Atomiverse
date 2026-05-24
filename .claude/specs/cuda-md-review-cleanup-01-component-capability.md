---
title: CUDA MD Review Cleanup — link 01 — ComponentBase capability-typed auto-injection
slug: cuda-md-review-cleanup-01-component-capability
status: done
created: 2026-05-24
parent: cuda-md-review-cleanup.md (Cluster 1, 5 findings)
chain: cuda-md-review-cleanup-* (link 1/9)
---

# CUDA MD Review Cleanup — link 01 — ComponentBase capability-typed auto-injection

## Summary

在 `cuda::ComponentBase` 上加 `reads()` 和 `needs_neighbor_list()` 两个虚函数,统一 CUDA 端组件的 capability 申报;在 5 个具体类(`VirialStress` / `Angle` / `Dihedral` / `NeighborList` / `NoseHooverNPT`)上正确声明 `role()` / `capabilities()` / `reads()` / `needs_neighbor_list()`;重写 `cuda::md::MD::install_infra_`,用 capability 位匹配替换所有 `name() == "..."` 字符串代理和 `dynamic_cast<const cuda::Pair*>`。`AsyncBridge` 派生类(Thermo / RDF / Mulliken / Recorder / DensityRecorder)各自声明 `role()`。本 link 是 chain 的结构根 —— link 03 的零-sync 测试和后续 link 都依赖正确的 `install_infra_` 注入。

## Domain basis

契约源:`CLAUDE.md ## Architecture § Component tagging`(每个 CUDA MD 组件必须声明 `role()` + `capabilities()`)+ `## Architecture § Categorical setter API`(`set_pair()` 调用时按 capability 自动注入 VirialStress / NeighborList,而非按名字)。CPU 同胞类已经全部按 typed `role()` 注册;CUDA 端必须对齐。

参考 CPU 侧落地:
- `cpu/md/virial_stress.h:58` — `role()==Observer`、`capabilities()==Stress`
- `cpu/observer/rdf.h:47` — `role()==Observer`
- `cpu/recorder/thermo.h:27` — `role()==Recorder`
- `cpu/md/md.cpp` 的 `InfraScan` 用 typed predicate,而非字符串

## Design

完整 5 项发现见 inventory `cuda-md-review-cleanup.md` Cluster 1。本 link 增量设计:

1. `cuda::ComponentBase` 新增:
   ```cpp
   virtual Capability reads() const noexcept { return Capability::None; }
   virtual bool needs_neighbor_list() const noexcept { return false; }
   ```
2. `Capability` enum **不** 新增 `Neighbor` 位 —— 用显式 `needs_neighbor_list()` query 更清晰(neighbor list 不是 frame 字段)。
3. `install_infra_` 改造规则:
   - VirialStress 自动注入条件:存在 component 满足 `(c.reads() & Capability::Stress) != Capability::None`
   - NeighborList 自动注入条件:存在 component 满足 `c.needs_neighbor_list() == true`
   - 删除全部 `comp.name() == "virial_stress"` / `comp.name() == "neighbor_list"` 字符串比较和 `dynamic_cast<const cuda::Pair*>(...)`
4. `cuda::Pair` 基类(在 `potentials/pair.cuh`)重载 `needs_neighbor_list()==true`,所有 pair-style 子类自动继承。

## Files to create or modify

- `include/atomiverse/cuda/component.cuh` — 加两个虚函数声明
- `include/atomiverse/cuda/observer/stress.cuh` — `cuda::VirialStress::role()` + `capabilities()`
- `include/atomiverse/cuda/potentials/angle.cuh` — `capabilities()` 加 `| Stress`
- `include/atomiverse/cuda/potentials/dihedral.cuh` — `capabilities()` 加 `| Stress`
- `include/atomiverse/cuda/potentials/pair.cuh` — `needs_neighbor_list()==true`
- `include/atomiverse/cuda/neighbor_list.cuh` — `role()` 声明(本身是 infra,可声明 `Structural`)
- `include/atomiverse/cuda/md/nose_hoover_npt.cuh` — `reads()==Stress`
- `include/atomiverse/cuda/md/md.cuh` — 重写 `install_infra_`
- `include/atomiverse/cuda/recorder/thermo.cuh` / `recorder.cuh` / `density_recorder.cuh` — 各自 `role()==Recorder`
- `include/atomiverse/cuda/observer/rdf.cuh` / `mulliken.cuh` — 各自 `role()==Observer`
- (新)`tests/cuda/md/test_install_infra_capability.cu` — 验证自动注入按 capability 触发

## Tasks

- [ ] Add `reads()` + `needs_neighbor_list()` virtuals to `cuda::ComponentBase` (`include/atomiverse/cuda/component.cuh`)
- [ ] Add `role()==Observer`, `capabilities()==Stress` to `cuda::VirialStress` (`include/atomiverse/cuda/observer/stress.cuh:22`)
- [ ] Add `Stress` to `cuda::Angle::capabilities()` (`include/atomiverse/cuda/potentials/angle.cuh:19-21`)
- [ ] Add `Stress` to `cuda::Dihedral::capabilities()` (`include/atomiverse/cuda/potentials/dihedral.cuh:19-21`)
- [ ] Override `needs_neighbor_list()==true` on `cuda::Pair` base (or each pair subclass)
- [ ] Add `reads()==Stress` to `cuda::NoseHooverNPT` (`include/atomiverse/cuda/md/nose_hoover_npt.cuh:84`)
- [ ] Declare `role()` on each AsyncBridge-derived recorder/observer: `cuda::Thermo`, `cuda::RDF`, `cuda::Mulliken`, `cuda::Recorder`, `cuda::DensityRecorder`
- [ ] Rewrite `cuda::md::MD::install_infra_` to use `reads()` + `needs_neighbor_list()` + `capabilities()` (drop `name()=="..."` and `dynamic_cast<const cuda::Pair*>`) (`include/atomiverse/cuda/md/md.cuh:243-256`)
- [ ] Add test `tests/cuda/md/test_install_infra_capability.cu`: build a `cuda::md::MD` with `PairLJCut` + `NoseHooverNPT`, assert Plan contains a component with `role()==Observer && (capabilities() & Stress)` and a `needs_neighbor_list()==true` consumer triggered NeighborList injection
- [ ] Verify existing `tests/cuda/md/test_cuda_md_components.cu` still passes
- [ ] `grep -rn` static check: no `name() == "virial_stress"` / `name() == "neighbor_list"` / `dynamic_cast<const cuda::Pair\*>` in `include/atomiverse/cuda/md/`

## Testing strategy

1. **Build** — `cmake --build build --target atomiverse_cuda` 干净通过
2. **Regression** — `ctest --test-dir build -R cuda_md` 全部通过(无行为变化预期)
3. **New typed-injection test** — `test_install_infra_capability.cu` 用 `find_if(plan.components, [](auto& c){ return c.role()==Observer && (c.capabilities() & Capability::Stress) != Capability::None; }) != plan.components.end()` 验证 VirialStress 注入;`needs_neighbor_list()` 同样
4. **Static grep** — 见 ac-006

## Out of scope

- `launch()` 纯度违例修复 → link 03
- 多 GPU halo-exchange 自动注入 → cluster 13(整个 chain out of scope)
- 文件重定位(`cuda/observer/stress` → `cuda/md/virial_stress`)→ link 08
- CPU 端 `InfraScan` 改动 → CPU 已有正确实现
