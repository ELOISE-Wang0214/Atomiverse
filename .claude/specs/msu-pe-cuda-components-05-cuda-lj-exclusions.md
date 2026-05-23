---
title: CUDA PairLJCut special-bonds 1-2/1-3/1-4 exclusions (chain link 5/8)
status: approved
created: 2026-05-23
parent: msu-pe-validation.md (Phase G.8)
chain: msu-pe-cuda-components-* (link 5/8)
---

# CUDA PairLJCut special-bonds 1-2/1-3/1-4 exclusions (chain link 5/8)

## 综述 (Summary)

为 `cuda::PairLJCut` 增加 LAMMPS 风格的 `special_bonds` 排除支持:暴露 `set_special_bonds(double s12, double s13, double s14)` 接口,在 `build()` 阶段从 `cuda::Frame` 的 `bonds`/`angles`/`dihedrals` 拓扑里枚举每个原子的 1-2/1-3/1-4 排除伙伴并上传到 device,在 LJ 核内对每对候选邻居用二分查找排除伙伴列表(与 CPU 版语义和算法一致)。空排除列表时核走零开销快路径,等价于现有未掩码的核。该补丁让 CUDA 端能够正确模拟键合聚合物 —— 这是父规范 `msu-pe-validation` Phase G.8 的硬正确性缺陷,未修复前任何 CUDA 上跑的 MSU PE 链都会被静默腐败。

## Domain basis

LAMMPS `special_bonds` 约定(参见 `lammps-src/src/pair_lj_cut.cpp` 与 `lammps-src/src/special.cpp`):对每个键合原子对 (i,j),若 j 是 i 的 1-2 / 1-3 / 1-4 拓扑伙伴,pair 能量与力按 `special_lj[1..3]` 缩放;`0.0` 表示完全排除该对,`1.0` 表示不做任何修改。MSU 非晶 PE 协议使用 Dreiding 联合原子模型 (UA),约定 `special_bonds lj/coul 0.0 0.0 0.0`,即 1-2/1-3/1-4 LJ 贡献被完全排除(DOI: 10.1021/j100389a010,Dreiding force field)。CPU 端 `cpu::PairLJCut::set_special_bonds` 已经严格实现该约定(`include/atomiverse/cpu/potentials/pair_lj_cut.h`、`src/cpu/potentials/pair_lj_cut.cpp` L44–L122),本规范在 CUDA 端做语义等价移植。

单位:epsilon [eV],sigma [Å],cutoff [Å],缩放因子 s12/s13/s14 无量纲(典型取值 0.0 全排除或 1.0 不排除)。

## Design

新增字段(在 `cuda::PairLJCut` 上):

- `double s12_ = 1.0, s13_ = 1.0, s14_ = 1.0`(host-side 配置)。
- `int n_excl_total_ = 0`(host-cached,排除对总数;==0 ⇒ 快路径)。
- 设备端 CSR 排除表:
  - `int* d_excl_offsets_ = nullptr`(device,长度 `n_atoms + 1`,前缀偏移)。
  - `int* d_excl_partners_ = nullptr`(device,长度 `n_excl_total_`,每段升序排序的 partner 原子 *local index*)。

排除列表用 **device-local atom index** 而非 global tag —— 与 CPU 用 tag 不同,因为 CUDA 路径目前是单卡(链 5/8 是单 GPU 正确性补丁,多 GPU/MPI 不在范围内),且 `cuda::Frame` 拓扑 `bonds.atomi/atomj`、`angles.atomi/atomk`、`dihedrals.atomi/atoml` 本身就是 local index 形式。这避免了为每个 candidate pair 在 kernel 内做 tag 间接,命中显存带宽。

生命周期:

- 构造:`PairLJCut(epsilon, sigma, cutoff)` 不变;新增 `set_special_bonds(s12, s13, s14)` 链式返回 `PairLJCut&`(CUDA 类不需要右值重载,CPU 那对右值版只为 fluent builder;CUDA 端组件按值移入 `on(...).add(...)`,单一左值 setter 就够)。
- `build(ctx)`:基类目前没有 `build` override;新增 `build()` override 在该阶段调用 `rebuild_exclusions_device(ctx.frame)` —— host 端按 CPU 同一算法枚举对、排序去重,然后 `cudaMalloc` + `cudaMemcpy` 把 CSR 上传到 device。任意 `s12/s13/s14 < 1.0`(按 LAMMPS 约定取 `== 0.0` 才进入完全排除,但允许部分缩放)触发构表;全 1.0 ⇒ `n_excl_total_ = 0`,跳过任何 alloc。
- `compute(ctx)`:把 `d_excl_offsets_` / `d_excl_partners_` / `n_excl_total_` 透传给 kernel;`n_excl_total_ == 0` 时传 `nullptr` 并由 kernel 走快路径。
- `destroy(ctx)`:新增 override,`cudaFree` 设备端 CSR。

Kernel 改造(`lj_force` in `lj_kernel.cuh` / `pair_lj_cut.cu`):

- 新增参数 `const int* __restrict__ excl_offsets, const int* __restrict__ excl_partners, work_t s12, work_t s13, work_t s14`(保持向后兼容:当 `excl_offsets == nullptr` 时整段排除/缩放逻辑被编译期/运行期短路)。
- 当前 LAMMPS-style 1-2/1-3/1-4 三档共享同一 device 表 —— 这里的简化是:MSU PE 协议要求三档全 0,缩放因子整体压缩为单个 `s_excl` 字段(对每对要么完全跳过、要么按 `s_excl` 倍缩放)。若调用者用混合缩放(如 (0, 0, 0.5)),按 CPU 实现一样用 *合并* 排除集 + 单一缩放是不够的。**作用域决定**:本链遵循 CPU 端语义(CPU 也把 `s12_/s13_/s14_` 折叠成一张"任意非 1 → 加入排除集"),CUDA 端按相同折叠规则实现,但只支持 *三个 s 相等* 或 *三个 s 同为 0* 两种合法组合(覆盖 MSU PE 用例);其它组合在 `set_special_bonds` 中 `throw std::invalid_argument`。这一限制在 doxygen 中写明,链 7 (CUDA validation) 不会触发。
- 排除查询:每个 thread 处理一个 pair (i,j):
  - 读 `int beg = excl_offsets[i], end = excl_offsets[i+1];`
  - `if (end > beg) { binary_search(excl_partners + beg, excl_partners + end, j); if (found) { if (s == 0) return; else energy *= s; f_over_r *= s; } }`
  - 仅查 i → partners 一次足够(对称构建保证 j 的 partner 列表也包含 i)。

新增的 `Pair` 基类 capability bit:本链不动 `pair.cuh` 的 `Capability` 枚举 —— `set_special_bonds` 是 PairLJCut 私有 API,capability 不变。

## Files to create or modify

- `include/atomiverse/cuda/potentials/pair_lj_cut.cuh` — 增加 `set_special_bonds` setter、`build/destroy` override 声明、device CSR 指针字段、私有 `rebuild_exclusions_device` 辅助。
- `include/atomiverse/cuda/potentials/lj_kernel.cuh` — `lj_force` 签名新增 `excl_offsets/excl_partners/s_excl` 参数;保持快路径分支。
- `src/cuda/potentials/pair_lj_cut.cu` — 实现 setter / build / destroy / 改写 launch;实现新版 `lj_force` kernel。
- `tests/cuda/potentials/test_pair_lj_special_bonds.cu` (new) — 4 原子链、3 键的单元测试套件。
- `tests/cuda/potentials/CMakeLists.txt` — 注册新测试文件。

## Tasks

- [ ] Write failing tests for cuda::PairLJCut special-bonds API and 4-atom-chain exclusion correctness (tests/cuda/potentials/test_pair_lj_special_bonds.cu)
- [ ] Extend lj_force kernel signature with excl_offsets/excl_partners/s_excl parameters and per-pair binary-search masking (include/atomiverse/cuda/potentials/lj_kernel.cuh, src/cuda/potentials/pair_lj_cut.cu)
- [ ] Implement set_special_bonds setter and validate allowed (s12,s13,s14) combinations on cuda::PairLJCut (include/atomiverse/cuda/potentials/pair_lj_cut.cuh, src/cuda/potentials/pair_lj_cut.cu)
- [ ] Implement build() override that enumerates 1-2/1-3/1-4 partners from frame topology, sorts/dedups, and uploads CSR to device (src/cuda/potentials/pair_lj_cut.cu)
- [ ] Implement destroy() override that frees device CSR and zero-init pointers (src/cuda/potentials/pair_lj_cut.cu)
- [ ] Wire compute() to pass exclusion CSR through to lj_force and take nullptr fast path when n_excl_total_ == 0 (src/cuda/potentials/pair_lj_cut.cu)
- [ ] Add doxygen for set_special_bonds, build, destroy with unit annotations and supported s-combination contract (include/atomiverse/cuda/potentials/pair_lj_cut.cuh)
- [ ] Verify CUDA energy against CPU PairLJCut on the 4-atom chain (full-exclusion and no-exclusion configs) within 1e-5 eV (tests/cuda/potentials/test_pair_lj_special_bonds.cu)
- [ ] Register test in tests/cuda/potentials/CMakeLists.txt and run full check + test suite

## Testing strategy

- **Happy path — 4-atom chain (0,0,0) full exclusion:** 4 共线原子在 σ 附近,3 键(0-1, 1-2, 2-3)、2 角(0-1-2, 1-2-3)、1 二面角 (0-1-2-3)。CUDA `set_special_bonds(0.0, 0.0, 0.0)` 后 LJ 能量应严格为 0(所有 pair 都被排除)。
- **Happy path — 4-atom chain (1,1,1) no exclusion:** 同样几何,`set_special_bonds(1.0, 1.0, 1.0)` (或不调用 setter),CUDA LJ 能量与不带特殊键的等价 CPU `PairLJCut` 计算在 `1e-5 eV` 内一致。
- **Fast-path regression:** 一个 256 原子无键合 LJ 流体盒,不调用 `set_special_bonds`,CUDA 能量/力与改造前 baseline (record-replay 或同位 CPU 比对) bit-equivalent;确认 `n_excl_total_ == 0` 时 kernel 不读 CSR、`d_excl_offsets_` 保持 `nullptr`。
- **Mixed-scaling rejection:** `set_special_bonds(0.0, 1.0, 0.0)` 或类似非法组合调用应抛 `std::invalid_argument`。
- **Domain validation:** CUDA vs CPU 跨后端能量比对 — 同一 4-原子链同一 (s12,s13,s14),CPU `cpu::PairLJCut::set_special_bonds` 与 CUDA `cuda::PairLJCut::set_special_bonds` 给出的总能量差 < `1e-5 eV` (相对 < `1e-4`)。这是 Phase G.8 的退出门。
- **Lifecycle:** `build → compute → destroy → build → compute` 重复一次,验证 device CSR 没有泄漏(`cudaFree` 后 `nullptr`)且第二次 build 重新分配。
- Full `ctest -L cuda` 在 GH200 节点上通过。

## Out of scope

- 三个 s 不全相等也不全为 0 的混合缩放组合 —— 本链拒绝;后续若 LAMMPS 输入需要可作为独立链 `cuda-pair-lj-mixed-special-bonds` 扩展(须改 kernel 为查表式三档缩放)。
- 多 GPU / CUDA+MPI 下的跨 rank tag-based 排除 —— 整个 MSU PE 父规范明确 CUDA 是单卡范围(fixed scope decision #1)。
- 排除支持扩展到其他 CUDA pair styles(`pair_lj_long`, `pair_coul` 等)—— 那些 pair style 在 CUDA 端尚不存在;本链只动 `pair_lj_cut`。
- Kernel SoA 优化(如 warp-cooperative binary search)—— 本链优先正确性;性能优化推迟到 Phase F 性能链。
