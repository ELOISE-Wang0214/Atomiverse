---
title: "arch-refactor-07-molrs — cuda::RDF 迁移 molrs streaming accumulator，删除最后的 C++ 分析数学"
status: code-complete
date: 2026-07-10
chain: arch-refactor (07/09)
depends_on: [arch-refactor-02-halo]
---

# arch-refactor 07/09 — cuda::RDF 迁移 molrs streaming accumulator，删除最后的 C++ 分析数学

## Summary

CPU 侧 molrs 收敛（2026-07-02）之后，全仓库唯一残留的分析数学重复是 legacy g(r)
归一化 `rdf_histogram.{h,cpp}`（~116 行，`rdf_histogram.cpp:5-8` 挂
TODO(roykid, 2026-07-02)），它活着只因 `cuda::RDF`
（`include/atomiverse/cuda/observer/rdf.cuh:141,146`）尚未迁移。本篇把
`cuda::RDF` 从"D2H pair-distances → host 手写直方图/归一化"改为"D2H
positions + per-frame cell → 喂 `cpu::compute::RDFAccumulator`（molrs
streaming accumulator 的 CUDA-safe PIMPL wrapper，`cpu/observer/trajectory_metrics.h`
既有）"，与 CPU `ATV::RDF` 观察者行为对齐；唯一 g(r) 文件发射器
`rdf_write_curve` 从 `rdf_histogram.{h,cpp}` 搬入新的最小 CUDA-safe 文件
`cpu/observer/rdf_curve.{h,cpp}`（保持"cuda 头可 include、不拉 CPU async
观察者栈、无 molrs/`rust/cxx.h`"的既有隔离性质）；随后删除
`rdf_histogram.{h,cpp}` 并加全树 grep 门禁。完成后 "engine measures, molrs
analyzes" 在两个后端上无例外成立，molrs accumulator 成为 g(r) 的单一真源。
附带完成两件治理事项：(a) 在 `.claude/specs/molrs-analysis-ffi.md` 加交叉标注
（本篇关闭其 Task 2 的 CUDA 尾巴；lammps_data 替换前置条件显式记为跨仓库
out-of-scope）；(b) 澄清 build/ 产物治理（`.gitignore:40` 已忽略 `build/`，
本地生成的 corrosion `bridge.h` 过期只是本地陈旧构建树问题，非入库问题——
以 `git ls-files` 复核并记录）。

## Domain basis

径向分布函数（units：r 为 Å，g(r) 无量纲）：

g(r) = ⟨n(r, r+Δr)⟩ / n_id(r)，其中理想气体壳层计数
n_id(r) = ρ · (4π/3)[(r+Δr)³ − r³] · N/2（半对列表 → 归一化时乘 2），
ρ = N/V。legacy `rdf_write`（`src/cpu/observer/rdf_histogram.cpp:43-66`）实现
的正是该式；molrs `RdfAccumulator`（per-frame LinkCell + ideal-gas
shell-volume normalization）实现同一数值定义，且 molrs 侧已单测
streaming == batch **bit-identical**（见 `.claude/notes/notes.md`
[2026-07-02] 条目与 `cpu/observer/trajectory_metrics.h` 头注释）。bin 中心
约定两侧一致：r_i = r_min + (i + 0.5)·Δr，本仓库 r_min = 0。

参考文献：
- Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., OUP 2017,
  DOI: 10.1093/oso/9780198803195.001.0001（§2.6 g(r) 直方图估计与壳层归一化）。
- Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., 2002,
  DOI: 10.1016/B978-0-12-267351-1.X5000-7（§4.4 RDF 采样）。

迁移后 C++ 两个后端均不再包含任何 binning/归一化数学——binning 与归一化的
唯一实现在 molrs（`molrs::rdf_accumulator_new` FFI，经
`cpu::compute::RDFAccumulator` PIMPL 暴露）。精度注记：`cuda::Frame` 位置为
`real_t`（默认 `double`），拓宽为 `double` 后喂 accumulator；默认构建下
CPU/CUDA 在相同帧序列上 g(r) 逐 bin 一致（见 Testing strategy）。

## Design

**`cuda::RDF`（`include/atomiverse/cuda/observer/rdf.cuh`，header-only，
AsyncBridge 系 deferred device-ring-freeze 观察者，形状由 arch-refactor-02
固化）迁移：**

- **捕获载荷改变**：`freeze_slot` 不再 D2D 冻结 pair `dist_sq` 缓冲，改为在
  compute stream 上 D2D 冻结帧位置数组 x|y|z（3 × n_atoms × `real_t`/slot）
  与 per-frame cell（9 × double/slot，支持 NPT——现状仅在 build 时读一次
  box，为已知缺陷）；`d2h_slot` 在 defer stream 上搬到 pinned host slot。
  ring 深度沿用 `AsyncBridge::kRingDepth`（=3）。
- **host 回调**：`on_data_ready` 把 pinned `real_t` 位置拓宽为 double、打包
  blocked x|y|z，调 `acc_.accumulate(positions, cell)`；被拒帧（degenerate
  box / shape mismatch，返回 false）计数，`on_destroy` 时 `ATV_WARN` 汇报
  （对齐 CPU 观察者的 dropped-frames 告警行为）。
- **finalize**：`on_destroy` 调 `acc_.finalize()` 取归一化 g(r)，经
  `cpu::rdf_write_curve` 写文件（r_min=0）；accumulator 为空（0 个被接受帧）
  则不写文件（对齐 CPU 观察者的早退行为）。构造函数直接构造
  `cpu::compute::RDFAccumulator{n_bins, r_max, 0.0}` → 非法参数
  `std::invalid_argument` fail-fast（no-fallback，对齐 CPU）。
- **解除 pair-list 依赖**：删除 `on_build` 里的
  `ctx.runtime.flat_pairs_requested = true`（`rdf.cuh:59`）与
  `can_run` 的 `pair_list` 检查、`neighbor_list.cuh` include。附带性能收益：
  RDF 是 flat-pair 导出的观察者侧唯一请求方，迁移后 NeighborList 的逐步 flat
  导出只剩 Coulomb 请求路径。生命周期不变：`build()` 一次性分配 ring，
  `launch()` 零分配，`destroy()` 释放。
- **所有权**：accumulator 状态（O(n_bins)）由观察者持有，PIMPL 内部才接触
  molrs；`rdf.cuh` 只 include CUDA-safe 头（`trajectory_metrics.h`、
  `rdf_curve.h`），`molrs.hpp` 仍是 CXX FFI 唯一入口、不出现在任何 `.cuh`。

**`rdf_write_curve` 去处（两案取一，本 spec 定案 A）：**

- **A（采纳）**：新建最小文件对 `cpu/observer/rdf_curve.{h,cpp}`，仅含
  `rdf_write_curve`（声明/实现原样搬移，std-only、无 molrs），CPU
  `src/cpu/observer/rdf.cpp:10,52` 与迁移后的 `cuda/observer/rdf.cuh` 共用。
  保持原 `rdf_histogram.h` 头注释声明的隔离性质：cuda 头可 include 而不拉
  CPU async 观察者栈。
- B（否决）：并入 `cpu/observer/rdf.h`/`rdf.cpp`——会迫使 `.cuh` include CPU
  `AsyncObserver`/xtensor 栈，违反 CUDA-safe 隔离，否决。

**删除 + 门禁**：删 `rdf_histogram.{h,cpp}`（`rdf_accumulate`/`rdf_write`
无替代者，直接消失；TODO 同步兑现），`CMakeLists.txt:157` 源表项换成
`rdf_curve.cpp`；全树 grep 门禁见 Testing strategy。`cpu/observer/rdf.h:10-12`
的 legacy-histogram 指引段落删除；`src/cuda/neighbor/neighbor_list.cu`
（:242,324,551,665）、`src/cuda/neighbor/nlist_impl.cuh`（:102）、
`include/atomiverse/cuda/neighbor/neighbor_list.cuh`（:43）、
`include/atomiverse/cuda/component/context.cuh`（:42）中"flat-pair 消费者 =
cuda::RDF"的陈旧注释改为只提 Coulomb（注释级 sweep，无行为改动）。

**跨仓库标注**：`.claude/specs/molrs-analysis-ffi.md` 加两行交叉标注——
(a) 其 Task 2（RDF via molrs）的 CUDA 尾巴由本篇（arch-refactor-07-molrs）
关闭；(b) lammps_data 替换的前置条件（见 Out of scope）。

**build/ 治理澄清**：`.gitignore:40` 已含 `build/`；本地
`build/corrosion_generated/cxxbridge/molrs_cxx/include/molrs_cxx/bridge.h`
缺 `frame_atomic_numbers` / `am1_bcc_assign_frame_from_base`（二者在
molrs-cxxapi/src/bridge.rs:28,39 已声明、被
`src/cpu/semiempirical/am1_bcc_charge_assigner.cpp:33,89` 调用）——即本地
构建树陈旧，重新 configure+build 由 corrosion 再生成即可。以
`git ls-files build/` 复核：若为空（预期），在提交信息里记录澄清即可，
不改 `.gitignore`；若有历史误入库文件（ignore 规则晚于 add），执行
`git rm -r --cached build/` 移出索引（保留本地文件）。

## Files to create or modify

- include/atomiverse/cuda/observer/rdf.cuh
- include/atomiverse/cpu/observer/rdf_curve.h (new)
- src/cpu/observer/rdf_curve.cpp (new)
- include/atomiverse/cpu/observer/rdf_histogram.h（删除）
- src/cpu/observer/rdf_histogram.cpp（删除）
- src/cpu/observer/rdf.cpp
- include/atomiverse/cpu/observer/rdf.h
- CMakeLists.txt（:157 源表项）
- src/cuda/neighbor/neighbor_list.cu（注释 sweep）
- src/cuda/neighbor/nlist_impl.cuh（注释 sweep）
- include/atomiverse/cuda/neighbor/neighbor_list.cuh（注释 sweep）
- include/atomiverse/cuda/component/context.cuh（注释 sweep）
- tests/cuda/observer/test_rdf_molrs.cu (new)
- tests/CMakeLists.txt（注册 test_cuda_observer_rdf_molrs）
- .claude/specs/molrs-analysis-ffi.md（交叉标注）

## Tasks

- [x] Write failing CUDA parity test for the migrated cuda::RDF
      (tests/cuda/observer/test_rdf_molrs.cu (new)) and register it as
      test_cuda_observer_rdf_molrs in tests/CMakeLists.txt: LJ lattice +
      VelocityVerlet{0.0}（静止轨迹，每帧位置相同 → g(r) 与被接受帧数无关），
      跑 cuda::md::MD 数步，解析 g(r) 输出文件，逐 bin 对比 host 侧直接喂
      cpu::compute::RDFAccumulator 同一帧得到的参考曲线（|Δ| ≤ 1e-6，
      文件按 %.6f 格式化），并断言最近邻晶格距离所在 bin 非零；测试同时 grep
      断言 rdf.cuh 不再 include rdf_histogram.h，确保 RED
- [x] Implement include/atomiverse/cpu/observer/rdf_curve.h (new) +
      src/cpu/observer/rdf_curve.cpp (new): move rdf_write_curve verbatim out
      of rdf_histogram.{h,cpp}（std-only、CUDA-safe、Doxygen 带单位），update
      the include in src/cpu/observer/rdf.cpp and swap the source entry in
      CMakeLists.txt:157（rdf_histogram.cpp → rdf_curve.cpp）
- [x] Implement the cuda::RDF migration in
      include/atomiverse/cuda/observer/rdf.cuh: ring 冻结 positions x|y|z +
      per-frame cell（替代 pair dist²），on_data_ready 拓宽/打包后喂
      cpu::compute::RDFAccumulator，on_destroy 走 acc_.finalize() +
      cpu::rdf_write_curve（空 accumulator 不写文件；被拒帧计数 + WARN），
      构造函数构造 accumulator（非法参数 fail-fast），删除
      flat_pairs_requested、pair-list can_run 依赖与 neighbor_list.cuh
      include，Doxygen per doc.style with units
- [x] Remove include/atomiverse/cpu/observer/rdf_histogram.h and
      src/cpu/observer/rdf_histogram.cpp, delete the legacy-histogram
      paragraph in include/atomiverse/cpu/observer/rdf.h, and verify the grep
      gate: `rg -n 'rdf_histogram|rdf_accumulate|\brdf_write\(' include src
      tests examples benchmarks CMakeLists.txt` returns nothing
- [x] Update stale flat-pair comments that name cuda::RDF as consumer
      (src/cuda/neighbor/neighbor_list.cu:242,324,551,665,
      src/cuda/neighbor/nlist_impl.cuh:102,
      include/atomiverse/cuda/neighbor/neighbor_list.cuh:43,
      include/atomiverse/cuda/component/context.cuh:42) to name Coulomb only
- [x] Add cross-annotations to .claude/specs/molrs-analysis-ffi.md: Task 2
      的 CUDA 尾巴由 arch-refactor-07-molrs 关闭；lammps_data.cpp 替换前置
      条件（molrs 修 atom-style 误检测 + 加 read FFI）显式记录，structure
      factor/DCD/TRR/PDB 仍留该 spec
- [x] Verify build/ governance: `git ls-files build/` 为空则在提交信息记录
      澄清（.gitignore:40 已忽略；本地 corrosion bridge.h 过期 = 陈旧构建树，
      重新 configure+build 再生成）；若有 tracked 文件则
      `git rm -r --cached build/`（.gitignore 无需改动）
- [ ] Verify the CUDA backend on the GPU node via sbatch（CLAUDE.md GPU 流程：
      `#!/bin/bash -l` + `module load GPU/buildenv-gcccuda/2026.03-cu13.0`，
      build-cuda + `ctest --test-dir build-cuda --output-on-failure`，含新
      parity test 与既有 MD.RDFOutput）
- [x] Run full check + test suite（CPU login node，ATV_BUILD_CUDA=OFF，
      `cmake --build build -j$(nproc)` + `ctest --test-dir build
      --output-on-failure`）

## Testing strategy

- **Happy path（GPU 节点）**：新 `test_cuda_observer_rdf_molrs`——静止轨迹
  （`VelocityVerlet{0.0}`，满足 DynamicsDriver 的 Force-capability 校验又保证
  逐帧位置恒等，从而 g(r) 对 drop_oldest 掉帧不敏感）下，CUDA 观察者输出
  与 host 侧 `RDFAccumulator` 参考曲线逐 bin 一致（默认 double 构建，
  |Δ| ≤ 1e-6，受 %.6f 文件格式限制）；既有 `MD.RDFOutput`
  （tests/cuda/md/test_md.cu:217）不改断言直接通过（文件头
  `# r [Angstrom]   g(r)` 与 50 行数据格式不变——发射器未变，仍是
  `rdf_write_curve`）。
- **Parity 链（molrs 单一真源）**：CPU 观察者 ≡ accumulator 由既有
  tests/cpu/observer/test_rdf_observer.cpp 保证；CUDA 观察者 ≡ accumulator
  由新测试保证；二者传递即得同轨迹 CPU/CUDA parity，无需跨后端联跑。
- **Edge cases**：0 个被接受帧 → 不写文件；degenerate box / shape mismatch
  → `accumulate` 返回 false、计数并在 destroy 时 WARN（CUDA 侧新对齐行为，
  可在新测试中以非法 cell 分支覆盖或由 molrs 侧既有单测覆盖，二选一在实现
  时定）；非法构造参数（n_bins<1、r_max≤0）→ `std::invalid_argument`。
- **Domain validation**（science.required）：静止晶格的 g(r) 在最近邻晶格
  距离 bin 处非零峰（物理位置正确）；归一化定义与 Allen & Tildesley §2.6
  一致由 molrs 侧单测（streaming==batch bit-identical）承载。
- **删除门禁**：`rg -n 'rdf_histogram|rdf_accumulate|\brdf_write\(' include
  src tests examples benchmarks CMakeLists.txt` 输出为空（`rdf_write_curve`
  以 `\brdf_write\(` 词形排除）。
- **CUDA async 既有护栏**：`test_hot_path_no_sync.cu` /
  `test_async_bridge_tearing.cu` 继续通过（迁移不得引入 step-loop 同步）。
- **全量**：CPU login 节点全量 ctest；CUDA 判据一律标注需 GPU 节点
  （sbatch，见 CLAUDE.md）。

## Out of scope

- **molrs 侧 LAMMPS reader 修复与 `lammps_data.cpp` 替换（跨仓库）**：
  `src/cpu/io/lammps_data.cpp`（~259 行手写 parser）是**有意分歧**——其
  :5-11 注释记录 molrs LAMMPS reader 的 atom-style 自动检测会把 MSU
  `atom_style=molecular` 文件误路由到 charge 分支（mol_id 误标为 type）。
  替换前置条件：molrs 侧修复检测启发式 + bridge 增加 LAMMPS read FFI +
  MSU PE benchmark 回归；已在 .claude/specs/molrs-analysis-ffi.md 交叉标注，
  本篇不动该文件。
- **structure factor / DCD / TRR / PDB readers**：留在
  molrs-analysis-ffi.md（其 Task 4/5），本篇只关闭其 Task 2 的 CUDA 尾巴。
- **CUDA 侧 MSD/VACF 观察者**：CUDA 后端现无此二观察者，不借本篇新增。
- **`long long` histogram → `BigInt` 对齐**：随 rdf_histogram 删除而失效
  （TODO 原文的第二半句），无事可做。
- **NPT + cuda::RDF 的端到端验证**：per-frame cell 捕获使其成为可能，但
  e2e 组合验证不在本篇。
- **08-python / 09-harness 的任何内容**（Python 绑定暴露、harness 接线）。
