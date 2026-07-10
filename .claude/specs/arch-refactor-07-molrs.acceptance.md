---
spec: arch-refactor-07-molrs
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "cuda::RDF holds no analysis math; feeds cpu::compute::RDFAccumulator"
    type: code
    pass_when: "include/atomiverse/cuda/observer/rdf.cuh includes atomiverse/cpu/observer/trajectory_metrics.h and atomiverse/cpu/observer/rdf_curve.h, holds a cpu::compute::RDFAccumulator member fed from on_data_ready, and contains no histogram vector, no shell-volume/normalization arithmetic, and no reference to rdf_accumulate/rdf_write or rdf_histogram.h."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "cuda::RDF drops pair-list dependency and flat_pairs_requested"
    type: code
    pass_when: "include/atomiverse/cuda/observer/rdf.cuh no longer sets ctx.runtime.flat_pairs_requested, no longer includes atomiverse/cuda/neighbor/neighbor_list.cuh, and its freeze_slot/d2h_slot copy frame positions (x|y|z) plus a per-frame 9-element cell instead of pair dist_sq; can_run no longer checks ctx.runtime.pair_list."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "rdf_write_curve is the single g(r) emitter, in a CUDA-safe file pair"
    type: code
    pass_when: "include/atomiverse/cpu/observer/rdf_curve.h + src/cpu/observer/rdf_curve.cpp exist, declare/define only rdf_write_curve (std-only; no molrs.hpp, no rust/cxx.h), and both src/cpu/observer/rdf.cpp and include/atomiverse/cuda/observer/rdf.cuh call cpu::rdf_write_curve for g(r) file output; CMakeLists.txt compiles rdf_curve.cpp and no longer lists rdf_histogram.cpp."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "rdf_histogram.{h,cpp} deleted with repo-wide grep gate clean"
    type: code
    evaluator_hint: "rg -n 'rdf_histogram|rdf_accumulate|\\brdf_write\\(' include src tests examples benchmarks CMakeLists.txt"
    pass_when: "include/atomiverse/cpu/observer/rdf_histogram.h and src/cpu/observer/rdf_histogram.cpp do not exist, and the hint command returns no matches (rdf_write_curve call sites excluded by the \\brdf_write\\( word shape)."
    status: verified
    last_checked: 2026-07-10
  - id: ac-005
    summary: "Stale flat-pair comments naming cuda::RDF updated to Coulomb-only"
    type: code
    pass_when: "src/cuda/neighbor/neighbor_list.cu, src/cuda/neighbor/nlist_impl.cuh, include/atomiverse/cuda/neighbor/neighbor_list.cuh and include/atomiverse/cuda/component/context.cuh contain no comment naming cuda::RDF as a flat-pair-export or pair-distance consumer, and include/atomiverse/cpu/observer/rdf.h no longer references the legacy histogram pair."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "CUDA/host RDF parity test passes on the GPU node"
    type: runtime
    evaluator_hint: "requires GPU node — sbatch + module load GPU/buildenv-gcccuda/2026.03-cu13.0; ctest --test-dir build-cuda -R test_cuda_observer_rdf_molrs"
    pass_when: "test_cuda_observer_rdf_molrs (tests/cuda/observer/test_rdf_molrs.cu, registered in tests/CMakeLists.txt) passes on the GPU node: static dt=0 LJ-lattice trajectory, parsed cuda::RDF g(r) file matches the host-side cpu::compute::RDFAccumulator reference per bin within 1e-6 (default double build), and the nearest-neighbor lattice-distance bin is non-zero."
    status: pending
  - id: ac-007
    summary: "Existing CUDA RDF e2e and async guards stay green on the GPU node"
    type: runtime
    evaluator_hint: "requires GPU node — ctest --test-dir build-cuda --output-on-failure"
    pass_when: "Full build-cuda ctest passes on the GPU node with MD.RDFOutput (tests/cuda/md/test_md.cu) unmodified in its assertions (header line '# r [Angstrom]   g(r)' + 50 data lines), and test_cuda_hot_path_no_sync / test_cuda_async_bridge_tearing still pass."
    status: pending
  - id: ac-008
    summary: "CPU full check + test suite green on the login node"
    type: runtime
    evaluator_hint: "login node — cmake --build build -j$(nproc) && ctest --test-dir build --output-on-failure (ATV_BUILD_CUDA=OFF)"
    pass_when: "The full CPU build and ctest suite pass with zero failures, including the unchanged tests/cpu/observer/test_rdf_observer.cpp (CPU observer ≡ accumulator leg of the parity chain)."
    status: pending
  - id: ac-009
    summary: "molrs-analysis-ffi.md cross-annotated; build/ tracking clarified"
    type: code
    pass_when: ".claude/specs/molrs-analysis-ffi.md states that Task 2's CUDA tail is closed by arch-refactor-07-molrs and records the lammps_data.cpp replacement preconditions (molrs atom-style detection fix + LAMMPS read FFI) as cross-repo out-of-scope, keeping structure factor/DCD/TRR/PDB in that spec; and `git ls-files build/` returns empty (after `git rm -r --cached build/` if anything was tracked), with .gitignore still containing the build/ rule."
    status: verified
    last_checked: 2026-07-10
out_of_scope:
  - "molrs 侧 LAMMPS reader 修复与 lammps_data.cpp 替换（跨仓库前置）"
  - "structure factor / DCD / TRR / PDB readers（molrs-analysis-ffi spec）"
  - "CUDA 侧 MSD/VACF 观察者新增"
  - "NPT + cuda::RDF e2e 组合验证"
---

# Acceptance — arch-refactor-07-molrs

- **ac-001–ac-003** 固化 "engine measures, molrs analyzes"：CUDA 观察者只做
  freeze/D2H/打包，binning + 归一化唯一实现在 molrs（经
  `cpu::compute::RDFAccumulator` PIMPL），g(r) 文件格式唯一发射器
  `rdf_write_curve` 存活于 CUDA-safe 的 `rdf_curve.{h,cpp}`。
- **ac-004** 是本篇的核心删除门禁：legacy 归一化对 `rdf_accumulate`/`rdf_write`
  连同 `rdf_histogram.{h,cpp}` 消失，全树 grep 为空。
- **ac-006–ac-007 需 GPU 节点**（登录节点无 CUDA toolchain，走 CLAUDE.md 的
  sbatch 流程）；ac-008 在 CPU login 节点执行。静止轨迹设计使 parity 判据对
  异步 drop_oldest 掉帧不敏感（每帧相同 → g(r) 与被接受帧数无关）。
- **ac-009** 收口两件治理事项：跨仓库 out-of-scope 的显式交叉标注，以及
  build/ 产物治理澄清（本地 corrosion `bridge.h` 过期属陈旧构建树，重新
  构建再生成；索引中不得残留 build/ 文件）。

> **审计注记（2026-07-10，本地 macOS）**：ac-001–ac-005/ac-009 本地实证——迁移后 `rdf.cuh`
> 无任何分析数学，删除 grep 门与新 shell 门 `test_cuda_rdf_molrs_migration` 全绿（门内禁词以
> 拼接构造，避免 ac-004 的全树 grep 命中门自身），`git ls-files build/` 为空且 .gitignore 规则
> 保留。ac-008 套件半边本地已跑：134/136，仅 2 个 stash 实证既有 macOS SCF 失败
>（spin_resolution/lindep，Linux 参考常数），零失败判据 owed Linux CI login 节点。
> ac-006/ac-007 owed GH200 GPU 节点（sbatch）。
>
> **迁移中实证的既有隐患（超本篇范围，已留 TODO 注释于 coul_long.cu + neighbor_list.cu）**：
> `cuda::CoulLong` 的 cluster-CSR 输入只在 `flat_pairs_requested` 置位时才被填充，而该 flag
> 自 Verlet 力路径迁移后唯一的 in-tree true-setter 就是旧 cuda::RDF——本篇迁移后无人置位。
> cl_* 指针无条件发布（非空、指向 size-1 占位 buffer），CoulLong 的空指针守卫因此放行并读取
> 从未填充的 CSR；`launch_compute_force_mask` 从未被调用（:375 旧注释已修）。预期 GPU 节点
> 跑 test_ewald_coul / test_pme_coul 会暴露该问题——它先于本篇存在（迁移前也只有"恰好挂了
> RDF observer"的组合才意外可用），修复需给 Coulomb 路径独立的 cluster-CSR 请求（follow-up）。
