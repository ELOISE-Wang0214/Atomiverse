---
spec: arch-refactor-03-neighbor
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "Old neighbor trio scrubbed from the tree"
    type: code
    pass_when: "rg -n 'needs_neighbor_list|provides_neighbor_list|neighbor_cutoff_hint' include/ src/ tests/ .claude/architecture.md returns zero matches."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "CPU opt-in protocols exist with the specified shape"
    type: code
    pass_when: "include/atomiverse/cpu/component/neighbor_protocol.h defines NeighborConsumer (pure virtual Real neighbor_cutoff() const, doxygen documents Å and the 0-means-no-cutoff semantic) and NeighborProvider; cpu::ComponentBase declares as_neighbor_consumer()/as_neighbor_provider() returning nullptr by default and no other neighbor-related virtual."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "CUDA mirror is same-shaped and closes the cutoff gap"
    type: code
    pass_when: "include/atomiverse/cuda/component/neighbor_protocol.cuh mirrors both protocols (double neighbor_cutoff), cuda::Pair overrides as_neighbor_consumer() to expose its ctor cutoff, cuda::NeighborList overrides as_neighbor_provider(); rg dynamic_cast over include/atomiverse/cpu/driver/cpu_driver.h and include/atomiverse/cuda/md/md.cuh returns zero matches."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "Protocol unit tests pass on both backends"
    type: runtime
    pass_when: "test_install_infra_hook, test_cpu_md_set_api and the CUDA test built from tests/cuda/md/test_install_infra_capability.cu pass, asserting: default accessors nullptr; PairLJCut consumer with neighbor_cutoff()==ctor cutoff; NeighborList provider; PairMorse and NoseHooverNPT return nullptr consumers; CoulLong consumer reports cutoff_."
    status: pending
  - id: ac-005
    summary: "Missing-NL throw and injection behavior unchanged after migration"
    type: runtime
    pass_when: "With a pair evaluator set and no neighbor list, MD/Minimize/Evaluate/MC build() throws naming the offending component and suggesting set_neighbor; with set_neighbor the same plan builds — asserted by the existing set-API test binaries passing without weakened assertions."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "MPI halo auto-injection parity preserved"
    type: runtime
    evaluator_hint: "requires ATV_HAS_MPI build (buildenv-gcc/2026.03-mpich)"
    pass_when: "tests/mpi/test_halo_exchange.cpp, tests/mpi/test_lj_energy_bruteforce.cpp and tests/mpi/test_coul_energy_bruteforce.cpp pass under the MPI build; the LJ set_pair path still injects HaloSend{1.1*pair_cutoff}+HaloRecv exactly once (explicit-halo dedup intact)."
    status: verified
    last_checked: 2026-07-10
  - id: ac-007
    summary: "PE stage-5 hot-path performance gate"
    type: performance
    evaluator_hint: "marker: pe-stage5 (benchmarks/md/pe/performance/polyethylene_tension.cpp)"
    pass_when: "On the same CPU node/config, median wall time per step of polyethylene_tension (identical --scale, >=3 runs) is within 3% of the pre-change dev baseline and the timing report's Neigh category share differs by <=2 percentage points; git diff for src/cpu/neighbor/cell_list.cpp is empty."
    status: pending
  - id: ac-008
    summary: "Full check + test suite green on both backends"
    type: runtime
    pass_when: "CPU: cmake build + full ctest pass on the login node (Release, ATV_BUILD_TESTS=ON, ATV_BUILD_RHF=ON). CUDA: build-cuda + ctest pass in a GH200 batch job with ATV_BUILD_CUDA=ON."
    status: pending
out_of_scope:
  - "cell_list.cpp 重建算法优化（cell-list-omp-finalize spec）"
  - "Context::pair_list 归属（arch-refactor-05-context）"
  - "first-consumer→max-over-consumers cutoff 语义升级"
  - "torch 后端"
---

# Acceptance — arch-refactor-03-neighbor

ac-001–003 锁死 clean-break 的接口形状（旧三虚函数零残留、双协议 + 双访问器、CUDA 同形状且补齐 cutoff、探测零 dynamic_cast）；ac-004–006 锁死行为回归面（协议断言、四 driver 缺邻表抛错、set_pair 触发的 NeighborList/Halo 注入行为不变，含显式 halo 去重与 CoulLong 的既定增量）；ac-007 是链级热路径门禁（PE stage-5 计时 + cell_list.cpp 零改动）；ac-008 是全量收尾。

## AC-001 — Old neighbor trio scrubbed from the tree

clean-break 验证：`needs_neighbor_list` / `provides_neighbor_list` / `neighbor_cutoff_hint` 三个旧符号在 include/、src/、tests/ 与 .claude/architecture.md 中零命中（含 doc 注释）。

## AC-002 — CPU opt-in protocols exist

neighbor_protocol.h 的协议形状与 Design 一致；ComponentBase 上只剩两个可空访问器，无其他 neighbor 相关虚函数。

## AC-003 — CUDA mirror same-shaped

CUDA 协议镜像 + cuda::Pair 首次上报 cutoff；探测路径零 dynamic_cast。

## AC-004 — Protocol unit tests

两后端协议单测按 Design 迁移映射逐组件断言。CUDA 判据需 GH200 批作业验证。

## AC-005 — Missing-NL throw unchanged

四类 driver 缺邻表抛错行为与文案契约（点名 offender + set_neighbor 指引）不减弱。

## AC-006 — MPI halo parity

MPI 构建下三个 bruteforce/halo 测试全绿；LJ 路径 halo 注入次数与尺寸公式逐行不变；CoulLong-only 场景的注入增量按 Design 以断言固化。

## AC-007 — PE stage-5 performance gate

同节点同配置 ≥3 次取中位数，step 墙钟与迁移前基线差 ≤3%，Neigh 份额差 ≤2 个百分点；cell_list.cpp diff 为空。

## AC-008 — Full suite green

CPU 登录节点全量 ctest + GH200 批作业 CUDA 构建与 ctest。

> **审计注记（2026-07-10，本地 macOS）**：ac-004 的 CPU 腿（test_install_infra_hook、
> test_cpu_md_set_api）全绿，CUDA 腿（test_install_infra_capability.cu）需 GH200 编译
> 运行——整条判据留 pending owed GPU 节点。ac-007 owed /mol:bench（集群基线）。ac-008
> CPU 侧 109/111（2 个 stash 实证既有 macOS 失败，见 01 台账），CUDA 侧需 GPU 节点——
> 留 pending owed Linux CI + GH200。原 spec 测试文件中的 has-member 删除门 static_assert
> 因其字面 token 会触发 ac-001 grep 门禁而移除；删除性由 ac-001 grep + 全部迁移调用点
> 可编译共同保证。
