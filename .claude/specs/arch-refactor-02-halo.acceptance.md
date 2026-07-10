---
spec: arch-refactor-02-halo
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "cpu::ComponentBase drops all eight halo virtuals, gains one query hook"
    type: code
    pass_when: "rg -c 'provides_halo_send|provides_halo_recv|comm_forward|comm_reverse|pack_forward|unpack_forward|pack_reverse|unpack_reverse' include/atomiverse/cpu/component/component.h returns no match, and the header declares exactly one halo-related virtual: `virtual HaloParticipant* halo_participant()` defaulting to nullptr. The neighbor trio (needs/provides_neighbor_list, neighbor_cutoff_hint) and the core set (name/role/capabilities/reads/scheduling/frequency/memory_bytes/default_bindings/build/launch/destroy) remain byte-comparable in signature."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "cuda::ComponentBase drops CommField trio, gains halo_participant + memory_bytes"
    type: code
    pass_when: "rg -c 'CommField|comm_forward_fields|comm_reverse_fields' include/atomiverse/cuda/component/component.cuh returns no match, and the header declares `virtual HaloParticipant* halo_participant()` (default nullptr) and `virtual std::size_t memory_bytes() const` (default 0)."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "Unified HaloParticipant protocol headers exist on both backends"
    type: code
    pass_when: "include/atomiverse/cpu/component/halo_participant.h and include/atomiverse/cuda/component/halo_participant.cuh both exist and declare member-for-member identical names: Endpoint{none,send,recv}, halo_endpoint(), CommField{data,width}, comm_forward_fields(), comm_reverse_fields(), bind_participants() — differing only in the field pointer tier (StorageReal* on CPU vs real_t* on CUDA). The .cuh compiles CUDA-free."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "Extra-field forward/reverse round-trip works across ranks via protocol"
    type: runtime
    evaluator_hint: "ctest -R halo_exchange (tests/mpi, >=2 ranks)"
    pass_when: "New MPI test cases in tests/mpi/test_halo_exchange.cpp pass on >=2 ranks: a test-only HaloParticipant declaring width=1 and width=3 fields observes (a) ghost rows of each forward field equal to their owners' values after a step, including across periodic faces, and (b) reverse-field ghost contributions accumulated exactly once onto owner rows."
    status: verified
    last_checked: 2026-07-10
  - id: ac-005
    summary: "User-registered halo endpoints detected via protocol, never double-injected"
    type: runtime
    evaluator_hint: "ctest -R install_infra"
    pass_when: "Test in tests/cpu/test_install_infra_hook.cpp passes: a plan with a manually added HaloSend/HaloRecv pair plus set_pair builds with exactly one send endpoint and one recv endpoint (probed through halo_participant()->halo_endpoint()), and a field-only participant (Endpoint::none) does not suppress auto-injection."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "Infra scan stays capability-typed — no dynamic_cast, no name() matching"
    type: code
    pass_when: "rg -c 'dynamic_cast|name\\(\\) ==' include/atomiverse/cpu/driver/cpu_driver.h returns no match; scan_pending_infra_ derives has_halo_send/has_halo_recv from halo_participant()/halo_endpoint(); InfraScan boolean field names are unchanged and include/atomiverse/cpu/md/md.h + include/atomiverse/cpu/minimizer/minimize.h show no diff versus the pre-spec tree."
    status: verified
    last_checked: 2026-07-10
  - id: ac-007
    summary: "Full build + default ctest gate + MPI parity suite green"
    type: runtime
    evaluator_hint: "cmake --build build && ctest --test-dir build --output-on-failure"
    pass_when: "cmake --build build succeeds and the default ctest gate passes with zero failures, including the tests/mpi multi-rank parity suite (test halo_exchange / lj_energy_bruteforce / bonded_energy_bruteforce / pe multirank targets) under their existing tolerances, unchanged."
    status: pending
  - id: ac-008
    summary: "PE stage-5 MD per-step performance holds after vtable relayout"
    type: performance
    evaluator_hint: "benchmarks/md/common/regression_gate.py"
    pass_when: "The MSU-PE stage-5 benchmark (protocol of .claude/notes/cpu_md_perf_assessment_2026-06.md) rerun on the same host and rank count as a pre-refactor baseline captured from the parent commit, and benchmarks/md/common/regression_gate.py reports PASS comparing the post-refactor stage-5 timing CSV against that baseline."
    status: pending
out_of_scope:
  - "neighbor trio virtuals (needs/provides_neighbor_list, neighbor_cutoff_hint) — arch-refactor-03-neighbor"
  - "Context/BuildContext restructuring — arch-refactor-05-context"
  - "CUDA MultiGPU consumption of declared CommFields (no in-tree CUDA declarant yet)"
  - "driver memory-summary consumption of cuda memory_bytes — arch-refactor-04-engines"
  - "docs / CLAUDE.md architecture narrative — arch-refactor-09-harness"
---

# Acceptance — arch-refactor-02-halo

## ac-001 — CPU 基类删八连、留一钩子

component.h:116-120（provides_halo_send/recv）、:184-186（comm_forward/reverse）、:190-206（pack/unpack 四件）全部删除，替之以单个 `halo_participant()` 查询虚。核心集与 neighbor 三件套签名不动（三件套归 03）。grep 即可判定，无需读 spec。

## ac-002 — CUDA 基类删字段三件、补 memory_bytes

component.cuh:191-201 的 `CommField` 结构体与两个 fields 方法删除；新增与 CPU 同名的 `halo_participant()` 钩子和默认 0 的 `memory_bytes()`（声明对齐，消费归 04）。

## ac-003 — 两后端协议头同形

统一形状是本篇核心交付：成员名逐一相同，仅指针精度档不同（CPU `StorageReal*` / CUDA `real_t*`）。CUDA 头维持 CUDA-free 纪律（同 component.cuh）。

## ac-004 — 声明字段多 rank 往返

对协议的行为性 RED 测试：forward 幽灵刷新（含 PBC 面）与 reverse 恰好一次累加。width=3 用例覆盖连续打包布局 `data[i*width+k]`。本征 fx/fy/fz reverse 通路（eV/Å）不在此判据内——由既有 tests/mpi 回归覆盖（ac-007）。

## ac-005 — 端点探测与防重复注入

手注册端点经 `halo_endpoint()` 被扫描识别，`set_pair` 的 MPI 自动注入不产生第二对 HaloSend/HaloRecv；纯字段参与者不冒充端点。

## ac-006 — 扫描纪律保持

2026-07-02 note 记录的既有纪律（无 dynamic_cast、无 name() 匹配）在协议化后原样成立；`InfraScan` 布尔字段名不变使 md.h/minimize.h 消费点零 diff——该零 diff 本身即判据的一部分。

## ac-007 — 全量门 + 领域奇偶性

clean-break 后全树可构建；多 rank vs 单 rank 能量/力暴力对照维持既有容差，证明 halo 重构未扰动物理。

> **审计注记（2026-07-10，本地 macOS 实施机）**：ac-007 的"零失败"在本机受两个
> stash 实证的既有环境性失败（test_cpu_scf_spin_resolution / test_cpu_scf_lindep，
> macOS/Accelerate vs Linux/OpenBLAS 位一致参考常数，见 01-bindings 台账注记）阻挡，
> 其余 109/111 全绿含全部 tests/mpi 奇偶性套件（n1/n2/n8）——留 pending 待 Linux CI
> 落章。ac-008 需集群基线 + regression_gate.py，owed /mol:bench。

## ac-008 — PE stage-5 性能门禁

链级不变量：组件基类虚表布局改动必须以 stage-5 MSU-PE（10 000 UA 原子）基准回归收尾；基线取本 spec 父提交在同机同 rank 数的运行，`regression_gate.py` PASS 为可观察判定。plan.h 构建期缓存与 backend.cpp:264 零开销路径不变是该判据成立的前提（见 spec Testing strategy）。
