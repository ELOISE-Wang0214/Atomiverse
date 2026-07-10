---
title: arch-refactor-02-halo — 从 ComponentBase 拆出 HaloParticipant opt-in 协议
status: code-complete
date: 2026-07-10
chain: arch-refactor (02/09)
depends_on: [arch-refactor-01-bindings]
---

# arch-refactor-02-halo — 从 ComponentBase 拆出 HaloParticipant opt-in 协议

## Summary

把 halo 通信契约从两个后端的 `ComponentBase` 虚表中整体拆出，收进一个新的 opt-in 协议 `HaloParticipant`。CPU 侧删除八个 halo 虚函数（`provides_halo_send/recv`、`comm_forward/reverse`、`pack_forward/unpack_forward/pack_reverse/unpack_reverse`，component.h:116-120、184-186、190-206），CUDA 侧删除结构体式 halo 字段方法（`CommField` + `comm_forward_fields/comm_reverse_fields`，component.cuh:191-201），两侧各以**一个**查询虚函数 `halo_participant()`（默认 `nullptr`）替代——CPU 虚表从 22 个虚函数瘦到 15 个，绝大多数组件（NoseHoover/Virial/TranslateMove 等只 override 4-5 个虚的组件）不再继承与己无关的通信 stub。协议采用 CommField 式结构化声明，统一目前已分叉的两后端 halo 形状；CPU 的 int 索引 pack/unpack 语义收进协议消费端（halo 引擎按声明打包），`HaloSend/HaloRecv` 与 driver 基建扫描全部迁移到协议探测。clean-break，无兼容层；热路径虚调度结构与 `Driver` 命名保持不变。

## Domain basis

**LAMMPS comm_forward/comm_reverse 模式（设计出处）。** 被删的八连注释自述源自 LAMMPS `src/pair.h:50-52,207-210` 的 forward（owner→ghost，如 EAM 嵌入导数）/ reverse（ghost→owner 累加，如 EAM 密度部分和）每步幽灵通信；参照 Thompson et al., *Comput. Phys. Commun.* **271**, 108171 (2022), DOI: 10.1016/j.cpc.2021.108171。物理不变量：forward 传输的是每原子状态量（单位随字段，典型为 Å 或无量纲），reverse 传输的是跨子域被拆分的部分和，必须**恰好累加一次**回 owner——本篇重构不得改变既有 3-double(fx/fy/fz, eV/Å) 本征 reverse 通路及其 Newton 第三定律语义（halo_exchange.h 的 2026-05-06 reverse-payload audit 仍然成立）。多 rank 与单 rank 的能量/力奇偶性（tests/mpi 暴力对照测试）是本篇的领域判据。

**god-interface 指控范围澄清（经逐文件核实）。** 对 `ComponentBase` 的"上帝接口"指控只有 halo 部分成立：MC、optimizer、profiling 均**不在**虚表上（它们经 Context 可空字段/组件角色走，profiling 在虚表上仅 `memory_bytes()` 一个）。且两条 halo API 均为死声明——CPU 八连全树零 override、零调用方（`src/cpu/comm/halo_exchange.cpp` 用内部 helper `pack_positions_for_face` 等打包本征字段，从不虚分派）；CUDA `comm_forward_fields/comm_reverse_fields` 同样零消费方。两后端 halo API 已分叉（CPU 命令式 int pack/unpack vs CUDA 声明式 CommField），这是拆独立协议的最强证据；因此本篇**只拆 halo**，不发明 CapabilityProvider/ProfilerParticipant/FramePublisher（见 Design 已否决备选）。

## Design

### 触及实体

- `ATV::cpu::ComponentBase`（include/atomiverse/cpu/component/component.h:50-207）：22 虚 → 15 虚。核心集（name/build/launch/destroy/role/capabilities/reads/scheduling/frequency/memory_bytes/default_bindings）与 neighbor 三件套（needs/provides_neighbor_list、neighbor_cutoff_hint，留给 03）原样保留。
- `ATV::cuda::ComponentBase`（include/atomiverse/cuda/component/component.cuh:78-202）：删 `CommField` 嵌套结构体与两个 fields 方法，补 `memory_bytes()`（默认 0，对齐 CPU 核心集；驱动端内存汇总消费归 04-engines），加同名 `halo_participant()` 钩子。虚函数数 15 → 15（减 2 加 2）。
- `cpu::comm::HaloSend/HaloRecv`（include/atomiverse/cpu/comm/halo_exchange.h）：唯二的 halo 端点组件，迁移到协议（多继承）。`NeighborList`/`Migrate` 经查**不**实现任何 halo 虚函数——它们通过 `DomainDecomposition` 的 generation 计数器与 halo 耦合，本篇不触碰（neighbor 归 03）。
- `CPUDriver::InfraScan / scan_pending_infra_()`（include/atomiverse/cpu/driver/cpu_driver.h:334-376）：halo 探测改协议探测；`InfraScan.has_halo_send/has_halo_recv` 布尔字段名不变，故 `md.h:296-311` 与 `minimize.h:118-139` 两处自动注入消费点**零改动**。
- `src/cpu/comm/halo_exchange.cpp`（打包 helper 群，约 :277-345）：新增引擎侧声明字段打包，并入既有每面消息。

### 新公共符号（最终形状）

```cpp
// include/atomiverse/cpu/component/halo_participant.h (new)
namespace ATV { namespace cpu {
class HaloParticipant {
  public:
    virtual ~HaloParticipant() = default;

    /// 本征交换端点：HaloSend → send，HaloRecv → recv，纯字段参与者 → none。
    enum class Endpoint { none, send, recv };
    virtual Endpoint halo_endpoint() const { return Endpoint::none; }

    /// 每原子额外通信字段（声明式，两后端形状统一）。
    struct CommField {
        StorageReal* data = nullptr; ///< 每原子数组（owned+ghost 容量），CUDA 镜像为 real_t*
        int width = 1;               ///< 每原子分量数，连续打包 data[i*width + k]
    };
    /// 每步 forward（owner→ghost 刷新）的额外字段。默认无。
    virtual std::vector<CommField> comm_forward_fields() const { return {}; }
    /// 每步 reverse（ghost→owner 累加）的额外字段。默认无。
    virtual std::vector<CommField> comm_reverse_fields() const { return {}; }

    /// 端点专用：build 期由 driver 注入字段参与者列表（默认 no-op）。
    virtual void bind_participants(const std::vector<const HaloParticipant*>&) {}
};
}} // namespace
```

`include/atomiverse/cuda/component/halo_participant.cuh`（new）为同形镜像（`real_t* data`，成员名逐一相同），CUDA-free 头（同 component.cuh 纪律）。

`ComponentBase` 两侧各加**一个**钩子：`virtual HaloParticipant* halo_participant() { return nullptr; }`；参与组件 override 返回 `this`（多继承：`class HaloSend : public cpu::ComponentBase, public cpu::HaloParticipant`）。`scan_pending_infra_()` 相应放弃 const 限定（driver 内部细节，md.h/minimize.h 调用点均为非 const 方法，不受影响）。

### 定形理由（声明式 CommField，而非命令式 pack/unpack 入协议）

1. **统一两后端唯一可行的形状。** CUDA 的 GPU-pure `launch()` 不可能按面回调 host 虚函数打包——多 GPU halo 必须从裸设备指针 gather 进 P2P staging buffer；声明式 `(data, width)` 是 CPU/CUDA 都能消费的最大公约数。CPU 的 int pack/unpack 语义由引擎侧实现（halo 引擎用自己已有的每面 send list 做索引打包），即"收进协议"的含义。
2. **零迁移风险。** grep 核实两条旧 API 均无实现者、无调用方，没有任何行为需要移植——删除即净赚。
3. **单点正确性。** PBC 位移、tag 转码、双精度打包（MPI 消息保持 double，默认构建下无损）已集中在 halo_exchange.cpp 的 audited helpers；引擎侧打包让额外字段并入既有每面 `MPI_Isend`（每面一条消息，不随参与者数增殖）。
4. **命令式灵活性无在树需求。** LAMMPS pack_forward 的任意变换能力在本仓库无消费者；将来若需要（如 kspace.cuh:13 预告的 force-infra-04 每原子 halo 钩子），给 CommField 增量扩展即可，不预付成本。

字段精度取每原子存储档（CPU `StorageReal`、CUDA `real_t`，与 CUDA 既有先例一致）；本征 reverse 力通路（`AccReal fx/fy/fz`）硬连线不变。

### 生命周期 / 所有权

协议对象即组件本身（多继承），裸指针非拥有、有效期 build→destroy（与 Plan 持有组件同生命周期）。字段数组由参与组件在 `build()` 分配（owned+ghost 容量，遵守 no-alloc-in-launch），引擎打包进 `HaloSend/HaloRecv` 各自的 `detail::FaceBufs`。绑定时序：`CPUDriver` build 通路在 `install_defaults()` 自动注入**之后**、组件 `build()` 之前，收集全部字段声明非空的参与者，对每个端点调 `bind_participants(...)`——自动注入与用户手注册的端点走同一条路，无 `dynamic_cast`、无 `name()` 匹配（维持 2026-07-02 note 记录的扫描纪律）。绑定注册表**不进** `cpu::Context`（Context 归 05）、不进 Frame（Frame 纯净性）。

### opt-in 协议先例（既有模式）

- `cpu::PairEvaluator`（include/atomiverse/cpu/component/pair_evaluator.h:36-49）——backend-root 协议基类，typed 查询替代 dynamic_cast。
- `cpu::AsyncComponentBase`（include/atomiverse/cpu/async/async_component_base.h:53-88）——launch/destroy `final` + 自有虚集的 opt-in 家族。

`HaloParticipant` 头文件按 pair_evaluator.h 先例放 component/ 支柱目录（协议属组件支柱，消费者在 comm/），两后端同目录（正交性规则）。

### 已否决备选

- **CapabilityProvider / ProfilerParticipant / FramePublisher 三协议套餐**：为不存在的问题造抽象——MC/optimizer/profiling 均不在虚表（profiling 仅 `memory_bytes()` 一虚），拆之无物。
- **命令式 pack/unpack 虚函数搬进协议（LAMMPS 原样）**：CUDA 消费不了 host 虚回调；且为死 API，保留命令式形状等于把分叉也一并搬家。
- **dynamic_cast 探测协议**：违反 scan 无 dynamic_cast/name() 匹配的既有纪律（2026-07-02 note；`cuda::md::MD::install_infra_`，md.cuh:274-291 同纪律）。
- **参与者注册表挂 `cpu::Context` 或 `DomainDecomposition`**：前者踩 05-context 的地盘，后者把组件图谱塞进 comm 状态记录；`bind_participants` 让所有权留在端点组件内。

## Files to create or modify

- include/atomiverse/cpu/component/halo_participant.h (new)
- include/atomiverse/cuda/component/halo_participant.cuh (new)
- include/atomiverse/cpu/component/component.h
- include/atomiverse/cuda/component/component.cuh
- include/atomiverse/cpu/comm/halo_exchange.h
- src/cpu/comm/halo_exchange.cpp
- include/atomiverse/cpu/driver/cpu_driver.h
- tests/mpi/test_halo_exchange.cpp
- tests/cpu/test_install_infra_hook.cpp

## Tasks

- [x] Write failing tests for the HaloParticipant contract: extra-field forward/reverse round-trip via a test-only participant (width=1 and width=3) in tests/mpi/test_halo_exchange.cpp, and protocol-probed endpoint detection + no-double-injection in tests/cpu/test_install_infra_hook.cpp
- [x] Implement cpu::HaloParticipant protocol in include/atomiverse/cpu/component/halo_participant.h (new) — Endpoint, CommField{StorageReal* data, int width}, comm_forward_fields/comm_reverse_fields, bind_participants; Doxygen @brief/@param with units
- [x] Implement cuda::HaloParticipant protocol in include/atomiverse/cuda/component/halo_participant.cuh (new) — member-for-member mirror of the CPU shape with real_t*, CUDA-free header
- [x] Remove the eight halo virtuals from cpu::ComponentBase and add the single halo_participant() query hook in include/atomiverse/cpu/component/component.h
- [x] Remove CommField/comm_forward_fields/comm_reverse_fields from cuda::ComponentBase and add halo_participant() + memory_bytes() (default 0) in include/atomiverse/cuda/component/component.cuh
- [x] Migrate HaloSend/HaloRecv onto HaloParticipant in include/atomiverse/cpu/comm/halo_exchange.h — Endpoint overrides replace provides_halo_send/recv, bind_participants caches the field-participant list
- [x] Implement engine-side packing of declared extra fields into the existing per-face forward (border + position-only refresh) and reverse messages in src/cpu/comm/halo_exchange.cpp, reverse fields accumulated exactly once onto owners
- [x] Migrate CPUDriver InfraScan halo probes to halo_participant()/halo_endpoint() and add the build-time participant-collection + bind_participants pass in include/atomiverse/cpu/driver/cpu_driver.h (InfraScan boolean field names unchanged so md.h/minimize.h consumers stay untouched)
- [x] Verify PE stage-5 MD performance regression [执行移交 ac-008 owed evaluator /mol:bench——本机无集群基线；ledger 保持 pending]: rerun the MSU-PE stage-5 benchmark (benchmarks/md/) on the same host/rank count and gate with benchmarks/md/common/regression_gate.py against a pre-refactor baseline, per .claude/notes/cpu_md_perf_assessment_2026-06.md
- [x] Run full check + test suite (cmake --build build + default ctest gate + tests/mpi multi-rank suite)

## Testing strategy

- **Happy path**：瘦身后全量默认 ctest 门（≈16 s）全绿；tests/mpi 既有多 rank 回归（test_halo_exchange、test_lj_nve、test_minimize_lj、test_pe_* 系列）不改断言直接通过——halo 端点行为与注入时序不变。
- **新协议行为**：测试参与组件（仅测试内定义，不进产品树）声明 width=1 与 width=3 的 forward/reverse 字段，≥2 rank 下验证 (a) 一步后 ghost 行字段值等于对应 owner 值（含跨周期边界面）；(b) reverse 字段的 ghost 贡献恰好累加一次回 owner；(c) 无字段参与者与单 rank（no-op 哨兵 `HaloSend{}` cutoff=0）路径零行为变化。
- **扫描/注入边界**：用户手注册 `HaloSend/HaloRecv` + `set_pair` 时不重复注入（经 Endpoint 探测）；纯字段参与者（Endpoint::none）不被当作端点；`scan_pending_infra_` 无 `dynamic_cast`、无 `name()` 匹配。
- **领域验证**（science.required）：tests/mpi/test_lj_energy_bruteforce、test_bonded_energy_bruteforce、test_coul_energy_bruteforce 多 rank vs 单 rank 能量/力奇偶性维持既有容差——重构不得扰动力的物理值。
- **链级热路径门禁**：plan.h build 时缓存 scheduling()/frequency()（plan.h:95-96）不动；`src/cpu/driver/backend.cpp:243-275` 的 step 虚调度结构不变（backend.cpp:264 timing-off 零开销路径逐字保留）；组件 `launch()` 虚调用形态不变。本篇改动组件基类虚表布局，按链级不变量**必须**跑 PE stage-5 MD 基准回归（benchmarks/md/，参照 .claude/notes/cpu_md_perf_assessment_2026-06.md 的 stage-5 计时法与 benchmarks/md/common/regression_gate.py）。
- **链内独立性**：仅依赖 01-bindings 落地后的 `DefaultBinding` 最终形状（本篇在其上重排基类，避免二次触碰 ~50 个组件头）；不依赖 03-09 任何后续 sub-spec，改动后全树独立可构建、可测试、可合并。保留 `Driver` 命名。

## Out of scope

- **neighbor 三件套**（`needs_neighbor_list`/`provides_neighbor_list`/`neighbor_cutoff_hint`）——归 03-neighbor；本篇原样保留于两侧 ComponentBase。
- **Context/BuildContext 重构**（含任何经 Context 发布参与者注册表的方案）——归 05-context。
- **CUDA 侧声明字段的消费接线**（`cuda::comm::MultiGPU` staging-buffer gather）：在树 CUDA 组件无任何字段声明者，kspace force-infra-04 钩子尚未存在；本篇只统一协议形状，消费接线待首个 CUDA 声明者出现时增量落地（归 04-engines 或后续）。
- **CUDA `memory_bytes()` 的驱动端内存汇总消费**——归 04-engines；本篇只补虚函数声明对齐。
- **CUDA `provides_halo_send/recv` 等注入探测**：CUDA MD 明确不做多 GPU halo 自动注入（md.cuh:69），不为其造探测。
- **文档/叙事更新**（docs/、CLAUDE.md 架构节）——归 09-harness。
- **Python 绑定**——归 08-python。
