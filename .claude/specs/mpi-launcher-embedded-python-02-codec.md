---
title: cpu/task 基础设施集群 — TaskSpec 数据模型（MD + AM1 单点）、无 MPI 字节编解码与任务物化器
status: approved
created: 2026-07-15
revised: 2026-07-15
---

# cpu/task 基础设施集群 — TaskSpec 数据模型（MD + AM1 单点）、无 MPI 字节编解码与任务物化器

## Summary

本 spec 是 `mpi-launcher-embedded-python` 五环链的第 02 环，交付一个新的后端根基础设施集群 `cpu/task/`（对标 `comm/`、`async/`、`replica/` 的"内聚多文件基础设施自成根目录"规则）：一个带版本号、带**任务种类判别子**的内存态 `TaskSpec`，v1 覆盖两类任务——经典 MD 任务（精确枚举今天 Python 绑定 `python/src/bindings/bind_md.cpp` 能表达的九个可绑定组件及其构造参数、显式/默认两条注册路径、步数预算）与 **AM1 半经验单点任务**（一条镜像 `cpu::semiempirical::AM1Config` 的评估器参数、无步数预算）；一对无 MPI 依赖的字节缓冲 pack/unpack 序列化器（版本头 + 任务种类字节校验，未知版本/未知种类/未知组件 id 抛异常）；以及一个 rank 无关的物化器 `TaskMaterializer`，按任务种类分别重建出与直建等价的 `cpu::md::MD` 驱动或 `cpu::Evaluate + AM1Evaluator` 单点管线，并以**一个与任务种类无关的消费入口 `run(frame)`** 供第 03 环的 WorkerLoop 使用。这是纯内部线协议——永不落盘、永不面向用户；Frame 数据不进入 TaskSpec（Frame 经 01-bcast 环独立广播）。AM1 梯度与 AM1 驱动的 MD 明确不在范围内（`AM1Evaluator` 的能力集为 `Energy|Charge|DensityMatrix`，无 `Force`）。

## Domain basis

不适用——本 spec 为纯基础设施（进程内数据模型 + 线协议），不含物理方程；AM1 的物理内容全部在既有 `cpu::semiempirical::NDDOSolver`（本环只搬运其配置参数）。参数单位一律沿用工程单位制（eV、Å、fs、K；AM1 收敛阈值为无量纲密度范数与能量 [Hartree] 容差，随字段记录在 Doxygen 注释中）。

## Design

**任务种类（`task_spec.h`）。** `enum class TaskKind : std::uint8_t { MD = 0, SinglePoint = 1 }`。判别子进入 v1 线格式（链尚未实现，v1 从第一天起就同时包含两类，`kWireVersion` 保持 1 不递增）。不变量由**两个校验构造器**承载（OOP 硬规则——构造期建立不变量，而非事后 validate）：

- `TaskSpec(std::vector<ComponentEntry>, std::int64_t n_steps)` → `kind() == MD`。校验 `n_steps >= 1`，且条目中**不得**出现 evaluator-only 参数（`AM1EvaluatorParams`）——违规抛 `std::invalid_argument` 并点名组件。
- `TaskSpec(AM1EvaluatorParams)` → `kind() == SinglePoint`。恰好一条条目、无显式绑定（binding 恒为 `nullopt`，物化器负责固定绑定，见下）、无步数预算（内部 `n_steps` 固定为 1，不暴露给单点语义）。

解码路径经由这两个构造器重建对象，因此线上层面的违规（SinglePoint 携带 ≠1 条目、携带非 AM1 条目、MD 携带 AM1 条目）同样在 unpack 时抛错。

**数据模型（`task_spec.h`）。** MD 侧逐一对照 `bind_md.cpp` 今天绑定的构造器，每个组件一个参数结构体（普通聚合 + 相等比较，物理校验仍由真实组件构造器负责——不发明组件本身没有的规则）：

| 参数结构体 | 镜像的构造器（bind_md.cpp 现状） |
|---|---|
| `PairLJCutParams{epsilon, sigma, cutoff}` | `cpu::PairLJCut(double, double, double)`（Python 仅绑定单类型构造器；逐类型 `LJTypeCoeff` 构造器未绑定，不入 v1） |
| `NeighborListParams{cutoff, skin}` | `cpu::md::NeighborList<>(cutoff, skin=-1)`（`every/delay/check` 未绑定，不入 v1） |
| `BondHarmonicParams{K, r0}`（`vector<double>` ×2） | `cpu::BondHarmonic` |
| `AngleHarmonicParams{K, theta0}` | `cpu::AngleHarmonic` |
| `DihedralMultiHarmonicParams{Cn}`（`vector<vector<double>>`） | `cpu::DihedralMultiHarmonic` |
| `VelocityVerletParams{dt}` | `cpu::md::VelocityVerlet` |
| `InitVelocityParams{temperature, seed}` | `cpu::md::InitVelocity`（`dist` 未绑定 → 物化时用默认 `gaussian`） |
| `NoseHooverParams{T_start, tau_T, dt}` | `cpu::md::NoseHoover` 三参构造器（`set_T_stop`/`set_run_span` 未绑定，不入 v1） |
| `ThermoParams{freq}` | `cpu::Thermo(int)` |
| **`AM1EvaluatorParams{charge, multiplicity, max_iter, density_tol, energy_tol}`** | **`cpu::semiempirical::AM1Evaluator(AM1Config)`——逐字段镜像 `AM1Config`（= `NDDOConfig`，`cpu/semiempirical/nddo_solver.h`），唯 `method` 不入线：由 variant 备选项本身蕴含 AM1，且 `AM1Solver::normalize_` 无条件钉死它** |

`ComponentParams = std::variant<以上十种>`——`AM1EvaluatorParams` **追加为第 10 备选项**，前九个 `component_id`（由 variant index 派生）保持稳定，故 `kWireVersion = 1` 不动。"不可枚举的组件"在类型系统层面直接不可表达，这就是"rank 0 上 TaskSpec 构造期即抛错"的实现形态。`ComponentEntry{ComponentParams params; std::optional<ExplicitBinding> binding}` 记录注册路径：`nullopt` = `md.add(comp)`；有值 = `md.on(phase).add(comp, order)`。`ExplicitBinding{ATV::md::Phase phase; int order}` 在内存中保持**类型化的** Phase 枚举（stage-typing 规则）；SinglePoint 任务的评估器绑定不进入数据模型——v1 固定为 `evaluate::Phase::evaluate`，由物化器落实，避免为一个无自由度的选择开线格式字段。不存任何 Frame 数据（Frame 纯度，Frame 走 01-bcast）。

**编解码（`task_codec.h`）。** 自由函数 `cpu::task::pack(const TaskSpec&) -> std::vector<std::byte>` 与 `cpu::task::unpack(const std::byte*, std::size_t) -> TaskSpec`——先例是 `comm::pack_exchange`/`unpack_exchange`（现存公开、无 MPI 的序列化器），序列化是边界操作而非组件工厂，不触犯 no-factory 指令。线格式：`{magic 'ATVT', u32 version, u8 task_kind, u32 n_entries, [entry...], i64 n_steps}`；entry = `{u16 component_id, u8 has_binding, (u32 phase_slot, i32 order)?, 逐字段参数}`，vector 以 `u64` 长度前缀。**Phase 编码走 `ATV::slot_of(phase)`**（全库唯一合法的枚举→槽位转换，`test_stage_typing_scrub` 强制）；解码先以 `slot_count<md::Phase>()` 做范围校验再重建枚举，越界抛错。magic/版本不匹配、`task_kind` 越出枚举、`component_id` 越出枚举、缓冲截断均抛 `std::runtime_error`，消息带上违规值本身；种类不变量违规经校验构造器抛 `std::invalid_argument`。同一 MPI 作业内所有 rank 运行同一二进制（第 04 环的单一可执行文件），故采用本机字节序并在头文件注明。整个 `cpu/task/` 不包含 `<mpi.h>`，纯 `tests/` 可测。

**物化器（`task_materializer.h`）。** 类 `TaskMaterializer`（OOP 校验构造器、无 make_* 自由函数）：构造器接收 `const TaskSpec&`，按 `kind()` 分派——

- **MD**：持有 `cpu::md::MD` 成员，按 entry 顺序逐条重放注册（`binding` 为空走 `md_.add(comp)`，有值走 `md_.on(entry.binding->phase).add(comp, entry.binding->order)`），随后 `md_.build()`。刻意不用类目 setter：Python 绑定今天走的就是 `on(phase).add` 路径，逐调用重放才是"与绑定构建出的同一对象"；驱动自身的注册校验原样生效。
- **SinglePoint**：持有 `cpu::Evaluate` 成员，`eval_.on(evaluate::Phase::evaluate).add(semiempirical::AM1Evaluator{cfg})` 后 `eval_.build()`。**不用 `cpu::single_point` 便捷函数**：它在调用内部构造临时 `Evaluate` 并即时 build，把构建推迟到 run 时刻，破坏物化器"构造完成即不变量成立（可直接 run 的驱动）"的既定形态，也不给 03 环 build-once/run-per-task 的余地；显式 `Evaluate + on(evaluate).add` 才是最小且与既有 evaluator 绑定惯例（`AM1Evaluator` 刻意无 `default_bindings()`）一致的路径。

**消费接口（03 环的边界，显式声明）**：唯一与任务种类无关的入口是 `void TaskMaterializer::run(cpu::Frame&)`——MD 内部转发 `md_.run(frame, n_steps_)`，SinglePoint 转发 `eval_.run(frame, 1)`。第 03 环的 WorkerLoop 以 `mat.run(frame)` 消费，**取代**上一版的 `driver()` + `n_steps()` 二元接口（两类驱动是不同的 `CPUDriver<P>` 实例化，没有可共同暴露的驱动类型；预算属于任务种类的内部事务）。另暴露 `TaskKind kind() const` 供协议层记账。rank 无关：不含任何 rank 分支。

**结果契约（单点）。** SinglePoint 运行后 `AM1Evaluator` 写 `frame.scalars.energy` [eV] 与逐原子电荷 `frame.atoms.charge` [e]（能力集 `Energy|Charge|DensityMatrix`）——两者经第 03 环既有的 Frame gather 原路返回：已核实 charge 在 `comm::pack_exchange` 的交换 schema 内（`cpu/comm/atom_pack.h`：x/y/z/v/f/**mass/charge**/type/tag/image），本环无需新增回传通道，只在物化器头文件的 Doxygen 中写明该契约。

**交叉引用。** 并行起草的独立 spec `python-api-namespaces-am1` 把 AM1 绑到 Python；本 spec 与 Python 层无依赖，但 `AM1EvaluatorParams` 与它必须镜像同一个 `AM1Config`——单一事实源是 `cpu/semiempirical/nddo_solver.h` 的 `NDDOConfig`，两边各自对照它，不互相引用。

**构建注册。** 三个 `src/cpu/task/*.cpp` 加入根 `CMakeLists.txt` 的 CPU 源列表（`src/cpu/mc/transaction.cpp` 同段，约 L186）；两个测试经 `tests/CMakeLists.txt` 的 `atv_add_cpp_test` 注册为 `test_cpu_task_codec`、`test_cpu_task_materializer`。物化器测试目标需要 `ATV_BUILD_RHF`-无关的 AM1（半经验不受 RHF 开关约束，随核心库常在）。

## Files to create or modify

- `include/atomiverse/cpu/task/task_spec.h` (new)
- `include/atomiverse/cpu/task/task_codec.h` (new)
- `include/atomiverse/cpu/task/task_materializer.h` (new)
- `src/cpu/task/task_spec.cpp` (new)
- `src/cpu/task/task_codec.cpp` (new)
- `src/cpu/task/task_materializer.cpp` (new)
- `tests/cpu/test_task_codec.cpp` (new)
- `tests/cpu/test_task_materializer.cpp` (new)
- `CMakeLists.txt`
- `tests/CMakeLists.txt`

## Tasks

- [ ] Write failing tests for the TaskSpec codec (MD all-9-component roundtrip equality; SinglePoint AM1 roundtrip equality incl. kind + all five AM1 fields; version mismatch throws; unknown task_kind throws; unknown component id throws; truncated buffer throws; `n_steps < 1` MD ctor throws; MD spec carrying an AM1 entry throws; wire-decoded SinglePoint with ≠1 entries or a non-AM1 entry throws) in `tests/cpu/test_task_codec.cpp`, registered as `test_cpu_task_codec` via `atv_add_cpp_test` in `tests/CMakeLists.txt`
- [ ] Implement the TaskSpec data model (`TaskKind`, ten `*Params` structs incl. `AM1EvaluatorParams` mirroring `AM1Config` + `operator==`, `ComponentParams` variant with AM1 appended as the 10th alternative, `ComponentEntry`/`ExplicitBinding`, the two validating `TaskSpec` ctors, `kWireVersion = 1`) in `include/atomiverse/cpu/task/task_spec.h` + `src/cpu/task/task_spec.cpp`, and register `src/cpu/task/*.cpp` in the root `CMakeLists.txt`
- [ ] Implement `cpu::task::pack`/`unpack` (versioned header + `task_kind` byte, `slot_of`-based phase encoding with `slot_count` range check on decode, decode-through-validating-ctors, value-naming error messages) in `include/atomiverse/cpu/task/task_codec.h` + `src/cpu/task/task_codec.cpp`
- [ ] Write failing tests for `TaskMaterializer` (MD: materialized driver runs 5 steps on a ~10-atom LJ frame via `run(frame)` and matches a directly-built `cpu::md::MD` bitwise on `scalars.energy` and all positions, bonded + thermostat variant included; SinglePoint: materialized AM1 task on the H2 fixture from `tests/cpu/semiempirical/test_am1_scaffold.cpp` matches a directly-built `cpu::Evaluate` + `AM1Evaluator` bitwise on `frame.scalars.energy` and every `atoms.charge` entry) in `tests/cpu/test_task_materializer.cpp`, registered as `test_cpu_task_materializer` in `tests/CMakeLists.txt`
- [ ] Implement `TaskMaterializer` (ctor dispatches on `kind()`: MD replays registrations per entry then `build()`; SinglePoint binds `AM1Evaluator` at `evaluate::Phase::evaluate` on a `cpu::Evaluate` member then `build()`; exposes the task-kind-agnostic `run(cpu::Frame&)` plus `kind()`) in `include/atomiverse/cpu/task/task_materializer.h` + `src/cpu/task/task_materializer.cpp`
- [ ] Add docstring per doxygen with units (`@brief`/`@param` [eV, Å, fs, K; AM1 tolerances dimensionless / Hartree]/`@return` on every public symbol in `cpu/task/`, incl. the single-point results contract: energy [eV] on `scalars.energy`, charges [e] on `atoms.charge`)
- [ ] Verify against the fidelity contract: field-by-field diff of the materializer's replay sequence vs `bind_md.cpp`'s bindable set, and of `AM1EvaluatorParams` vs `NDDOConfig` in `include/atomiverse/cpu/semiempirical/nddo_solver.h` (every ctor arg represented, nothing unbindable representable, `method` deliberately absent)
- [ ] Run full check + test suite

## Testing strategy

- **正常路径（`test_cpu_task_codec`）：** (a) MD 任务——含全部九种 MD 组件（三个默认绑定条目 + 六个显式 `on(phase).add(comp, order)` 条目，覆盖非零 order 与至少两个 `md::Phase` 值）与 `n_steps`，pack→unpack 后逐字段相等（含 vector 参数、`optional` 绑定、phase、order、n_steps、kind）。(b) SinglePoint 任务——非默认的 `AM1EvaluatorParams`（如 charge=-1、multiplicity=2、max_iter=50、非默认双容差）往返后 kind 与全部五字段逐一相等。
- **边界/异常：** 版本字段被篡改 → 抛错且消息含违规版本号；`task_kind` 字节越出枚举 → 抛错且消息含违规值；`component_id` 越出枚举 → 抛错且消息含违规 id；缓冲截断 → 抛错；MD 构造 `n_steps = 0` → 抛 `std::invalid_argument`；MD 条目列表混入 `AM1EvaluatorParams` → 构造期抛错点名组件；手工构造的线缓冲让 SinglePoint 携带两条条目或一条非 AM1 条目 → unpack 抛错（解码经校验构造器）；空 entry 列表的 MD 任务可编解码（物化时由驱动自身的组件校验报错，不在 codec 层拦截）。
- **物化器等价（`test_cpu_task_materializer`）：** (a) MD——~10 原子 LJ 小体系（PairLJCut + NeighborList + VelocityVerlet + InitVelocity 固定种子），物化驱动经 `run(frame)` 与手写直建驱动各跑 5 步，`frame.scalars.energy` 与全部 x/y/z **逐位相等**；第二用例加 BondHarmonic + NoseHoover。(b) AM1 单点——H2 定长分子（复用 `test_am1_scaffold.cpp` 的 fixture 形态），物化的 SinglePoint 任务 vs 直建 `cpu::Evaluate` + `on(evaluate::Phase::evaluate).add(AM1Evaluator{同 cfg})`，`frame.scalars.energy` 与全部 `atoms.charge` 逐位相等。**关于 tests/ 纪律的诚实说明**：H2 的 AM1 SCF 在既有单元测试中即完整收敛（`tests/cpu/semiempirical/test_am1_scaffold.cpp` 已对 H2/CH4/甲醇跑 `AM1Solver::run` 到收敛，H2 ≤4 次迭代、毫秒级）——微型 AM1 求解在树内的地位与几步小体系 LJ MD 相同，是既定规范，不属于三层铁律针对的"完整 SCF/MD 求解"（那针对的是 Gaussian/PW/GPW 量级与平衡化跑长）；故等价检查留在 `tests/`，不下放 04 环回归。
- **既有门禁回归：** `test_stage_typing_scrub` 必须在新增 codec 源上保持通过（phase 编码只走 `slot_of`，无手写 `static_cast<int>`）。

## Out of scope

- MPI 传输本身：Frame 广播是 01-bcast 环，`RUN_TASK`/`SHUTDOWN` 协议标签、WorkerLoop 与 rank-0 派发是 03-protocol 环；本环产物只是它们搬运的字节负载。
- **AM1 梯度与 AM1 驱动的 MD/优化**：`AM1Evaluator` 无 `Capability::Force`，SinglePoint 之外的 AM1 任务在类型上不可表达；将来 AM1 出力时随 `kWireVersion` 递增进入。
- 其他方法的评估器参数（SCF/PW/GPW/RHF）：v1 单点只收 AM1；扩张同样走版本递增。
- Python 绑定钩子、`app/atomiverse/` 可执行文件与 CPython 嵌入（04-embed）、wheel 打包（05-packaging）；AM1 的 Python 命名空间绑定属于并行 spec `python-api-namespaces-am1`，与本环互不依赖（共同镜像 `NDDOConfig`）。
- 今天绑定集之外的任何 MD 组件或参数：`NoseHooverNPT`/`FixLangevin`/`Barostat`/`Deform`、`PairLJCut` 逐类型构造器、`NeighborList` 的 `every/delay/check`、`InitVelocity` 的 `dist`、`NoseHoover` 的 `set_T_stop`/`set_run_span`。
- 任何落盘格式、版本协商/向后兼容解码（单版本、不匹配即抛）、跨字节序可移植性（同作业同二进制）。
- CUDA/torch 后端的 task 集群（本链全程 CPU 后端）。
