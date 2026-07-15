---
title: "MPI task protocol: RUN_TASK/SHUTDOWN + WorkerLoop + TaskDispatcher (mpi-launcher-embedded-python 03)"
status: approved
created: 2026-07-15
revised: 2026-07-15
---

# MPI task protocol: RUN_TASK/SHUTDOWN + WorkerLoop + TaskDispatcher（链 03-protocol）

## Summary

本篇是 `mpi-launcher-embedded-python` 五链中的第 03 环：在 `cpu/task/` 基础设施层落地纯原生、零 Python 的任务粒度 MPI 协议。rank 0 上的 `TaskDispatcher` 在每次任务提交时做**一次**任务粒度广播（TaskSpec 字节走链 02 的 codec，Frame 走链 01 的 `cpu::comm::bcast_frame`），随后**全部 rank（含 rank 0）**进入既有的原生 MPI 计算路径（02 的 `TaskMaterializer` 构造即物化+build，经其任务类型无关入口 `run(frame)` 执行——MD 与 AM1 单点同一条协议路径，库内 `Backend` 照常探测 `MPI_Initialized` 并做域分解）；跑完后 rank 0 按 `gather_to_root` 纪律收回结果 Frame，worker rank 回到 `WorkerLoop` 的阻塞等待，直到收到 SHUTDOWN 才退出。协议只有 RUN_TASK 与 SHUTDOWN 两个动词（用户明确否决 CREATE_SYSTEM/CREATE_MODEL 这类细粒度同步），带版本握手护栏，错误经全归约状态传播而非死锁，同一会话内可顺序派发多个任务。整环可从 `tests/mpi/` 的原生 main 完整测试；链 04 的可执行程序与 CPython 嵌入是本协议的消费者。

## Design

**前置依赖（链内契约，接口名以落地头文件为准）。** 本篇消费链 01 的 `cpu::comm::bcast_frame`（整 Frame 广播）与链 02 的 `cpu/task/` TaskSpec + pack/unpack codec + MD materializer（TaskSpec → driver + 步数预算）。注意：当前 `dev` 工作树中尚无这两环的产物（`include/atomiverse/cpu/task/` 目录与 `bcast_frame` 符号均不存在）——实施本篇前必须先合入 01/02；实现时以它们落地的实际签名为准，不得在本篇内重造。

**线协议（wire protocol）。** 所有集体操作走构造时传入的 communicator，root 固定为 rank 0：

- `enum class Command : std::int32_t { RUN_TASK = 1, SHUTDOWN = 2 };`（从 1 起，零初始化缓冲不构成合法命令）。
- `constexpr std::int32_t kProtocolVersion = 1;`
- 每条命令一个定长头 `CommandHeader{ std::int32_t version; std::int32_t command; std::int64_t payload_bytes; }`，一次 `MPI_Bcast(MPI_BYTE)`。同一二进制、同一 ABI 是唯一受支持模式（版本握手即为此的廉价护栏）。
- **RUN_TASK 序列**（每个 rank 严格同序参加每个集体调用）：
  1. Bcast `CommandHeader{version, RUN_TASK, spec_bytes}`；
  2. Bcast TaskSpec 字节（02 codec 的 pack 产物）；
  3. `cpu::comm::bcast_frame`（链 01）；
  4. 全 rank 走**同一条**执行路径：unpack → 构造 02 的 `TaskMaterializer`（构造即物化 + build），本地异常一律捕获、折算为本地状态码，**不得**在集体点之间提前逃逸；
  5. **setup 状态全归约**（`MPI_Allreduce`, MAX）：任一 rank 非零 → 所有 rank 同步跳过 run 与 gather，rank 0 的 `dispatch()` 抛出并携带原因，worker 回到阻塞——锁步保持；
  6. 全零则 `materializer.run(frame)`（02 的任务类型无关唯一入口；步数预算/单点语义是任务种类的内部事务，协议层无 per-kind 分支），本地异常同样捕获；
  7. **run 状态全归约**：全零才进入第 8 步；
  8. `gather_to_root(frame, dd, 0)`（`cpu/comm/gather.h`）→ rank 0 上 `to_frame(gathered)` 拼出结果 Frame；scalars（energy/kinetic 等）取 rank 0 本地 post-run frame（`Reduce` 已在 `end_of_step` 把全局归约值写回每个 rank 的 `frame.scalars`）。若现行 `to_frame` 未回填 `image` 列，在本篇内补齐（`GatheredFrame` 已携带该数据）。
- **SHUTDOWN 序列**：仅头广播；`WorkerLoop::run()` 返回，调用方（04 的 main）随后自行 `MPI_Finalize`。

**版本握手。** `TaskDispatcher`（rank 0）与 `WorkerLoop`（rank > 0）的构造函数对称地参加同一集体握手：rank 0 Bcast 自己的协议版本，各 rank 与本地版本比较，比较结果经一次 `MPI_Allreduce` 汇总；不匹配 → **所有 rank** 抛 `std::runtime_error`，错误信息点名两个版本号（Allreduce 是握手的最后一个集体点，抛出后仍锁步）。构造函数带 `std::int32_t version = kProtocolVersion` 参数作为测试注入缝。每条命令头里的 `version` 字段是纵深防御：worker 侧不匹配视为协议损坏，致命退出。

**两个类（OOP 校验构造，无 make_*/工厂自由函数）。**

- `cpu::task::WorkerLoop`：以 communicator 构造（校验：MPI 已初始化、comm 非 `MPI_COMM_NULL`、本 rank ≠ 0，否则抛出）；构造即参加版本握手。`run()` 阻塞循环：收头 → RUN_TASK 走上述共享执行路径后**回到阻塞**（任务级失败不退出会话，只有协议级损坏才致命）→ SHUTDOWN 时干净返回。SHUTDOWN 之前绝不自行退出。
- `cpu::task::TaskDispatcher`：以 communicator 构造（校验：本 rank == 0）；`dispatch(const TaskSpec&, Frame)` 广播 + 本地参与同一执行路径 + gather + 返回结果 Frame（**Frame 是唯一进出货币**）；`shutdown()` 广播 SHUTDOWN。
- 头广播之后的执行路径（unpack → materialize → build → run → 双状态归约 → gather）由 `cpu::task::internal` 下的一个共享入口承载，恰好两个调用方（dispatcher 与 worker）——符合"第二个调用者才提取"的规则。

**关停顺序契约。** 两个类都**不拥有** communicator（不 dup、不 free），析构函数中**没有任何 MPI 调用**；每个任务的 driver 都是执行路径内的局部对象，run 返回即析构（`Backend::shutdown()` → `dd.release_comm()`，见 `cpu/comm/domain.h:181-191`），worker 回到阻塞时不残留任何 cart_comm。因此 04 的 main 持有的 `MPI_Finalize` 顺序契约（Backend 先释放 cart_comm、async executor 探测 `MPI_Finalized`）由构造上保证可满足，本层不新增任何 finalize 时序义务。

**锁步纪律的边界（如实声明）。** 本协议保证**协议级**失败（版本不匹配、TaskSpec unpack 失败、materialize/build 抛出）锁步且经状态归约传播；`driver.run()` **中途**的集体调用内部失败（今天的原生 MPI 跑一样会遇到）不在本协议的恢复范围内——与现状持平，不倒退。

**单 rank / 非 MPI 构建。** comm size == 1 时无 worker 存在，`dispatch()` 退化为本地直接执行（广播与 gather 均为平凡路径）；非 MPI 构建下语义与单 rank 相同（对齐 `gather_to_root` 的既有约定）。`WorkerLoop` 在 rank 0 或单 rank 下构造即抛（校验构造，无回退值）。

## Files to create or modify

- `include/atomiverse/cpu/task/task_protocol.h` (new) — Command / CommandHeader / kProtocolVersion / WorkerLoop / TaskDispatcher 声明（目录本身由链 02 创建）
- `src/cpu/task/task_protocol.cpp` (new) — 握手、命令环、共享执行路径、dispatch/gather 实现
- `tests/mpi/test_task_protocol.cpp` (new) — 五场景 MPI 测试（自带 MPI main，gtest 惯例同 `tests/mpi/test_lj_nve.cpp`）
- `tests/CMakeLists.txt` — `ATV_USE_MPI` 块内注册 `atv_add_mpi_exec(test_mpi_task_protocol …)` + n=1/2/4 的 `atv_add_mpi_run`
- `CMakeLists.txt` — 根源文件清单加入 `src/cpu/task/task_protocol.cpp`（比照 225–234 行 comm 源的显式列法）

## Tasks

- [ ] Write failing tests for the task protocol (tests/mpi/test_task_protocol.cpp (new)): (a) 派发一个 ~10 原子 LJ MD 小任务（数步、dt 1 fs 量级），rank 0 结果 Frame 与同一 TaskSpec 不经协议、全 rank 直跑 + gather 的参照逐字段比对（位置/速度容差 1e-10，能量相对 1e-10；n=1 变体即"与串行跑等价"）；(b) 干净 SHUTDOWN，所有 rank 返回并以 0 退出；(c) 同一会话顺序两个任务各自匹配参照；(d) 版本不匹配（构造注入 version+1）→ 所有 rank 抛 `std::runtime_error` 且不死锁；(e) 非法 TaskSpec → `dispatch()` 在 rank 0 抛出、worker 存活、随后一个合法任务成功、干净关停
- [ ] Register the MPI test in tests/CMakeLists.txt (`atv_add_mpi_exec` + `atv_add_mpi_run` n=1/2/4, `ATV_USE_MPI` 块内)
- [ ] Implement Command, CommandHeader, kProtocolVersion and the WorkerLoop / TaskDispatcher declarations in include/atomiverse/cpu/task/task_protocol.h (new)，全部公有符号带 doxygen `@brief`/`@param`（含单位：位置 Å、速度 Å/fs、能量 eV、预算 steps）/`@return`
- [ ] Implement the shared post-broadcast execution path (unpack → materialize → build → run → 双状态 Allreduce → gather_to_root/to_frame) in src/cpu/task/task_protocol.cpp (new)，本地异常一律折算状态码、集体调用严格同序；若 `to_frame` 缺 `image` 回填则在 src/cpu/comm/gather.cpp 一并补齐
- [ ] Implement WorkerLoop in src/cpu/task/task_protocol.cpp：校验构造 + 握手，`run()` 阻塞命令环，RUN_TASK/任务级错误后回到阻塞，SHUTDOWN 才返回，析构无 MPI 调用
- [ ] Implement TaskDispatcher in src/cpu/task/task_protocol.cpp：校验构造 + 握手，`dispatch(TaskSpec, Frame)` 一次任务粒度广播（header → spec bytes → bcast_frame）+ 本地同路径执行 + 返回结果 Frame（scalars 取 rank 0 post-run frame），`shutdown()` 广播 SHUTDOWN
- [ ] Add src/cpu/task/task_protocol.cpp to the root CMakeLists.txt source list
- [ ] Verify lockstep and equivalence at n=1/2/4 via `ctest -L mpi`（含错误路径不超时、退出码为 0）
- [ ] Run full check + test suite

## Testing strategy

- **happy path**：任务粒度派发 ≡ 同一 TaskSpec 的直跑 MPI 计算（协议只增删广播/物化/收集，不得改变物理结果）；n=1 变体即串行等价，n=2/4 变体覆盖真实分解。系统保持 ~10 原子、数步，秒级完成（`tests/mpi/` 既有 mpirun ctest 基建：`mpi` label、`RESOURCE_LOCK mpi_singleton`）。
- **边界**：同一会话两个顺序任务（worker 回环）；干净 SHUTDOWN 全员退出 0（同时佐证关停顺序契约——finalize 处无 cart_comm 残留、无析构集体调用）；每个 TEST 自成一个"握手→任务→关停"会话，多个会话串行也随之被覆盖。
- **错误路径锁步**：版本不匹配（不对称注入，恰好检验非对称失败下的锁步）→ 所有 rank 抛出、握手后仍对齐；TaskSpec unpack/materialize 失败 → setup 状态归约令全员同步跳过 run/gather，rank 0 得到非零传播，会话存活可继续派发。
- **域验证**：本篇为纯基础设施，无新物理；物理正确性由"派发 ≡ 直跑"这一等价性检查承载。链的唯一 regression 案例按裁决随 04 交付，本篇不加 regression。

## Out of scope

- 链 04（`app/atomiverse/` 可执行 + rank 0 CPython 嵌入）与链 05（wheel 打包）——本篇只保证纯原生可测。
- 细粒度对象同步动词（CREATE_SYSTEM/CREATE_MODEL 等）——用户已明确否决，协议永远只有 RUN_TASK/SHUTDOWN。
- `driver.run()` 中途集体调用失败的容错/恢复（与现状原生 MPI 跑持平）。
- 跨二进制异构部署支持（版本握手只做守门，不做协商降级）。
- communicator 的 dup/隔离与多会话并发；本篇是单会话、顺序派发。
- regression 案例（链裁决：随 04 一并交付）。
- 链 01/02 接口的任何改动——本篇是消费者；若其签名不敷使用，回到对应链修，不在此篇内绕。
