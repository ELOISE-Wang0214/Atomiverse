---
title: MPI launcher — rank-0 embedded CPython executable (chain 04-embed)
status: approved
created: 2026-07-15
revised: 2026-07-15
---

# MPI launcher — rank-0 embedded CPython executable (chain 04-embed)

## Summary

本 spec 是已批准的 5 环链 `mpi-launcher-embedded-python` 的第 4 环：新增顶层 `app/atomiverse/` 可执行程序。用户以 `mpirun -np N atomiverse atv.py` 启动：全部 rank 由 `main.cpp` 统一 `MPI_Init_thread(MPI_THREAD_MULTIPLE)`；rank 0 嵌入 CPython 并以 `__main__` 身份执行 `atv.py`，rank 1..N-1 进入第 03 环的 `WorkerLoop`、完全不触碰 Python。分发缝（dispatch seam）落在 `python-api-namespaces-am1` 重构后的绑定层，覆盖**两种任务**：脚本中每次 `MD.run()`（TaskKind::MD）与每次 AM1 单点调用（TaskKind::SinglePoint，携 AM1 参数），在 launcher 托管且 world size>1 时经由第 02 环（本日并行修订版）的 `TaskSpec` + 第 03 环的 `TaskDispatcher` 做一次任务粒度广播（TaskSpec + Frame，复用第 01 环 `cpu::comm::bcast_frame`），无任何落盘 IR；能量与逐原子电荷经既有 Frame gather 通道返回。普通 `python atv.py` 串行场景下该缝是逐字节等价的 no-op；本次修订同时修复评审指出的子进程继承缺陷——`LauncherGate` 的标记改为 **PID 限定**且"不一致即降级为 inactive、绝不抛出"。脚本正常结束或抛异常时都向 worker 广播 SHUTDOWN、把非零退出码传给 mpirun。交付整条链**唯一**的 regression 案例 `regressions/launcher/lj_parity/`（LJ MD 串并行 parity），其编排脚本同案内新增 AM1 单点 parity 腿。

## Domain basis

本特性是**进程模型**（启动、嵌入、分发、关停），不引入任何新物理。regressions/ 的"必须对照有出处的文献值"规则针对的是*物理结果*——它防止的是用我们自己的输出做 golden 快照来自证。而本链可观测的契约本身就是**串行/并行等价**：launcher 分发路径的正确性定义即"与串行 python 托管的同一脚本产生相同的物理结果"。两条腿：

- **LJ MD parity**（本案例的主体）：同一 `atv_lj.py`，`python` 串行托管 vs `mpirun -np 2 atomiverse` 托管，能量与位置紧公差一致；串行腿的 LJ MD 物理已由既有 `regressions/md/lj_argon`（文献值）背书——parity 案例不重复验物理，只验进程模型。公差取双精度下 2-rank 域分解求和序差异的量级：逐采样点总能量相对偏差 ≤ 1e-9，200 步后位置偏差 ≤ 1e-8 Å。
- **AM1 单点 parity 腿**（同案例编排脚本的附加腿）：同一 `atv_am1.py` 串行 vs 2-rank launcher 托管。AM1 单点是**确定性**计算（无轨迹混沌、无积分序放大），故公差取到收敛判据地板：总能量相对偏差 ≤ 1e-10，逐原子电荷绝对偏差 ≤ 1e-10 e。AM1 物理本身由 semiempirical 特性自己的 regression 背书，此腿仍只验进程模型。

**为何折叠进同一个 regression 案例而不另立一案**：一特性一案例的铁律以*交付的特性*为单位，而本环交付的特性是 **launcher 这一个进程模型**——MD 与 SinglePoint 是它分发的两种任务负载，不是两个特性。两腿共享同一可执行文件、同一 gate、同一 SHUTDOWN 路径与同一编排脚本；拆成两案会把一个特性的验收拆成两份重复的进程编排。

**运行环境（regression 依赖，必须显式声明）**：Python 侧全部工作与测试在 venv `~/work/.x86_64`（x86_64 登录节点）中进行；该环境须能同时 `import atomiverse`（`ATV_BUILD_PYTHON=ON` 构建产物在 `PYTHONPATH` 上）与 `import molpy`（`~/work/molcrafts/molpy`）。案例另需 `mpirun` 与本环的 `atomiverse` 可执行文件；以 `module load buildenv-gcc/2026.03-mpich` + `export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"` 运行，整案 ≤ 60 s。**解释器发现**：本环显式假定构建期 Python（CMake `FindPython` 找到的解释器与其默认 `PyConfig` home）；wheel 内的解释器发现/重定位推迟到第 05 环。

## Design

**前置依赖（诚实声明）**：本环消费四个先行 spec 的产物——`cpu::comm::bcast_frame`（01）、`cpu/task/` 下的 `TaskSpec` + codec + 物化器（02，**本日并行修订**：`TaskSpec` 增加 `TaskKind {MD, SinglePoint}` 判别子 + `AM1EvaluatorParams`（镜像 `cpu::semiempirical::AM1Config` = `NDDOConfig`：`charge, multiplicity, max_iter, density_tol, energy_tol`，method 不入线），物化器接口改为任务类型无关的 `run(frame)`）、`task_protocol.h` 的 RUN_TASK/SHUTDOWN + `WorkerLoop` + `TaskDispatcher`（03）、以及 **`python-api-namespaces-am1`**（Python 包重构为深层规范命名空间 `atv.potentials.pair.LJCut` / `atv.md.integrator.VelocityVerlet` / `atv.semiempirical.AM1Evaluator` + 顶层再导出，删除旧平铺布局，新增 `bind_semiempirical.cpp` 与可从 Python 调用的 AM1 单点入口）。**这些符号在当前 dev checkout 中尚未全部出现**；本环实现必须在上述四者合入 dev 之后启动，并以合入后的真实头文件/绑定文件路径为准。

**`app/atomiverse/` 可执行程序（新顶层）**：
- `main.cpp` 拥有 `MPI_Init_thread(MPI_THREAD_MULTIPLE)` / `MPI_Finalize` 括号——与库的 adopt-only 契约一致（库只探测 `MPI_Initialized`：`src/cpu/driver/backend.cpp:82-89`），与 `tests/mpi/*.cpp` 的既有先例一致。流程：`LauncherArgs` 校验 argv（脚本缺失/不可读 → 抛出；mpirun 下各 rank argv 相同，校验确定性一致，故所有 rank 同步走非零退出、干净 `MPI_Finalize`，不挂起）→ rank 0 构造 `PythonHost` 运行脚本、**无条件**（正常返回或脚本异常皆然）广播 SHUTDOWN → 其余 rank 运行 `WorkerLoop` 直到收到 SHUTDOWN → 全体 `MPI_Finalize`，rank 0 的脚本退出码传给 mpirun。`Backend::shutdown()`-先于-`MPI_Finalize` 的顺序（`cpu/comm/domain.h:184`）自动满足：驱动器析构发生在 `Py_RunMain` 内部的解释器终结阶段，早于 SHUTDOWN 广播与 `MPI_Finalize`。
- `LauncherArgs`（header-only，校验构造函数）：`argc/argv` → 脚本路径；缺参、路径不存在、不可读均抛 `std::runtime_error`。不是 `parse_X()` 自由函数——是构造即校验的类。
- `PythonHost`（校验构造函数持 `LauncherArgs`）：在 `Py_InitializeFromConfig` **之前** `setenv("ATV_LAUNCHER", <本进程 PID 十进制>)`（见下），然后 `PyConfig` + `config.run_filename = script`，`Py_RunMain()` 以 `__main__` 执行脚本并返回退出码——脚本异常时 CPython 自行打印 traceback 并返回非零。仅此可执行目标链接 `libpython`（`find_package(Python 3.9 COMPONENTS Interpreter Development.Embed)`）；`libatomiverse` 与 `_atomiverse.so` 不新增任何链接项，wheel 的 `python atv.py` 串行路径保持逐字节不变。

**分发缝（为何不能是 C++ 全局钩子，为何不会泄漏进串行）**：`atomiverse` 静态库同时被 launcher 可执行文件与 `_atomiverse.so` 各自链接一份——launcher 一侧写入的 C++ 全局变量对扩展模块一侧**不可见**。真正跨越"可执行文件 ↔ 扩展模块"边界的进程级状态只有两种：进程环境变量，与共享的 `libmpi.so` 状态。缝仍为库内新类 `cpu::task::LauncherGate`（`include/atomiverse/cpu/task/launcher_gate.h`），但契约**本次修订重设计**：

- **标记 PID 限定**：`PythonHost` 写入的标记值是 rank 0 launcher 进程自己的 PID（唯一写入者，且在 `Py_Initialize` 之前）。`active()` ⇔ 标记存在 ∧ 标记值 == `getpid()` ∧ `MPI_Initialized` ∧ size>1。
- **不一致一律降级，绝不抛出**：标记存在但 PID 不匹配 → **静默 inactive**（这正是子进程形态：用户脚本 spawn 的子 python 继承了环境变量但拥有新 PID，`import atomiverse` 后照常走串行路径，零噪声、零抛出——修复评审指出的"子进程见 标记∧MPI 未初始化 即抛"缺陷）；标记 PID 匹配 ∧ MPI 未初始化 → inactive + **一次性** stderr 警告（`std::once_flag`；宿主状态真不一致，launcher 主流程中不可能出现，保留可诊断性）；size==1 → 静默 inactive。
- **为何选 PID 限定而不是 unsetenv 或裸警告**：`PythonHost` 无法在设标记之后再 `unsetenv`——gate 在扩展模块一侧于首次 `run()` 时才惰性构造读取环境，`PythonHost`（可执行文件一侧的另一份库副本）不可能知道扩展侧何时读完，过早 unset 会让 gate 永远读不到标记；而无 PID 的裸"警告 + inactive"会让使用 multiprocessing 池的用户脚本每个 worker 打一条警告，噪声正比于池大小。PID 限定使常见的子进程路径完全静默，把警告留给真正的宿主不一致。

普通 `python atv.py` 下标记必然不存在，代价是一次返回 `nullptr` 的 `getenv`，其后走与修改前完全相同的指令路径；mpi4py 托管的普通 python（本链 out of scope）同样因无标记而保持串行路径。

**绑定层改动（`python-api-namespaces-am1` 重构后的绑定层，两种任务）**：
- **MD（`python/src/bindings/bind_md.cpp`，重构后路径为准）**：诚实声明——绑定的驱动器持有从 Python 构造出的 C++ 组件对象，分发路径**不**内省已 build 的组件，而是在注册时同步记录：`nb::class_` 绑定对象改为 binding 内部 holder（`struct MdHandle { cpu::md::MD driver; <02 的 TaskSpec 镜像>; }`，Python 可见类名与全部方法签名不变，重构后的 Python 包层零改动），每个 `add_*` lambda 在转发给驱动器的同时向镜像追加对应 codec 条目。`run()` 中：`LauncherGate` 不活跃 → 原路径原样执行（含既有 `gil_scoped_release`）；活跃 → rank 0 将镜像 `TaskSpec`（TaskKind::MD）+ 借入的 molrs.Frame 交给 03 的 `TaskDispatcher`（GIL 同样在计算窗口释放），结果写回同一 molrs.Frame。凡 02 codec 不覆盖的已注册组件，在 launcher 下 `run()` 即抛 `std::runtime_error` 并点名组件类型（无回退）；串行下照常工作。
- **AM1 单点（`python/src/bindings/bind_semiempirical.cpp`，由 `python-api-namespaces-am1` 创建）**：该 spec 绑定的 AM1 单点入口（Evaluate 包装或 single_point 函数，以其落地形态为准）加同型分支：gate 不活跃 → 原串行单点路径；活跃 → rank 0 打包 TaskKind::SinglePoint 的 `TaskSpec`（`AM1EvaluatorParams` 镜像该入口收到的 `AM1Config`）+ Frame 交给 `TaskDispatcher`，所有 rank 经任务类型无关的物化器走 Evaluate 路径，能量与逐原子电荷经既有 Frame gather 通道返回并写回借入的 molrs.Frame。同一 gate、同一广播、同一 SHUTDOWN——AM1 不引入第二套协议。AM1 梯度 / AM1-MD 不在本环（见 Out of scope）。

**构建接线**：根 `CMakeLists.txt` 新增 `option(ATV_BUILD_LAUNCHER ... OFF)`；`ATV_BUILD_LAUNCHER=ON` 而 `ATV_USE_MPI=OFF` 或 `ATV_BUILD_PYTHON=OFF` 时 configure 期 `FATAL_ERROR`（沿用 torch 选项的失败先例，绝不静默产出不能工作的配置）；`add_subdirectory(app/atomiverse)`。**Target 命名**：库 target 已占用 `atomiverse` 名（根 `CMakeLists.txt:312`），可执行 CMake target 命名为 `atomiverse_launcher` 并设 `OUTPUT_NAME atomiverse`——磁盘上的二进制仍叫 `atomiverse`（与 05-packaging 的 install 规则对齐）。regression 接入 `regressions/CMakeLists.txt` 既有 `regressions` 自定义 target，但因需要 `mpirun` + venv Python 环境而不能复用 `atv_add_regression()`——单独的 custom target 执行 `reg_lj_parity.sh`（LJ 串行腿 → LJ 2-rank 腿 → AM1 串行腿 → AM1 2-rank 腿 → `compare.py` 判两组公差），且仅在 `ATV_BUILD_LAUNCHER=ON` 时注册。

## Files to create or modify

- `CMakeLists.txt` — 新增 `ATV_BUILD_LAUNCHER` 选项、配置期一致性 FATAL_ERROR、`add_subdirectory(app/atomiverse)`
- `app/atomiverse/CMakeLists.txt` (new) — 可执行目标 `atomiverse_launcher`（`OUTPUT_NAME atomiverse`），`Development.Embed`，仅此目标链接 libpython
- `app/atomiverse/main.cpp` (new)
- `app/atomiverse/launcher_args.h` (new) — header-only `LauncherArgs`
- `app/atomiverse/python_host.h` (new)
- `app/atomiverse/python_host.cpp` (new)
- `include/atomiverse/cpu/task/launcher_gate.h` (new)
- `src/cpu/task/launcher_gate.cpp` (new)
- `python/src/bindings/bind_md.cpp` — MdHandle holder + TaskSpec 镜像记录 + `run()` 分发分支（以 `python-api-namespaces-am1` 合入后的实际路径为准）
- `python/src/bindings/bind_semiempirical.cpp` — AM1 单点入口的分发分支（该文件由 `python-api-namespaces-am1` 创建）
- `tests/CMakeLists.txt` — 接入两个新单元测试（`launcher_gate.cpp` 一并加入库源列表）
- `tests/cpu/task/test_launcher_gate.cpp` (new)
- `tests/app/test_launcher_args.cpp` (new)
- `regressions/CMakeLists.txt` — 接入 launcher parity 案例（`ATV_BUILD_LAUNCHER` 门控）
- `regressions/launcher/lj_parity/atv_lj.py` (new) — 小 LJ MD 脚本（固定种子，分段 run 采样能量/位置）
- `regressions/launcher/lj_parity/atv_am1.py` (new) — 小分子 AM1 单点脚本（输出能量 + 逐原子电荷）
- `regressions/launcher/lj_parity/compare.py` (new) — 双腿公差比较器（MD 能量/位置 + AM1 能量/电荷）
- `regressions/launcher/lj_parity/reg_lj_parity.sh` (new) — 四腿编排脚本（LJ 串行/2-rank + AM1 串行/2-rank + 比较）

## Tasks

- [ ] Write failing tests for LauncherGate and LauncherArgs (tests/cpu/task/test_launcher_gate.cpp: 无标记→inactive；标记为**他进程 PID** ∧ MPI 未初始化→静默 inactive、**不抛**（子进程形态）；标记为**本进程 PID** ∧ MPI 未初始化→inactive + 一次性警告、**不抛**；tests/app/test_launcher_args.cpp: 缺参/不存在/不可读→throw、可读文件→ok；接入 tests/CMakeLists.txt) *(rework: gate 契约从"宿主不一致抛出"改为"PID 限定 + 降级不抛"，测试重写)*
- [ ] Implement LauncherGate in include/atomiverse/cpu/task/launcher_gate.h + src/cpu/task/launcher_gate.cpp (校验构造函数：读 ATV_LAUNCHER 并解析 PID + MPI_Initialized/size 探测；active() ⇔ 标记∧PID 匹配∧MPI 已初始化∧size>1；一切不一致降级 inactive，PID 匹配而 MPI 未初始化时打一次性 stderr 警告，绝不抛出) *(rework: 同上契约变更)*
- [ ] Implement LauncherArgs in app/atomiverse/launcher_args.h (header-only 校验构造函数)
- [ ] Implement PythonHost in app/atomiverse/python_host.h + python_host.cpp (Py_Initialize 前 setenv ATV_LAUNCHER=<本进程 PID>；PyConfig + Py_RunMain 以 __main__ 执行脚本；返回退出码，异常时 traceback 由 CPython 打印) *(rework: 标记值从 "1" 改为 PID)*
- [ ] Implement main.cpp in app/atomiverse/ (MPI_Init_thread(MPI_THREAD_MULTIPLE)/MPI_Finalize 括号；rank 0 = PythonHost + 无条件 SHUTDOWN 广播；rank>0 = 03 的 WorkerLoop；用法错误所有 rank 干净非零退出)
- [ ] Add ATV_BUILD_LAUNCHER wiring: 根 CMakeLists.txt 选项 + FATAL_ERROR 守卫 + add_subdirectory；编写 app/atomiverse/CMakeLists.txt（target `atomiverse_launcher`、`OUTPUT_NAME atomiverse`，仅可执行目标链接 Python3 Embed，libatomiverse 与 _atomiverse.so 不变）
- [ ] Implement dispatch seam for both task kinds in python/src/bindings/bind_md.cpp and bind_semiempirical.cpp (MD：MdHandle holder + 注册期 TaskSpec 镜像，run() 中 gate 活跃→TaskKind::MD 经 TaskDispatcher 分发（释放 GIL），codec 不覆盖的组件抛出并点名；AM1：单点入口 gate 活跃→TaskKind::SinglePoint（AM1EvaluatorParams）分发，能量+逐原子电荷经 Frame gather 写回；两处 gate 不活跃→原路径不变) *(rework: 缝落在重构后的绑定层并扩展到 SinglePoint)*
- [ ] Write regression case regressions/launcher/lj_parity/ (atv_lj.py 固定种子分段采样、atv_am1.py 小分子单点输出能量+电荷、compare.py 双腿公差、reg_lj_parity.sh 四腿编排；在 ~/work/.x86_64 venv + buildenv-gcc/2026.03-mpich + LD_LIBRARY_PATH 镜像下运行；接入 regressions/CMakeLists.txt 的 regressions target，ATV_BUILD_LAUNCHER 门控) *(rework: 同案例新增 AM1 单点 parity 腿)*
- [ ] Verify parity and serial no-op: LJ——串行 `python atv_lj.py` vs `mpirun -np 2 atomiverse atv_lj.py`，逐采样点能量相对偏差 ≤ 1e-9、末位置 ≤ 1e-8 Å；AM1——能量相对偏差 ≤ 1e-10、逐原子电荷 ≤ 1e-10 e；整案 ≤ 60 s；并确认无标记时两个 run 入口都走原路径 *(rework: 增加 AM1 腿判据)*
- [ ] Run full check + test suite

## Testing strategy

- **单元（tests/，秒级，绝无 E2E）**：`LauncherGate`——无环境标记时 `active()==false`；**标记为他进程 PID ∧ MPI 未初始化 → inactive、无警告、不抛**（子进程形态，评审缺陷的直接回归测试）；**标记为本进程 PID ∧ MPI 未初始化 → inactive + 恰好一次警告、不抛**（宿主不一致的防御分支）；`LauncherArgs`——缺 argv、不存在路径、不可读文件抛出，可读文件通过。完整 launcher 流程明确**不属于**单元层。
- **边界/失败路径（regression 编排脚本内附带断言，2-rank）**：脚本抛异常 → traceback 打印、mpirun 退出码非零、worker 被 SHUTDOWN 释放、不挂起；缺脚本参数 / 不可读脚本 → 所有 rank 非零干净退出、不挂起；用户脚本内 spawn 的子 python `import atomiverse` 后照常串行工作（继承的标记因 PID 不匹配而静默失效）。
- **领域校验（parity regression，本链唯一 regression 案例，双腿）**：LJ 腿——同一 `atv_lj.py` 串行 python 托管 vs 2-rank launcher 托管，逐采样点总能量相对偏差 ≤ 1e-9、200 步末位置偏差 ≤ 1e-8 Å；AM1 腿——同一 `atv_am1.py`，能量相对偏差 ≤ 1e-10、逐原子电荷绝对偏差 ≤ 1e-10 e（确定性单点，无轨迹混沌）。整案 ≤ 60 s，在 `~/work/.x86_64` venv（atomiverse + molpy 可共同导入）+ `module load buildenv-gcc/2026.03-mpich` + `LD_LIBRARY_PATH` 镜像下经 `regressions` target 运行、绝不注册进 ctest。
- **不回归**：默认配置（`ATV_BUILD_LAUNCHER=OFF`）构建与既有全套单元测试不变；串行 `python atv.py` 语义逐字节不变（gate 一次 `getenv` 后走原路径）；`python-api-namespaces-am1` 重构后的 Python 包层不因本环再改一字。

## Out of scope

- **高层 Python 门面**（`atv.read` / `atv.DFT` / `atv.Optimize`）——命名空间重构（`python-api-namespaces-am1`）之外的门面在本链任何一环都不新增。
- **AM1 梯度与 AM1-MD**——本环只分发 AM1 **单点**（能量 + 电荷）；AM1 力提供者及其 MD 接线是后续 spec。
- **mpi4py 共存**——mpi4py 托管的 python 下 gate 因无标记保持串行路径；主动协作不做。
- **按对象粒度的同步动词**（per-object sync verbs）——本链只有任务粒度一次广播。
- **CUDA 后端分发**、**AM1 之外的其他单点任务类型**（PW/SCF/GPW 单点可作为 TaskKind 的后续扩展环）。
- **解释器发现/重定位与 wheel 打包**——第 05 环（wheel ships bin/atomiverse + interpreter discovery）。
- 本环不修改 01–03 环已落地的 wire 协议或 codec 语义（02 的 TaskKind/AM1EvaluatorParams 扩展属于 02 自己的并行修订，不由本环实现）。
