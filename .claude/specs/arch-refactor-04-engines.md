---
title: arch-refactor-04-engines — 方法引擎归位：CUDA SCFSolver 脱离 Context、补 GPWSolver、EigenSolver 消歧
status: approved
date: 2026-07-10
chain: arch-refactor (04/09)
depends_on: [arch-refactor-02-halo]
---

# arch-refactor-04-engines — 方法引擎归位

## Summary

本篇把两处偏离"方法 = inner Solver（纯数值 kernel）+ Evaluator（组件 wrapper）"模式的代码拉回该模式，并消除一处三义命名。第一，`cuda::scf::SCFSolver` 的 `build/run/destroy` 现在吞下完整 `cuda::Context&`（含 pair_list、host_box、exclusion CSR、force_rows、outer_status 等它不该看见的 workflow 状态），实际只消费 `runtime.stream`——签名窄化为 `void* stream`（沿用 context.cuh "void* 保持头文件 CUDA-free" 的既有约定），六个 `internal::*` step helper 改吃新的 `internal::SolverContext`，`cuda::scf::SCFEvaluator` 负责从 Context 取 stream 传入；原先只靠注释承诺的"不碰 outer_status"（scf_solver.cuh:14-16）改由类型系统保证，并新增 grep 隔离门禁。第二，GPW 方法只有自由函数 `gpw_scf`，破坏 SCFSolver/PWSolver/NDDOSolver 全域守约的 *Solver 类 + *Evaluator 模式——补 `cpu::gpw::GPWSolver` 类（build 一次性校验、run 可重入求解），`gpw_scf` 降为 ≤30 行 shim（对齐 `pw_scf()` 先例），`GPWEvaluator` 改持 GPWSolver。第三，注册进 Hessian driver 的组件 `cpu::vib::EigenSolver` 与 `cpu::scf::internal::EigenSolver`、`cuda::scf::internal::EigenSolver` 两个纯 helper 同名异役，组件改名 `HessianDiagonalizer`（对齐既有 `NDDODiagonalizer` 命名先例），internal helper 不动。全程 clean-break，无兼容层；所有既有数值结果不得移动。

## Domain basis

本篇为纯架构重构，不引入新物理；约束条件是数值不变性，锚点为既有领域校验：

- GPW 方法本体（高斯基 KS + 均匀网格 Hartree/XC/长程局域赝势）不变，实现按 Lippert, Hutter, Parrinello, *Mol. Phys.* **92**, 477 (1997), DOI: 10.1080/00268979709482119 的 GPW 方案，由 CP2K 位一致锚点测试（`tests/cpu/gpw/test_cpu_gpw_{h2,ar,nacl,zno}_cp2k.cpp`）钉住；`gpw_scf` 函数体逐字移入 `GPWSolver::run`，能量 [Hartree] 必须位一致。
- 振动分析（TR 投影 + 质量加权 Hessian 对角化）为 PySCF `hessian/thermo.py::harmonic_analysis` 的逐字移植，由 `benchmarks/qc/vib/validation/`（vs PySCF 频率 [cm^-1]、热化学）钉住；改名不触算法。
- CUDA RHF 能量 [Hartree] 由 `tests/cuda/scf/test_rhf.cu` 既有 CPU/CUDA parity 容差钉住；stream 窄化不改任何 kernel 与调用次序。

## Design

### 1. cuda::scf::SCFSolver 解耦 cuda::Context

- **新增 `include/atomiverse/cuda/scf/internal/solver_context.cuh`**：定义 `ATV::cuda::scf::internal::SolverContext { void* stream; cuda::Frame& frame; int step; }`（step = SCF 迭代序号；build/destroy 传 0）。仅 include `atomiverse/cuda/frame/frame.cuh`（frame 支柱，架构允许），CUDA-free。
- **`SCFSolver` 签名窄化**：`build(void* stream, cuda::Frame&)` / `run(void* stream, cuda::Frame&)` / `destroy(void* stream, cuda::Frame&)`。头文件移除 `component/component.cuh` include（改 include frame.cuh + solver_context.cuh）；实现内 `BuildContext{runtime, frame}` / `LaunchContext{runtime, frame, 0, iter}` 改为 `internal::SolverContext{stream, frame, iter}`。
- **钦定 `void* stream`**（三选一决议）：与 `cuda::Context::stream` 的既有 void* 约定一致，.cu 内 `static_cast<cudaStream_t>` 模式不变。已否决：`cudaStream_t` 裸参（头文件被迫拉 cuda_runtime.h，破坏 CUDA-free 约定）；新造 `cuda::StreamHandle` 包装类型（为单一字段引入新抽象，若 05-context 分级后需要类型化句柄可再统一引入）。
- **六个 internal helper**（hamiltonian / eigen_solver / density_builder / convergence_check / rhf_fock / rhf_density）的 `build/step/destroy` 参数从 `const BuildContext&`/`const LaunchContext&`/`const DestroyContext&` 换为 `const internal::SolverContext&`；实现中 `ctx.runtime.stream` → `ctx.stream`，`ctx.frame` 用法不变（rhf_fock 读原子坐标、convergence_check/rhf_density 写 `frame.scalars->energy`，均为合法 frame 依赖）。
- **`SCFEvaluator`（组件 wrapper，职责不变）**：`solver_->build(ctx.runtime.stream, ctx.frame)` 等三处调用点改写；scf_evaluator.cuh 本身继续合法 include component.cuh（它就是组件）。
- **新增 grep 门禁 `tests/cuda/scf_solver_isolation.sh`**（CTest `test_cuda_scf_solver_isolation`，仿 `test_scf_phase_scrub` 无条件注册、纯 grep 不需 CUDA 构建）：`include/atomiverse/cuda/scf/` 与 `src/cuda/scf/` 下除 `scf_evaluator.cuh` / `scf_evaluator.cu` 外，不得出现 `component/component.cuh`、`component/context.cuh` include token，亦不得出现 `cuda::Context`、`BuildContext`、`LaunchContext`、`DestroyContext` token。
- scf_solver.cuh:14-16 及 scf_solver.cu 头注释中"never touches outer_status"的措辞更新为"由签名类型保证"。

### 2. 补 cpu::gpw::GPWSolver 类

- **公共面（对齐 PWSolver 形态，move-only）**：
  ```cpp
  class GPWSolver {
    public:
      explicit GPWSolver(GpwConfig config);
      // 禁拷贝、可移动（同 PWSolver）
      void build(const std::vector<GpwAtom>& reference);
      GpwResult& run(const std::vector<GpwAtom>& atoms,
                     const linalg::RealMatrix* guess_C = nullptr);
      const GpwResult& last_result() const;
      const GpwConfig& config() const;
  };
  ```
- **build() 一次性**：fail-fast 校验（atoms/basis 非空、`n_occ >= 1`、周期路径 cell 三边 > 0、`d3 == true` 时 `d3_param_file` 非空——违例抛 `std::invalid_argument`，no-fallback 规则），记录 `built_`。**run() 可重入**：现 `gpw_scf` 函数体逐字移入（网格、lattice images、S/T/Vpp、SCF 循环均按现状每 run 重建——与 PWSolver "几何相关结构每 run 重建"注释一致），结果写入 owned `result_` 并返回引用；`built_ == false` 时抛 `std::logic_error`（fail-fast，不静默自动 build）。
- **`gpw_scf` 降为 ≤30 行 shim**（签名不变，26 处调用点零改动）：构造瞬态 `GPWSolver{cfg}` → `build(atoms)` → `run(atoms, guess_C)` → 按值返回。对齐 `scf_solve()` / `pw_scf()` 先例。
- **`GPWEvaluator` 改持 solver**：成员 `std::unique_ptr<GPWSolver> solver_;`。`build()` 从 build-time frame 组装 `GpwConfig`（cell [Bohr]、n_occ = Σz_val/2、`ATV_GPW_OT` 环境变量检查移到 build——每 launch 重读环境变量本就无意义）并构造 + `solver_->build(...)`；`launch()` 按当前几何组装 atoms 后 `solver_->run(atoms)`（冷启动保持现状，不传 guess_C，数值零变化），能量/解析力写回逻辑不变。新增 guard：launch 时 frame box 与 build-time cell 不一致则抛 `std::runtime_error`（GPW 不支持变胞；原实现静默用新 box 属未定义行为，无在-tree 用例依赖）。
- **已否决备选（记录在案）**：豁免方案——保留自由函数 `gpw_scf` 为 GPW 的唯一面、改 CLAUDE.md 文档承认例外。否决理由：破坏 SCFSolver/PWSolver/NDDOSolver/AM1Solver（含未提交 semiempirical 新代码）全域一致的 *Solver 类 + *Evaluator 模式，且 CLAUDE.md:184 已按 `GPWSolver` 类记载该形态——文档是对的，代码欠账。

### 3. EigenSolver 三义消歧

- `cpu::vib::EigenSolver`（ComponentBase 组件，注册进 `cpu::vib::Hessian` driver 的 `diagonalize` phase）改名 **`cpu::vib::HessianDiagonalizer`**——名称直陈职责（对角化质量加权 Hessian），对齐既有 `NDDODiagonalizer` 先例；`name()` 返回值 `"vib_eigen_solver"` → `"hessian_diagonalizer"`。文件改名 `eigen_solver.{h,cpp}` → `hessian_diagonalizer.{h,cpp}`，根 CMakeLists.txt:202 源列表同步。
- `cpu::scf::internal::EigenSolver`、`cuda::scf::internal::EigenSolver` **不动**：处于 internal 命名空间、职责名准确（通用广义/对称本征求解），且受既有隔离门禁（`test_cpu_scf_solver_internal_isolation`、`tests/cuda/scf/test_scf.cu` 静态断言）保护。
- 调用点与文档引用同步：`src/cpu/vib/hessian.cpp:41`、`benchmarks/qc/vib/validation/test_tr_projection.cpp:40`，以及 vib 头文件 doc 注释中的交叉引用（mass_weighting.h / thermochem.h / hessian_workspace.h / hessian_result.h / normal_mode_writer.h / hessian.h）与 `docs/tutorial/ch14-vibrations.md`。

## Files to create or modify

CUDA SCF 解耦：
- include/atomiverse/cuda/scf/internal/solver_context.cuh (new)
- include/atomiverse/cuda/scf/scf_solver.cuh
- include/atomiverse/cuda/scf/internal/hamiltonian.cuh
- include/atomiverse/cuda/scf/internal/eigen_solver.cuh
- include/atomiverse/cuda/scf/internal/density_builder.cuh
- include/atomiverse/cuda/scf/internal/convergence_check.cuh
- include/atomiverse/cuda/scf/internal/rhf_fock.cuh
- include/atomiverse/cuda/scf/internal/rhf_density.cuh
- src/cuda/scf/scf_solver.cu
- src/cuda/scf/hamiltonian.cu
- src/cuda/scf/eigen_solver.cu
- src/cuda/scf/density_builder.cu
- src/cuda/scf/convergence_check.cu
- src/cuda/scf/rhf_fock.cu
- src/cuda/scf/rhf_density.cu
- src/cuda/scf/scf_evaluator.cu
- tests/cuda/scf_solver_isolation.sh (new)
- tests/cuda/scf/test_scf.cu
- tests/CMakeLists.txt

GPWSolver：
- include/atomiverse/cpu/gpw/gpw_solver.h
- src/cpu/gpw/gpw_solver.cpp
- include/atomiverse/cpu/gpw/gpw_evaluator.h
- src/cpu/gpw/gpw_evaluator.cpp
- tests/cpu/gpw/test_cpu_gpw_solver.cpp
- tests/cpu/gpw/test_cpu_gpw_evaluator.cpp

vib 改名：
- include/atomiverse/cpu/vib/hessian_diagonalizer.h (renamed from include/atomiverse/cpu/vib/eigen_solver.h)
- src/cpu/vib/hessian_diagonalizer.cpp (renamed from src/cpu/vib/eigen_solver.cpp)
- src/cpu/vib/hessian.cpp
- CMakeLists.txt
- benchmarks/qc/vib/validation/test_tr_projection.cpp
- include/atomiverse/cpu/vib/mass_weighting.h
- include/atomiverse/cpu/vib/thermochem.h
- include/atomiverse/cpu/vib/hessian_workspace.h
- include/atomiverse/cpu/vib/hessian_result.h
- include/atomiverse/cpu/vib/normal_mode_writer.h
- include/atomiverse/cpu/vib/hessian.h
- docs/tutorial/ch14-vibrations.md

## Tasks

- [ ] Write failing tests for cuda::scf stream isolation：新增 tests/cuda/scf_solver_isolation.sh（grep 门禁：include/atomiverse/cuda/scf/ 与 src/cuda/scf/ 除 scf_evaluator.{cuh,cu} 外无 component/component.cuh、component/context.cuh include token 及 cuda::Context/BuildContext/LaunchContext/DestroyContext token），在 tests/CMakeLists.txt 无条件注册为 test_cuda_scf_solver_isolation（仿 test_scf_phase_scrub），并在 tests/cuda/scf/test_scf.cu 增加 static_assert 断言 SCFSolver::build/run/destroy 可以 (void*, cuda::Frame&) 调用——当前均 RED
- [ ] Implement cuda::scf::SCFSolver stream narrowing：新增 include/atomiverse/cuda/scf/internal/solver_context.cuh（SolverContext{void* stream; cuda::Frame& frame; int step}），窄化 include/atomiverse/cuda/scf/scf_solver.cuh 与 src/cuda/scf/scf_solver.cu 的 build/run/destroy 为 (void* stream, cuda::Frame&) 并移除 component.cuh include，六个 internal 头（hamiltonian/eigen_solver/density_builder/convergence_check/rhf_fock/rhf_density.cuh）及对应 src/cuda/scf/*.cu 改吃 const internal::SolverContext&（ctx.runtime.stream → ctx.stream），src/cuda/scf/scf_evaluator.cu 三处调用点改传 ctx.runtime.stream，同步更新 doxygen（@param stream 语义 + outer_status 承诺改为类型保证）
- [ ] Write failing tests for cpu::gpw::GPWSolver：tests/cpu/gpw/test_cpu_gpw_solver.cpp 增加 GPWSolver 类路径测试（build+run 能量与 gpw_scf shim 位一致；run-before-build 抛 std::logic_error；n_occ/cell/d3 非法配置 build 抛 std::invalid_argument），tests/cpu/gpw/test_cpu_gpw_evaluator.cpp 增加 evaluator-持-solver 行为测试（能量不变 + launch 时 box 变化抛 std::runtime_error）——当前均 RED
- [ ] Implement cpu::gpw::GPWSolver in include/atomiverse/cpu/gpw/gpw_solver.h + src/cpu/gpw/gpw_solver.cpp：类面按 Design §2（ctor/build/run/last_result/config，move-only），gpw_scf 函数体逐字移入 run()，gpw_scf 降为 ≤30 行 shim（签名不变），doxygen 含单位（能量 [Hartree]、坐标 [Bohr]）
- [ ] Implement GPWEvaluator owning GPWSolver in include/atomiverse/cpu/gpw/gpw_evaluator.h + src/cpu/gpw/gpw_evaluator.cpp：build() 组装 GpwConfig（cell/n_occ/ATV_GPW_OT）并构造 + build 持有的 unique_ptr<GPWSolver>，launch() 组装当前几何 atoms 后 solver_->run(atoms)（冷启动、数值零变化）+ 变胞 guard，launch 内不再出现 gpw_scf 调用
- [ ] Implement vib rename EigenSolver → HessianDiagonalizer：include/atomiverse/cpu/vib/eigen_solver.h → hessian_diagonalizer.h、src/cpu/vib/eigen_solver.cpp → hessian_diagonalizer.cpp（类名/name() 同步），更新根 CMakeLists.txt 源列表、src/cpu/vib/hessian.cpp 注册点、benchmarks/qc/vib/validation/test_tr_projection.cpp，以及 mass_weighting.h/thermochem.h/hessian_workspace.h/hessian_result.h/normal_mode_writer.h/hessian.h 的 doc 注释交叉引用与 docs/tutorial/ch14-vibrations.md
- [ ] Verify CUDA SCF regression on GPU node：sbatch（CLAUDE.md GPU 节点流程，GPU/buildenv-gcccuda/2026.03-cu13.0）构建 ATV_BUILD_CUDA=ON + ATV_BUILD_RHF=ON，test_cuda_scf / test_cuda_scf_rhf / test_cuda_scf_device_integrals / test_cuda_scf_solver_isolation / test_cuda_scf_phase_scrub 全绿
- [ ] Verify CPU regressions：GPW 全套（ctest -R test_cpu_gpw，含 -DATV_RUN_E2E_TESTS=ON 的 CP2K 锚点）数值不变，examples/09_gpw.cpp 与 examples/10_gpw_znoedl.cpp 构建运行冒烟，vib 校验（-DATV_BUILD_BENCHMARKS=ON 下 val_vib_tr_projection / val_vib_frequencies_vs_pyscf / val_vib_thermochem_vs_pyscf）全绿
- [ ] Run full check + test suite：登录节点 CPU 配置（Release + RHF + PW，CUDA OFF）cmake 构建零警告新增，ctest --output-on-failure 全绿（含 test_scf_phase_scrub、test_cpu_scf_solver_internal_isolation、新 test_cuda_scf_solver_isolation 纯 grep 门禁）

## Testing strategy

- **Happy path**：CUDA SCF 三条 run 路径（synthetic / RHF host / RHF device integrals）经既有 tests/cuda/scf/ 全绿（GPU 节点）；GPWSolver 类路径与 gpw_scf shim 路径能量位一致；GPW 全套 13 个测试 + examples 09/10 数值不变；vib 频率/热化学 vs PySCF 校验不变。
- **Edge cases**：GPWSolver run-before-build 抛 std::logic_error；build 对非法 n_occ / cell / d3 配置抛 std::invalid_argument（no-fallback）；GPWEvaluator 变胞 guard 抛 std::runtime_error；CUDA SCFSolver 奇电子数构造 guard（既有）不受签名窄化影响。
- **隔离门禁（新增 + 既有回归）**：test_cuda_scf_solver_isolation（新，纯 grep，登录节点即可跑）；test_scf_phase_scrub、test_cpu_scf_solver_internal_isolation、tests/cuda/scf/test_scf.cu 的 not-ComponentBase 静态断言继续全绿——SCF/PW/GPW internal helper 永不继承 ComponentBase、永不注册 Plan 的链级不变量不动。
- **领域校验**（science.required）：GPW CP2K 锚点（H2/Ar/NaCl/ZnO）能量位一致；CUDA RHF 与 CPU 参考能量在既有容差内；vib vs PySCF 频率 [cm^-1] / 热化学在既有容差内。
- **构建门禁注意**：ATV_BUILD_CUDA 仅 GPU 节点可编译，CUDA 判据需 sbatch 验证（CLAUDE.md GPU 节点流程）；纯 grep 门禁与 CPU 侧全部判据在登录节点验证。

## Out of scope

- **Context 分级/瘦身**（cuda::Context 字段按消费者拆分）→ arch-refactor-05-context（本篇使 SCFSolver 先行脱离 Context，是 05 的前置受益项，但不动 Context 本体）。
- **CLAUDE.md / .claude/architecture.md 系统性表述核对**（GPWSolver 记载在本篇落地后自动为真；全量文档一致性修订）→ arch-refactor-09-harness。
- **CUDA SCF 梯度能力缺失**（cuda::scf::SCFEvaluator 无 Force capability）——既有能力缺口，与本篇解耦无关。
- **GPWEvaluator warm-start**（AIMD 传 guess_C）——改变数值行为，如需另立 spec。
- **cpu/cuda pw、semiempirical（NDDO/AM1 未提交代码）、torch backend**——不触碰。
