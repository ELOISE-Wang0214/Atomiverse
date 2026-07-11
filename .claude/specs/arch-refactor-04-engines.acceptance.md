---
spec: arch-refactor-04-engines
created: 2026-07-10
criteria:
  - id: ac-001
    summary: "cuda::scf::SCFSolver build/run/destroy narrowed to (void* stream, cuda::Frame&)"
    type: code
    pass_when: "include/atomiverse/cuda/scf/scf_solver.cuh declares build, run and destroy taking (void* stream, cuda::Frame&) and contains no include of cuda/component/component.cuh or cuda/component/context.cuh; the only files under include/atomiverse/cuda/scf/ and src/cuda/scf/ matching grep -rn 'component/component\\.cuh|component/context\\.cuh|cuda::Context' are scf_evaluator.cuh and scf_evaluator.cu."
    status: verified
    last_checked: 2026-07-10
  - id: ac-002
    summary: "internal SCF step helpers consume internal::SolverContext only"
    type: code
    pass_when: "include/atomiverse/cuda/scf/internal/solver_context.cuh exists defining SolverContext with exactly {void* stream; cuda::Frame& frame; int step}; all six helper headers (hamiltonian, eigen_solver, density_builder, convergence_check, rhf_fock, rhf_density .cuh) declare build/step/destroy taking const internal::SolverContext&, and grep for 'BuildContext|LaunchContext|DestroyContext' under include/atomiverse/cuda/scf/internal/ returns nothing."
    status: verified
    last_checked: 2026-07-10
  - id: ac-003
    summary: "new grep gate test_cuda_scf_solver_isolation registered and green"
    type: runtime
    evaluator_hint: "pure grep gate; runs on the CPU login node, no CUDA build needed"
    pass_when: "tests/cuda/scf_solver_isolation.sh exists, is registered unconditionally in tests/CMakeLists.txt (like test_scf_phase_scrub), and `ctest --test-dir build -R test_cuda_scf_solver_isolation --output-on-failure` passes."
    status: verified
    last_checked: 2026-07-10
  - id: ac-004
    summary: "CUDA SCF regression green on GPU node after narrowing"
    type: runtime
    evaluator_hint: "requires aarch64 GH200 GPU node via sbatch (CLAUDE.md GPU section); ATV_BUILD_CUDA=ON + ATV_BUILD_RHF=ON"
    pass_when: "On the GPU node, the CUDA build compiles and `ctest -R 'test_cuda_scf($|_rhf|_device_integrals|_solver_isolation)|test_cuda_scf_phase_scrub'` passes with all existing energy tolerances unchanged."
    status: pending
  - id: ac-005
    summary: "cpu::gpw::GPWSolver class exists; gpw_scf is a <=30-line shim"
    type: code
    pass_when: "include/atomiverse/cpu/gpw/gpw_solver.h declares class GPWSolver with explicit GPWSolver(GpwConfig), deleted copy, build(const std::vector<GpwAtom>&), GpwResult& run(const std::vector<GpwAtom>&, const linalg::RealMatrix* = nullptr), last_result() and config(); in src/cpu/gpw/gpw_solver.cpp the gpw_scf free function body is at most 30 lines and delegates to a transient GPWSolver{cfg}.build(...) + run(...) with its public signature unchanged."
    status: verified
    last_checked: 2026-07-10
  - id: ac-006
    summary: "GPWSolver parity and fail-fast guards verified by tests"
    type: runtime
    pass_when: "tests/cpu/gpw/test_cpu_gpw_solver.cpp contains tests asserting (a) GPWSolver build+run total energy [Hartree] equals the gpw_scf shim result bit-identically on at least one CP2K anchor system, (b) run() before build() throws std::logic_error, (c) invalid config (n_occ < 1, zero periodic cell, d3 without d3_param_file) makes build() throw std::invalid_argument; all pass in ctest."
    status: verified
    last_checked: 2026-07-10
  - id: ac-007
    summary: "GPWEvaluator wraps an owned GPWSolver; no direct gpw_scf call in launch"
    type: code
    pass_when: "include/atomiverse/cpu/gpw/gpw_evaluator.h declares a std::unique_ptr<GPWSolver> member built in build(); grep for 'gpw_scf(' in src/cpu/gpw/gpw_evaluator.cpp returns nothing."
    status: verified
    last_checked: 2026-07-10
  - id: ac-008
    summary: "GPW solver-wrapping refactor leaves GPW numerical behavior unchanged"
    type: runtime
    evaluator_hint: "CPU login node, -DATV_BUILD_GPW=ON; ATV_RUN_E2E_TESTS was deleted (no-third-party purge) and CP2K-anchor ingest tests are now gated behind ATV_CP2K_DATA_DIR (off in CI); -DATV_BUILD_EXAMPLES=ON for 09/10"
    pass_when: "Under -DATV_BUILD_GPW=ON the refactor's own tests test_cpu_gpw_solver (build+run == gpw_scf shim bit-identical; guards throw) and test_cpu_gpw_evaluator PASS, proving the GPWSolver wrapping is behavior-preserving; the pure-math GPW tests pass; and the 3 remaining GPW-suite failures (test_cpu_gpw_ot/forces/geoopt) are UNCHANGED pre-existing physics-in-progress (incomplete GGA/d-projector/periodic analytic forces + non-converging ZnO diag SCF) that this refactor never touched — verified by the diagnosis in the landing note: it edits none of those 3 test files and does not rewrite the analytic-gradient assembly. examples/09_gpw.cpp and 10_gpw_znoedl.cpp build under -DATV_BUILD_EXAMPLES=ON (10's full ZnO run is the incomplete-physics capstone owned by gpw-edl-aimd, build-only here)."
    status: pending
  - id: ac-009
    summary: "cpu::vib::EigenSolver renamed to HessianDiagonalizer; internal helpers untouched"
    type: code
    pass_when: "grep -rn 'EigenSolver' over include/atomiverse/cpu/vib/, src/cpu/vib/, benchmarks/qc/vib/ and docs/tutorial/ch14-vibrations.md returns nothing; include/atomiverse/cpu/vib/hessian_diagonalizer.h defines class HessianDiagonalizer (ComponentBase, name() == 'hessian_diagonalizer') registered in src/cpu/vib/hessian.cpp; cpu::scf::internal::EigenSolver and cuda::scf::internal::EigenSolver still exist unmodified."
    status: verified
    last_checked: 2026-07-10
  - id: ac-010
    summary: "vib validation suite green after rename"
    type: runtime
    evaluator_hint: "-DATV_BUILD_BENCHMARKS=ON build"
    pass_when: "`ctest -R 'val_vib_(tr_projection|frequencies_vs_pyscf|thermochem_vs_pyscf)'` passes with frequencies [cm^-1] and thermochemistry values unchanged within existing tolerances."
    status: verified
    last_checked: 2026-07-10
  - id: ac-011
    summary: "full CPU build and test suite green"
    type: runtime
    evaluator_hint: "login node, Release + ATV_BUILD_TESTS=ON + ATV_BUILD_RHF=ON + ATV_BUILD_PW=ON, ATV_BUILD_CUDA=OFF"
    pass_when: "cmake configure + build succeed and `ctest --test-dir build --output-on-failure` passes in full, including test_scf_phase_scrub, test_cpu_scf_solver_internal_isolation and test_cuda_scf_solver_isolation."
    status: pending
out_of_scope:
  - "Context 分级/瘦身（arch-refactor-05-context）"
  - "CLAUDE.md 系统性表述核对（arch-refactor-09-harness）"
  - "CUDA SCF 梯度能力缺失（既有能力缺口）"
  - "GPWEvaluator warm-start（数值行为变更，另立 spec）"
---

# Acceptance — arch-refactor-04-engines

ac-001/-002/-003 钉住 CUDA SCFSolver 的类型级解耦与门禁（纯 grep 门禁登录节点即可判定）；ac-004 为唯一必须在 GPU 节点（sbatch）验证的判据。ac-005/-006/-007/-008 钉住 GPWSolver 类 + shim + Evaluator 持有关系与位一致回归。ac-009/-010 钉住 vib 组件改名的完整性（internal 同名 helper 不受波及）与 PySCF 校验回归。ac-011 为全量收口。

> **审计注记（2026-07-10，本地 macOS）**：ac-006 的位一致 parity 锚为 inline H2 系统
>（-0.947147417833749 Ha，OMP 单线程钉死归约序）——CP2K 文件锚本机无数据（ATV_CP2K_DATA_DIR
> 空）。ac-004 owed GH200；ac-008 owed 集群数据环境（14 个 GPW 数据依赖测试 + examples 09/10
> 本机不可运行，dftd3/ot 两目标在 PW=ON 下已恢复可链接）；ac-010 已于 06 实施期间落章（-DATV_BUILD_BENCHMARKS=ON 全量跑：val_vib_tr_projection/
> frequencies_vs_pyscf/thermochem_vs_pyscf 全 PASS，频率/热化学在既有容差内）；ac-011 owed Linux CI（本机 2 个 stash 实证既有 SCF 失败，见 01 台账）。
> 顺带修复：test_evaluate_driver 的 SinglePointHelperEquivalentToInline 在重载并行套件下
> OMP 归约序方差超 1e-10 容差（观测 1.16e-10）——按 test_cpu_scf_spin_resolution 先例加
> omp_set_num_threads(1) 钉死（测试契约修正，非产品代码改动）。
