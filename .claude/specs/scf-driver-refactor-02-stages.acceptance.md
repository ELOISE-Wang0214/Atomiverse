---
slug: scf-driver-refactor-02-stages
criteria:
  - id: ac-001
    summary: State flows via ctor-wired named references — no Workspace type, no ctx slice
    type: code
    evaluator_hint: "grep: ScfWorkspace + inspect stage ctors"
    pass_when: |
      No ScfWorkspace type and no ctx.scf slice exist anywhere in the tree; every stage
      component constructor receives the named state it reads as const references
      (e.g. const RealMatrix& P_alpha) and the state it writes as mutable references
      (e.g. RealMatrix& F_alpha), plus discrete parameters; the isolated unit tests wire
      test-local matrices into those constructors. No DIIS history is shared state
      (it is a MixStage member).
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: FockStage reproduces HamiltonianBuilder bit-identically, config-free
    type: code
    evaluator_hint: "target: test_cpu_scf_stages (FockStage)"
    pass_when: |
      FockStage inherits cpu::Component, binds scf::Phase::fock via make_binding, and
      for a fixed workspace P produces F_alpha/F_beta/J/e_xc_grid/e_exact_exchange/
      int_n bit-identical to internal::HamiltonianBuilder on the same input; its
      constructor takes discrete parameters only (no SCFConfig token anywhere in the
      FockStage or its HamiltonianBuilder discrete entry); build() allocates, launch()
      allocates no fixed buffer.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: MixStage reproduces DIIS error + extrapolation + level shift; history internal
    type: code
    evaluator_hint: "target: test_cpu_scf_stages (MixStage)"
    pass_when: |
      MixStage binds scf::Phase::mix and, for a fixed workspace (F,P,S,X,iter),
      writes err_alpha/err_beta and the DIIS-extrapolated + level-shifted F_alpha/F_beta
      bit-identical to internal::diis_error + DIISMixer::step + internal::apply_level_shift
      across restricted (single-channel) and unrestricted (dual-channel) branches; the
      DIIS history is a member of MixStage (not shared state) and is cleared by a
      reset() the method driver calls at run start.
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: DiagonalizeStage privately holds the closed EigenSolver kernel
    type: code
    evaluator_hint: "target: test_cpu_scf_stages (DiagonalizeStage)"
    pass_when: |
      internal::EigenSolver does NOT inherit cpu::Component; DiagonalizeStage binds
      scf::Phase::diagonalize, holds EigenSolver as a private member, and for a fixed
      workspace F produces C/eig bit-identical to EigenSolver::step (n_dropped==0) or
      step_orthogonal (n_dropped>0), reusing the alpha result when restricted.
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: DensityStage reproduces density build + damping + energy for all branches
    type: code
    evaluator_hint: "target: test_cpu_scf_stages (DensityStage)"
    pass_when: |
      DensityStage binds scf::Phase::density and, for a fixed workspace C and old P,
      writes P_alpha/P_beta (factor 2 restricted), applies damping only when DIIS is
      inactive (fraction-of-old convention), and writes e_total/e_one/e_coulomb
      bit-identical to internal::DensityBuilder + the scf_solver.cpp energy formulas
      across all four branches (restricted/unrestricted × HF/KS), using old P for HF
      and new P for KS.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: ScfConverged is one Criterion component — launch measures, verdict judges
    type: code
    evaluator_hint: "target: test_cpu_scf_stages (ScfConverged)"
    pass_when: |
      ScfConverged (cpu/scf/stages/scf_converged.h) derives the Criterion kind
      (cpu/component/criterion.h), binds scf::Phase::check, and in launch() computes
      delta_energy_ (=+inf on iter 0) and diis_error_ (combined alpha/beta Frobenius)
      into its own plainly-named members, matching DualConvergenceCheck's dE/err_norm —
      it never touches LoopMetrics; its verdict() returns Verdict::Converged iff
      |delta_energy_|<tol_energy && diis_error_<tol_grad, else Verdict::Continue,
      judging its own members only (never the iteration counter); the component never
      writes loop status (verdict is driver-pulled).
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: New stages are config-free Components wired by named references via ctors
    type: code
    evaluator_hint: "grep: SCFConfig under cpu/scf/stages"
    pass_when: |
      Every constructor under include/atomiverse/cpu/scf/stages/ (scf_converged.h included)
      takes discrete parameters and named state references with zero SCFConfig /
      config-struct / Workspace-bag references (internal included); each stage inherits
      cpu::Component, declares its phase via make_binding, guards work with
      ctx.at(scf::Phase::…), and accesses state only through its ctor-wired references.
    status: verified
    last_checked: 2026-07-20
  - id: ac-008
    summary: New stage sources registered and documented
    type: docs
    evaluator_hint: "inspect: root CMakeLists.txt ATV_BUILD_RHF list + doxygen"
    pass_when: |
      All new src/cpu/scf/stages/*.cpp are listed in the root CMakeLists.txt ATV_BUILD_RHF
      source block; every new public symbol (the five stages incl. ScfConverged and
      their ctor reference parameters) carries a Doxygen @brief with units where applicable.
    status: pending
  - id: ac-009
    summary: Whole suite green in seconds; old SCF path bit-unchanged
    type: runtime
    evaluator_hint: "target: full ctest + existing tests/cpu/scf/*"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the
      whole suite in seconds; test_cpu_scf_stages passes and every pre-existing
      tests/cpu/scf/* test (test_scf_solver, test_uhf_solver, test_uks_solver, …) still
      passes with unchanged energies (old ScfLoop/SCFSolver remains authoritative).
    status: verified
    last_checked: 2026-07-20
---

# Acceptance criteria — scf-driver-refactor-02-stages

把五个内部 helper 加壳成绑定 `scf::Phase` 的 stage 组件、状态经构造器具名引用直连（无 Workspace、无切片）、收敛判据成为一个 Criterion 组件（launch 测量存自身真名成员 + verdict 判定），全程纯增量、逐位对拍。全部 `code`/`docs`/`runtime` 型，可逐条 yes/no 判定。

- **ac-001** —— 状态走构造器具名引用；全树无 `ScfWorkspace`/`ctx.scf`；DIIS 历史为 MixStage 私有。
- **ac-002 / ac-003 / ac-004 / ac-005** —— 四个计算 stage 各自对拍其包装 helper 逐位一致，并落实 config-free 构造、封闭 `EigenSolver`、四分支能量、Fock-DIIS 分支。
- **ac-006** —— `ScfConverged` 为单个 Criterion 组件：`launch()` 测量存自身真名成员（不碰 LoopMetrics）、`verdict()` 只判自身成员，裁决由 driver 拉取；置于 flag-gated 的 `cpu/scf/stages/`。
- **ac-007** —— "无 config、无 Workspace 袋"与组件契约（`ctx.at` + 构造器引用）的结构审查。
- **ac-008 / ac-009** —— 源登记 + 文档；整套单测秒级全绿且旧 SCF 路径能量逐位不变。

无 `scientific`/`performance`/`ui_runtime` 型：本环无外部 bench 文献验证（逐位对拍由 ctest 单测在库内完成，归 `code`）、无性能预算、无前端。
