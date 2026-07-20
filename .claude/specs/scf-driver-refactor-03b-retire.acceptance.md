---
slug: scf-driver-refactor-03b-retire
criteria:
  - id: ac-001
    summary: H2O/STO-3G RHF regression reproduces the cited literature energy
    type: runtime
    evaluator_hint: "target: regressions/scf/h2o_rhf_sto3g (ATV_BUILD_REGRESSIONS)"
    pass_when: |
      Built with -DATV_BUILD_REGRESSIONS=ON and run via the `regressions` target,
      regressions/scf/h2o_rhf_sto3g/reg_h2o_rhf_sto3g.cpp (an out-of-library consumer:
      only public cpu/scf/methods.h + frame.h, linking atomiverse) drives cpu::scf::RHF
      with STO-3G on H2O at the REFERENCE.md geometry (R(OH)=0.95 A, HOH=104.5 deg) and
      reproduces the cited E_RHF/STO-3G = -74.961754 Ha (LibreTexts Live Textbook of
      Physical Chemistry, Peverati, sec. 28.2; cross-checked vs PySCF 2.13; STO-3G origin
      Hehre/Stewart/Pople JCP 51, 2657 (1969)) within 1e-3 Ha in <=60 s; not part of the
      default ctest run.
    amended: "2026-07-20 — spec's original -74.9420 Ha/Szabo-Ostlund pairing was a
      misattribution (that number is Crawford Programming Project 3 water at stretched
      R=1.1 A/104.0 deg; verified via PySCF 2.13 reproducing all three numbers). Amended
      per the spec's own instruction to verify the cited digits at implementation time."
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: Gradient path decoupled from SCFConfig with identical forces
    type: code
    evaluator_hint: "target: test_gradient* / test_gradient_vs_libcint"
    pass_when: |
      The analytic gradient path is an internal class (internal::ScfGradient, constructor
      taking the discrete method parameters, held privately by the method; the FD oracle is
      a member of it or lives in its consuming test) and integral_derivatives takes discrete
      parameters — no SCFConfig type and no free helper function in the gradient path; the
      existing gradient tests pass with forces bit-identical to pre-refactor.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: SCFResultHandle loses config; SCFResult gains spin; cube writer migrated
    type: code
    evaluator_hint: "inspect: mean_field_method.h + scf_result.h + density_cube_writer.cpp"
    pass_when: |
      SCFResultHandle has no config field and config_sp_ is gone; SCFResult carries a
      SpinType spin set in the method epilogue; density_cube_writer reads result.spin
      (no handle_.config); test_density_cube, test_orbital_info_observer, and example
      16_density_cube compile and pass; orbital_cube_writer/orbital_info_observer are
      unchanged (they never read config).
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: SCFSolver/ScfLoop/SCFConfig/make_config fully removed
    type: code
    evaluator_hint: "grep: SCFSolver|ScfLoop|SCFConfig|make_config"
    pass_when: |
      An evaluator grep over include/, src/, examples/, tests/ finds zero live references
      to SCFSolver, ScfLoop, SCFConfig, or make_config, excluding: git history, the root
      cpu/scf/{scf_solver,scf_types}.h tombstone #error strings, and CUDA's own
      cuda::scf::SCFSolver plus incidental cross-reference comments. Expressed as an
      evaluator-run grep, not a committed script.
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: internal/scf_types.h keeps only SpinResolution; HamiltonianBuilder discrete-only; tombstones updated
    type: code
    evaluator_hint: "inspect: internal/scf_types.h + hamiltonian_builder.h + root tombstones"
    pass_when: |
      internal/scf_types.h defines SpinResolution (+ its includes) and no SCFConfig;
      HamiltonianBuilder exposes only the 02 discrete-parameter prepare/step (no SCFConfig
      overload); cpu/scf/scf_solver.h and cpu/scf/scf_types.h remain #error tombstones with
      text updated to point at RHF/UHF/RKS/UKS (methods.h) and SpinResolution respectively.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: Doc examples and delta_scf_reference migrated to the member API
    type: code
    evaluator_hint: "target: test_delta_scf/test_koopmans/test_level_shift/test_accelerators"
    pass_when: |
      vib/hessian.h and evaluate.h @code examples and tests/cpu/scf/delta_scf_reference.h
      use the RHF/UHF member API with no SCFConfig token; test_delta_scf, test_koopmans,
      test_level_shift, and test_accelerators pass with unchanged expected values.
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: Solver-internal tests retired; CUDA isolation bash gate deleted; coverage preserved
    type: code
    evaluator_hint: "inspect: tests/cpu/scf/stages + tests/CMakeLists"
    pass_when: |
      test_scf_solver_internals and test_scf_solver_internal_isolation are removed with
      their ctest registrations; the six meaningful helper assertions (DensityBuilder
      factor-2 + empty occupation, EigenSolver generalized eig, DualConvergence iter-0-false
      + both-criteria + reset) are present in the 02 stage test; tests/cuda/scf_solver_isolation.sh
      and its test_cuda_scf_solver_isolation registration are deleted; no new .sh gate is added
      and the stale rdf_molrs_migration.sh comment is swept.
    status: verified
    last_checked: 2026-07-20
  - id: ac-008
    summary: Surviving isolation concerns expressed as evaluator greps, not scripts
    type: code
    evaluator_hint: "grep: internal isolation properties"
    pass_when: |
      Evaluator greps confirm: (i) the 02 stage helpers (EigenSolver etc.) do not inherit
      cpu::Component; (ii) facade-public SCF headers (methods.h, mean_field_method.h) do not
      include cpu/scf/internal/*; (iii) the CUDA SCF solver takes its narrow SolverContext,
      not the driver Context. No committed .sh enforces any of these.
    status: verified
    last_checked: 2026-07-20
  - id: ac-009
    summary: Build wiring updated and the regression documented
    type: docs
    evaluator_hint: "inspect: root CMakeLists + tests/CMakeLists + regressions/CMakeLists + REFERENCE.md"
    pass_when: |
      The root CMakeLists no longer lists src/cpu/scf/scf_solver.cpp; tests/CMakeLists drops
      the three retired registrations; regressions/CMakeLists registers h2o_rhf_sto3g behind
      ATV_BUILD_REGRESSIONS on the regressions target; REFERENCE.md records the cited value,
      exact geometry, citation, and tolerance; new symbols carry Doxygen @brief.
    status: verified
    last_checked: 2026-07-20
    verified_by: human
  - id: ac-010
    summary: Whole suite green in seconds on the driver-ized methods
    type: runtime
    evaluator_hint: "target: full ctest"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the
      whole suite in seconds with no SCFSolver reference remaining; every tests/cpu/scf/*
      test runs on the driver-ized RHF/UHF/RKS/UKS.
    status: verified
    last_checked: 2026-07-20
    note: "Non-MPI suite 123/123 green; all 25 SCF tests green; reg_h2o_rhf_sto3g PASS.
      Pre-existing MPI halo/coul failures (missing Potential in HaloFieldParticipant
      single-rank fixture) are orthogonal to this SCF retirement."
---

# Acceptance criteria — scf-driver-refactor-03b-retire

删除已被 03a 证冗余的旧 SCF 机器、把最后的 `SCFConfig` 消费方迁到离散/成员 API、交付首个高斯 SCF 文献值回归。全部 `code`/`docs`/`runtime`。

- **ac-001** —— 文献值端到端回归（唯一 `runtime` 回归关口，`ATV_BUILD_REGRESSIONS`，库外消费方，≤60s）。
- **ac-002 / ac-003** —— 梯度路径与 handle 去 `SCFConfig`，力值/消费方行为不变。
- **ac-004 / ac-005** —— `SCFSolver`/`ScfLoop`/`SCFConfig`/`make_config` grep-count-zero（排除墓碑串与 CUDA 自有类）、`scf_types.h` 只留 `SpinResolution`、墓碑文案更新。
- **ac-006 / ac-007 / ac-008** —— 文档/参考测试迁移、求解器内部测试退役 + CUDA 隔离 bash 门删除（隔离关注改 evaluator grep，不设脚本）。
- **ac-009 / ac-010** —— 构建接线 + 回归文档；整套快套件秒级全绿。

无 `scientific`/`performance`/`ui_runtime`：按 evaluator 协议，库内 `regressions/` e2e 归 `runtime`（高精度 bench 才是 `scientific`，不在本仓）；无性能预算；无前端。
