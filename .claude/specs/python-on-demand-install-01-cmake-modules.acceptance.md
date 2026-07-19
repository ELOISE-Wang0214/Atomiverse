---
spec: python-on-demand-install-01-cmake-modules
created: 2026-07-19
criteria:
  - id: ac-001
    summary: Product-module options exist with classical defaults ON
    type: code
    pass_when: |
      Root configure defines ATV_BUILD_MD, MC, MINIMIZER, VIB, PATH,
      SEMIEMPIRICAL, POTENTIALS, NEIGHBOR, EVALUATE as options defaulting ON;
      existing ATV_BUILD_RHF/PW/GPW/KSPACE/CUDA defaults unchanged.
    status: verified
    last_checked: 2026-07-19
  - id: ac-002
    summary: Auto-enable DAG forces required dependencies
    type: runtime
    pass_when: |
      Architecture cmake test: MD=ON forces POTENTIALS+NEIGHBOR ON; PATH=ON
      forces MINIMIZER ON; KSPACE=ON forces POTENTIALS ON; no cycles.
    status: verified
    last_checked: 2026-07-19
  - id: ac-003
    summary: Invalid GPW/TORCH/LAUNCHER combos still FATAL_ERROR
    type: runtime
    pass_when: |
      Configure GPW without RHF or PW FATALs; TORCH without CUDA and LAUNCHER
      without MPI/PYTHON still FATAL.
    status: verified
    last_checked: 2026-07-19
  - id: ac-004
    summary: OFF modules omit domain sources from atomiverse
    type: runtime
    pass_when: |
      MD=OFF excludes src/cpu/md and src/cuda/md from atomiverse SOURCES;
      SEMIEMPIRICAL=OFF excludes semiempirical sources.
    status: verified
    last_checked: 2026-07-19
  - id: ac-005
    summary: Public ATV_HAS_* matches each enabled product module
    type: code
    pass_when: |
      Each enabled product gate emits PUBLIC ATV_HAS_*; disabled gates do not.
    status: verified
    last_checked: 2026-07-19
  - id: ac-006
    summary: nanobind only compiles bind_*.cpp for ON modules
    type: code
    pass_when: |
      bind_md/minimizer/pw/semiempirical sources are conditional; bind.cpp
      guards init with #ifdef ATV_HAS_*.
    status: verified
    last_checked: 2026-07-19
  - id: ac-007
    summary: Unit and regression targets register only for ON modules
    type: code
    pass_when: |
      tests/ and regressions/ CMake wrap domain targets in if(ATV_BUILD_*);
      default-all-ON still registers pre-spec classical tests.
    status: verified
    last_checked: 2026-07-19
  - id: ac-008
    summary: Regression module_compile_defs passes on full-classical defaults
    type: runtime
    pass_when: |
      regressions/architecture/module_compile_defs exits 0 within 60s under
      default classical ON; fails compile if required ATV_HAS_* missing.
    status: verified
    last_checked: 2026-07-19
  - id: ac-009
    summary: Default CI-shaped configure still builds and tests
    type: runtime
    pass_when: |
      Configure+build equivalent to tools/ci_cpu_unit.sh succeeds; no new
      link failures from gated objects.
    status: verified
    last_checked: 2026-07-19
  - id: ac-010
    summary: Docs on new public cmake surface
    type: docs
    pass_when: |
      cmake/AtvProductModules.cmake documents each product option, DAG edges,
      and always-on core in file-header comments.
    status: verified
    last_checked: 2026-07-19
out_of_scope:
  - Python CLI and pure wheel (02/03)
  - Physics algorithm changes
---

# Acceptance — python-on-demand-install-01-cmake-modules

完成态：产品模块选项与 DAG 正确、OFF 不编译、HAS 宏与 nanobind/测试注册一致，默认 monorepo 构建不回退。
