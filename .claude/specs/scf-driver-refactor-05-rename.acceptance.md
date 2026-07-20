---
slug: scf-driver-refactor-05-rename
criteria:
  - id: ac-001
    summary: CMake option + define + GPW dependency renamed
    type: code
    evaluator_hint: "inspect: CMakeLists.txt option/define/GPW checks"
    pass_when: |
      The CMake option is ATV_BUILD_HF (was ATV_BUILD_RHF), the compile define is ATV_HAS_HF
      (was ATV_HAS_RHF), and the GPW dependency check requires ATV_BUILD_HF; tests/CMakeLists.txt,
      examples/CMakeLists.txt, and tests/architecture/*.cmake use the new names.
    status: pending
  - id: ac-002
    summary: All source #ifdef/#define renamed on both backends, logic unchanged
    type: code
    evaluator_hint: "grep: ATV_HAS_RHF across include/ src/"
    pass_when: |
      Every ATV_HAS_RHF #ifdef/#define in CPU (cpu/observer/density.h) and CUDA
      (cuda/scf/*, cuda/observer/mulliken.cuh, cuda/recorder/density_recorder.cuh, src/cuda/scf/*)
      is renamed ATV_HAS_HF with no logic change.
    status: pending
  - id: ac-003
    summary: Scripts and packaging renamed
    type: code
    evaluator_hint: "grep: -DATV_BUILD_RHF in tools/slurm/python"
    pass_when: |
      tools/ci_cpu_unit.sh, tools/gpu_job.sh, slurm/{job_cpu_full,job_cpu_gpw_test,job_cpu_xc_smoke}.sh,
      and python/pyproject.toml use -DATV_BUILD_HF / ATV_BUILD_HF.
    status: pending
  - id: ac-004
    summary: Docs, examples, and tests renamed
    type: docs
    evaluator_hint: "grep: ATV_BUILD_RHF in docs/examples/tests"
    pass_when: |
      CLAUDE.md/AGENTS.md/README.md/.claude/architecture.md option tables + build recipes, docs/*,
      examples README + example/test build comments and #ifdefs all use ATV_BUILD_HF/ATV_HAS_HF.
    status: pending
  - id: ac-005
    summary: Zero old-token residue outside git history and this chain's specs
    type: code
    evaluator_hint: "grep: ATV_BUILD_RHF|ATV_HAS_RHF"
    pass_when: |
      An evaluator grep for ATV_BUILD_RHF or ATV_HAS_RHF over the tree (excluding git history and
      .claude/specs/*) returns zero matches.
    status: pending
  - id: ac-006
    summary: Full suite green under the renamed flag
    type: runtime
    evaluator_hint: "target: full ctest under ATV_BUILD_HF"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the whole
      suite in seconds, identical to the pre-rename run.
    status: pending
---

# Acceptance criteria — scf-driver-refactor-05-rename

纯机械改名，行为零变。`code`/`docs`/`runtime`。ac-001..004 分组改名；ac-005 旧 token 计数清零（排除 git 历史 + 本链 spec）；ac-006 新开关下全绿。无 `scientific`/`performance`/`ui_runtime`。
