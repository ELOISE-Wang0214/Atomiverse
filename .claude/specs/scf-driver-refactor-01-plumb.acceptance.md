---
slug: scf-driver-refactor-01-plumb
criteria:
  - id: ac-001
    summary: scf::Phase enum defined with fock/diagonalize/density/mix/check + COUNT
    type: code
    evaluator_hint: "inspect: include/atomiverse/facade/stage.h namespace ATV::scf"
    pass_when: |
      include/atomiverse/facade/stage.h defines namespace ATV::scf with
      enum class Phase whose ordered members are fock,mix,diagonalize,density,check
      (plus a COUNT sentinel; Fock-space DIIS precedes diagonalization per
      scf_solver.cpp:241-291); slot_of(scf::Phase::fock)==0 through
      slot_of(scf::Phase::check)==4 and slot_count<ATV::scf::Phase>()==5,
      with no hand-written static_cast<int> over the enum.
    status: pending
  - id: ac-002
    summary: PhaseKind::SCF added and kind_of/make_binding route scf::Phase to it
    type: code
    evaluator_hint: "inspect: core/binding.h PhaseKind + stage.h kind_of"
    pass_when: |
      PhaseKind (core/binding.h) has an SCF value, phase_kind_name(PhaseKind::SCF)
      returns "SCF", kind_of(scf::Phase::fock)==PhaseKind::SCF, and
      make_binding(scf::Phase::fock, 0).kind==PhaseKind::SCF.
    status: pending
  - id: ac-003
    summary: facade/loop.h and cpu/component/context.h are untouched (bags withdrawn)
    type: code
    evaluator_hint: "diff: facade/loop.h + cpu/component/context.h"
    pass_when: |
      Per the user rulings, this link adds NO fields to LoopMetrics (facade/loop.h is
      byte-for-byte unchanged) and NO ScfServices/ctx.scf slice or ScfWorkspace forward
      declaration (cpu/component/context.h is byte-for-byte unchanged). SCF state will
      flow via method-member matrices wired into stage constructors by reference (02/03a).
    status: pending
  - id: ac-005
    summary: CPUDriver<scf::Phase> instantiates via a PhaseTraits<scf::Phase> spec
    type: code
    evaluator_hint: "inspect: cpu/driver/cpu_driver.h PhaseTraits<ATV::scf::Phase>"
    pass_when: |
      cpu_driver.h contains a PhaseTraits<ATV::scf::Phase> specialization with
      eval_phase==fock, neighbor_phase==fock, end_phase==check, and
      CPUDriver<ATV::scf::Phase> instantiates in test_scf_phase_plumbing.
    status: pending
  - id: ac-006
    summary: scf_phase_scrub script and its ctest registration removed, no dangling refs
    type: code
    evaluator_hint: "grep: scf_phase_scrub across repo"
    pass_when: |
      tests/cuda/scf_phase_scrub.sh no longer exists; tests/CMakeLists.txt has no
      test_scf_phase_scrub add_test block; an evaluator grep for 'scf_phase_scrub'
      over include/, src/, tests/, examples/, tools/ returns zero live references
      (the removed stale comments in scf_solver_isolation.sh and pw_scf.cuh are gone).
      No new .sh gate is added.
    status: pending
  - id: ac-007
    summary: New public symbols documented; stage.h NOTE rewritten
    type: docs
    evaluator_hint: "inspect: doxygen @brief on new symbols"
    pass_when: |
      scf::Phase, kind_of(scf::Phase), ScfServices, PhaseTraits<scf::Phase>, and the
      two new LoopMetrics fields each carry a Doxygen @brief (units stated on delta_e/
      diis_error); the former stage.h:132-137 "deliberately no scf::Phase" NOTE is
      replaced by a block describing scf::Phase as the CPU Gaussian-SCF driver phase enum.
    status: pending
  - id: ac-008
    summary: Full unit-test suite green in seconds after the plumbing lands
    type: runtime
    evaluator_hint: "target: test_scf_phase_plumbing + full ctest"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` completes
      with the whole suite passing in seconds; test_scf_phase_plumbing passes and no
      pre-existing test regresses.
    status: pending
  - id: ac-009
    summary: Criterion component kind + driver-pulled verdict polling landed additively
    type: code
    evaluator_hint: "inspect: cpu/component/criterion.h + cpu_driver.h add_criterion"
    pass_when: |
      cpu/component/criterion.h defines Criterion deriving cpu::Component with a pure
      virtual Verdict verdict(const LoopState&); CPUDriver exposes add_criterion (owning)
      and add_criterion_ref (borrowing); after each iteration the driver polls registered
      criteria, composes verdicts by severity Stop > Converged > Continue, and reduces
      once over MPI together with the existing handler poll; components still cannot
      write loop status (ctx.loop borrow semantics unchanged); existing LoopHandler and
      add_handler are untouched; test_scf_phase_plumbing exercises a toy Criterion.
    status: pending
---

# Acceptance criteria — scf-driver-refactor-01-plumb

纯增量门面/核心/基础设施铺设 + 删除硬阻塞 scrub 门。全部为 `code`/`docs`/`runtime` 型，可由第三方逐条 yes/no 判定，无需重读规格。

- **ac-001 / ac-002 / ac-005** —— 阶段类型底座：枚举形状、`PhaseKind::SCF` 判别、`PhaseTraits` 特化，保证 03 能实例化 `CPUDriver<scf::Phase>`。
- **ac-003** —— 两个"袋"已按用户裁定撤销：`facade/loop.h` 与 `context.h` 逐字不动（状态走成员 + 构造器引用）。
- **ac-006** —— 硬阻塞 grep 门删除且无悬挂引用；以 evaluator 直接 grep 表达，不引用任何脚本作为强制手段，也不新增 `.sh` 门。
- **ac-007** —— 文档与代码注释同步（反转 stage.h 里"故意没有 scf::Phase"的旧结论）。
- **ac-008** —— 整套单测秒级全绿的 runtime 关口（本环无 `regressions/` 案例——SCF 驱动首次端到端跑能量在 03，该链级回归案例及其 runtime 验收由 03 交付）。
- **ac-009** —— `Criterion` 组件类别（终止判据即组件，决策⑨）+ driver 拉取式 verdict 轮询，纯增量落地；`LoopHandler` 全树迁移属 06。

无 `scientific`/`performance`/`ui_runtime` 型：本环无物理数值、无性能预算、无前端。
