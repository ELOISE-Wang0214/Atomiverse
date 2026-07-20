---
slug: scf-driver-refactor-04a-nest
criteria:
  - id: ac-001
    summary: cpu::Driver erased base promoted from SubDriver ops; CPUDriver implements it
    type: code
    evaluator_hint: "inspect: cpu/driver/driver.h + cpu_driver.h"
    pass_when: |
      cpu/driver/driver.h defines a non-template polymorphic cpu::Driver with
      virtual build_on(Backend&), run_nested_erased(Context&, Frame&, int), destroy_on(Backend&),
      default_budget(), and a virtual destructor; CPUDriver<P> derives and implements it, with
      build_on = compile()+build_plan, run_nested_erased = run_nested, destroy_on = destroy_plan
      (the sub_driver.h:94-115 ops).
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: PendingBinding + Plan node is a Component*/cpu::Driver variant threaded through build()
    type: code
    evaluator_hint: "inspect: facade/driver.h + plan.h"
    pass_when: |
      detail::PendingBinding and the compiled Plan/PlanBuilder node are a variant of
      Component* and DriverBinding{cpu::Driver*, budget}; Driver::build() (driver.h:219-225)
      compiles both node kinds into plan_; the driver-node alternative is device-parameterized
      so backends without nesting still compile.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: on(phase).add/add_ref driver overloads with ownership + budget default; bare add deleted
    type: code
    evaluator_hint: "inspect: facade/driver.h PhaseBuilder + owned_drivers_"
    pass_when: |
      on(phase).add(CPUDriver<InnerPhase>&&, order, budget) moves the child into owned_drivers_
      (vector<unique_ptr<cpu::Driver>>); on(phase).add_ref(CPUDriver<InnerPhase>&, order, budget)
      borrows; budget<0 uses child.default_budget(); the bare add(driver) is deleted with a
      use-on(phase) message; the driver overloads are SFINAE-disabled on backends without nesting.
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: CPU backend dispatches node kind; component probes skip driver-nodes
    type: code
    evaluator_hint: "inspect: backend.cpp + cpu_driver.h for_each_pending"
    pass_when: |
      cpu::Backend build_components/step/destroy_components dispatch driver-nodes to
      build_on/run_nested_erased/destroy_on and component-nodes to build/launch/destroy; the
      component-oriented walks (for_each_pending cpu_driver.h:276, bind_halo_participants_ :303,
      constraint scan :191, validate_long_range_contract_ :335, validate_components) skip
      driver-nodes (they carry no Role).
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: Native nesting runs the inner budget each outer iteration
    type: code
    evaluator_hint: "target: test_cpu_driver_nesting"
    pass_when: |
      In test_driver_nesting, an inner Minimize bound via outer.on(md::Phase::end_of_step).add_ref(inner,50)
      runs up to its budget on each outer iteration, behaving identically to the deleted SubDriver
      path, built once and looped per outer step (no per-iteration allocation).
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: Inner convergence cannot stop the outer loop (isolation preserved)
    type: code
    evaluator_hint: "target: test_cpu_driver_nesting (isolation)"
    pass_when: |
      With an inner loop that converges early, the outer loop still runs its full step budget;
      the assertion is preserved from the old test_sub_driver, backed by the unchanged run_nested
      ctx.loop re-point (cpu_driver.h:432-437) so the child verdict lands only on the child status_.
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: Borrowed inner driver is inspectable after the parent run
    type: code
    evaluator_hint: "target: test_cpu_driver_nesting (borrow)"
    pass_when: |
      An inner driver bound via add_ref exposes converged()/stop_reason() reflecting its last
      nested run after outer.run() returns (the non-owning lvalue-registration contract).
    status: verified
    last_checked: 2026-07-20
  - id: ac-008
    summary: validate_binding_kinds_ unchanged for components; driver-nodes bypass it correctly
    type: code
    evaluator_hint: "inspect: driver.h:270-284 + on() driver path"
    pass_when: |
      validate_binding_kinds_ (driver.h:270-284) still throws on a component whose default_bindings
      target another workflow (test_driver_binding_kind passes); a driver-node has no default_bindings
      and is bound only via on(parent-phase), so its slot PhaseKind is trivially the parent's and the
      child's internal phase enum never reaches the parent's kind gate.
    status: verified
    last_checked: 2026-07-20
  - id: ac-009
    summary: sub_driver.h deleted; no SubDriver token; SCF untouched; CUDA additive
    type: code
    evaluator_hint: "grep: SubDriver"
    pass_when: |
      cpu/driver/sub_driver.h and the old test_sub_driver.cpp are removed; an evaluator grep finds
      no SubDriver token under include/, src/, tests/ (excluding git history); MeanFieldMethod is
      still a Component via HF (no SCF change this link); the facade variant change is additive and
      the driver overloads SFINAE-guarded so CUDA facade compiles unchanged.
    status: verified
    last_checked: 2026-07-20
  - id: ac-010
    summary: Whole suite green in seconds
    type: runtime
    evaluator_hint: "target: full ctest"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the whole
      suite in seconds; test_cpu_driver_nesting passes and no existing driver/MD/minimizer/evaluate
      test regresses. No regressions/ case is added in 04a.
    status: verified
    last_checked: 2026-07-20
---

# Acceptance criteria — scf-driver-refactor-04a-nest

把 `SubDriver` 的适配职责提升进 Plan/后端的原生节点表示，删除公开包装类，子 verdict 隔离与误绑门不变。全部 `code`/`runtime`。

- **ac-001..004** —— `cpu::Driver` 擦除基、`PendingBinding`/`Plan` 变体、`on().add` driver 重载（拥有/借用/预算默认）、后端节点分派 + 组件探针跳过 driver 节点。
- **ac-005..007** —— 原生嵌套跑内层预算、内层收敛不能停外层、借用式 run 后可查。
- **ac-008** —— 误绑门对组件不变、driver 节点只经父 phase 绑定绕过之。
- **ac-009 / ac-010** —— `sub_driver.h` 删净、SCF 不动、CUDA 加法式；整套快套件秒级全绿。

无 `scientific`/`performance`/`ui_runtime`：纯基础设施、无端到端功能（本链回归在 03b）、无性能预算、无前端。
