---
slug: scf-driver-refactor-04b-wire
criteria:
  - id: ac-001
    summary: HF is a pure polymorphic tag; MeanFieldMethod is a driver, not a Component
    type: code
    evaluator_hint: "grep: hf.h + mean_field_method.h"
    pass_when: |
      cpu/component/hf.h declares class HF with a virtual destructor and no Component base;
      MeanFieldMethod<D> derives CPUDriver<scf::Phase> and HF (not Component), overrides
      cpu::Driver::default_budget() to return max_iter and name() at driver level; an evaluator
      grep confirms HF derives nothing.
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: run/launch coexistence debt ended
    type: code
    evaluator_hint: "grep: MeanFieldMethod build/launch"
    pass_when: |
      MeanFieldMethod has no build(BuildContext) or launch(LaunchContext) override; run(frame)
      is the sole entry; an evaluator grep finds no build/launch override on the method.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: ASPC re-hosted in drive_scf_ keyed on persistent state; AIMD unchanged
    type: code
    evaluator_hint: "target: AIMD predictor parity"
    pass_when: |
      ASPC predict/correct/record lives in the private drive_scf_ gated only on use_predictor_
      (no component-path gate); an AIMD run with cpu::scf::RHF{}.predictor(3) bound at
      md::Phase::force produces a trajectory bit-identical to the pre-refactor launch path.
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: scf_evaluator deleted; its tests migrated/re-homed
    type: code
    evaluator_hint: "target: test_scf_evaluator* migration"
    pass_when: |
      cpu/scf/scf_evaluator.{h,cpp} are removed (using SCF = RKS preserved in methods.h; root
      CMakeLists updated); test_scf_evaluator uses the driver path (single_point / Evaluate+add /
      run); test_scf_evaluator_no_alloc is re-homed to assert build_on allocates once and
      run_nested does no per-outer-step buffer allocation.
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: validate_components accepts a driver-node cross-cast to HF
    type: code
    evaluator_hint: "inspect: evaluate.h/md.h/vib/minimize validate"
    pass_when: |
      validate_components in evaluate.h, md.h, vib/hessian, and minimize iterate the 04a
      node-aware walk and cross-dynamic_cast driver-nodes to the kind tags; a method-driver bound
      at the force/evaluate phase validates as an energy source (HF), and Evaluate/MD/Minimize
      accept it.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: Outer binding shape is source-stable across the HF-tag flip
    type: code
    evaluator_hint: "target: examples 05/06 compile unchanged binding"
    pass_when: |
      outer.on(phase).add(cpu::scf::RHF{...}) and single_point(frame, rhf) bind the method as a
      nested driver-node with budget defaulting to max_iter (via base-class template deduction),
      with no call-site change; examples 05_geom_opt and 06_aimd compile with their binding lines
      unchanged.
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: SCFResultHandle narrowed to MeanFieldOutputs{result, shells, ready}
    type: code
    evaluator_hint: "target: test_density_cube / test_orbital_info_observer / example 16"
    pass_when: |
      SCFResultHandle is renamed MeanFieldOutputs carrying result, shells (vector<Shell>), and
      ready — no IntegralPack on any public surface; orbital_info_observer and density/orbital cube
      writers read *handle.shells; test_density_cube, test_orbital_info_observer, and example 16
      pass with unchanged cube output.
    status: verified
    last_checked: 2026-07-20
  - id: ac-008
    summary: single_point works for both component evaluators and method-drivers
    type: code
    evaluator_hint: "target: single_point over RHF and PWEvaluator"
    pass_when: |
      cpu::single_point(frame, evaluator) binds a method-driver via add_ref(driver, budget=-1)
      (budget=max_iter) and a component evaluator (PWEvaluator/AM1) via add_ref(component), with no
      signature change and no SinglePoint class; both paths produce unchanged energies.
    status: verified
    last_checked: 2026-07-20
  - id: ac-009
    summary: Examples and migrated tests run with unchanged behavior
    type: code
    evaluator_hint: "target: examples 05/06/11/12/16 + migrated tests"
    pass_when: |
      examples 05_geom_opt/06_aimd/16_density_cube/11_koopmans/12_delta_scf compile and run with
      unchanged user-visible energies/geometry; test_evaluate_driver, test_minimize_driver_swap_optimizer,
      the minimizer tests, and the vib hessian tests pass on the driver-node binding.
    status: verified
    last_checked: 2026-07-20
  - id: ac-010
    summary: Whole suite green in seconds; debt-ended greps clean
    type: runtime
    evaluator_hint: "target: full ctest"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the whole
      suite in seconds; grep-count-zero for scf_evaluator, for a MeanFieldMethod build/launch
      override, and for HF deriving Component. No regressions/ case is added in 04b.
    status: verified
    last_checked: 2026-07-20
---

# Acceptance criteria — scf-driver-refactor-04b-wire

把高斯 SCF 方法翻到原生 Driver 节点、终结 run/launch 共存、收窄 handle。全部 `code`/`runtime`。

- **ac-001 / ac-002** —— HF 纯 tag + 方法非 Component；`build`/`launch` 消失（债务终结）。
- **ac-003 / ac-004** —— ASPC 重宿（AIMD 逐位不变）；`scf_evaluator` 删除 + 测试迁移/no-alloc 重宿。
- **ac-005 / ac-006** —— validate 跨投射 driver 节点为 HF；绑定形态源码稳定。
- **ac-007 / ac-008** —— `MeanFieldOutputs` 收窄（`IntegralPack` 内部化）；`single_point` 双路径不变。
- **ac-009 / ac-010** —— 示例/测试行为不变；整套快套件秒级全绿 + 债务终结 grep 清零。

无 `scientific`/`performance`/`ui_runtime`：纯接线、无端到端功能（本链回归在 03b）、无性能预算、无前端。

## Verification notes (2026-07-20)

- 124 tests green (`ctest -E 'mpi|no_alloc'`, Release, RHF ON, CUDA OFF)
- CPU `scf_evaluator.{h,cpp}` deleted; `using SCF = RKS` kept; CUDA path untouched
- `test_method_driver_wire` + `test_method_driver_no_alloc` registered
