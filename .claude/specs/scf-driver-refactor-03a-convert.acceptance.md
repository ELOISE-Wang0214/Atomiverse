---
slug: scf-driver-refactor-03a-convert
criteria:
  - id: ac-001
    summary: RHF converged energy matches the SCFSolver reference to 1e-12 Ha
    type: code
    evaluator_hint: "target: test_cpu_scf_method_driver_parity (RHF)"
    pass_when: |
      For an H2/H2O-scale STO-3G molecule, cpu::scf::RHF driven by CPUDriver<scf::Phase>
      converges to last_result().energy within 1e-12 Ha of the in-tree SCFSolver reference
      on identical integrals; the test runs in the fast unit suite (seconds), constructs no
      out-of-library consumer, and reads no literature value.
    status: pending
  - id: ac-002
    summary: UHF converged energy matches the SCFSolver reference to 1e-12 Ha
    type: code
    evaluator_hint: "target: test_cpu_scf_method_driver_parity (UHF)"
    pass_when: |
      For a small open-shell case, cpu::scf::UHF driven by CPUDriver<scf::Phase> matches the
      in-tree SCFSolver reference energy within 1e-12 Ha (and spin_contamination to 1e-10)
      on identical integrals, in the fast unit suite.
    status: pending
  - id: ac-003
    summary: RKS converged energy matches the SCFSolver reference to 1e-12 Ha
    type: code
    evaluator_hint: "target: test_cpu_scf_method_driver_parity (RKS)"
    pass_when: |
      For a small closed-shell case, cpu::scf::RKS (PBE/MatchPySCF grid) driven by
      CPUDriver<scf::Phase> matches the in-tree SCFSolver reference energy within 1e-12 Ha
      on identical integrals and grid, in the fast unit suite.
    status: pending
  - id: ac-004
    summary: UKS converged energy matches the SCFSolver reference to 1e-12 Ha
    type: code
    evaluator_hint: "target: test_cpu_scf_method_driver_parity (UKS)"
    pass_when: |
      For a small open-shell case, cpu::scf::UKS driven by CPUDriver<scf::Phase> matches the
      in-tree SCFSolver reference energy within 1e-12 Ha on identical integrals and grid,
      in the fast unit suite.
    status: pending
  - id: ac-005
    summary: MeanFieldMethod is a CPUDriver<scf::Phase> wired config-free from members
    type: code
    evaluator_hint: "inspect: mean_field_method.h/.cpp ctor + drive_scf_"
    pass_when: |
      MeanFieldMethod<D> derives CPUDriver<scf::Phase>; its constructor wires the five 02
      stages on scf::Phase fock/mix/diagonalize/density/check via on(phase).add with
      parameters taken directly from members, and registers ScfConverged{tol,
      conv_tol_grad>0?conv_tol_grad:sqrt(tol)} via add_criterion (Criterion component,
      driver-pulled verdict); the drive_scf_ path contains no SCFConfig,
      SCFSolver, or ScfLoop token; drive_scf_ is a private member function (not a free
      function); no new make_* free function is introduced anywhere in the change.
    status: pending
  - id: ac-006
    summary: Staged run/launch coexistence is explicit, bounded, and routes one loop
    type: code
    evaluator_hint: "inspect: mean_field_method.h build/launch + drive_scf_"
    pass_when: |
      MeanFieldMethod<D> still derives HF/Component and retains build()/launch(); both the
      run() overloads and the Component launch() funnel into the same private drive_scf_
      (same CPUDriver<scf::Phase> loop), a molecule run via run(frame) and via launch (bound
      into Evaluate) yields identical energy, and an in-code comment marks run/launch
      coexistence as debt resolved in 04.
    status: pending
  - id: ac-007
    summary: mp2()/fci() stay members with unchanged results
    type: code
    evaluator_hint: "target: test_cpu_scf_mp2 / test_cpu_scf_fci"
    pass_when: |
      cpu::scf::RHF::mp2(frame) and ::fci(frame) remain members that, after
      drive_scf_(compute_eri=true), consume integrals_sp_ + last_result() and return results
      unchanged (existing test_cpu_scf_mp2 / test_cpu_scf_fci pass without value edits).
    status: pending
  - id: ac-008
    summary: Structure admits ROHF via one fock component + one preset file
    type: docs
    evaluator_hint: "design-review: final cpu/scf structure"
    pass_when: |
      Review of the final structure confirms adding ROHF would require only one new
      fock-phase component plus one new cpu/scf/rohf.{h,cpp} preset file (ctor defaults +
      resolve_occupation + name()), with zero edits to the CPUDriver<scf::Phase> loop, the
      other four methods, or the mix/diagonalize/density/check stages; each per-method file
      (rhf/uhf/rks/uks) holds only ctor defaults + resolve_occupation + name().
    status: pending
  - id: ac-009
    summary: Per-method file split, umbrella header, sources registered and documented
    type: docs
    evaluator_hint: "inspect: cpu/scf files + root CMakeLists"
    pass_when: |
      MeanFieldMethod<D> lives in cpu/scf/mean_field_method.{h,cpp}; RHF/UHF/RKS/UKS live in
      cpu/scf/{rhf,uhf,rks,uks}.{h,cpp}; methods.h is an umbrella including them (SCF=RKS
      alias kept); src/cpu/scf/methods.cpp is removed and the new .cpp are listed in the root
      CMakeLists ATV_BUILD_RHF block; every new public symbol carries a Doxygen @brief.
    status: pending
  - id: ac-010
    summary: Whole suite green in seconds; SCFSolver reference and its tests still pass
    type: runtime
    evaluator_hint: "target: full ctest + existing tests/cpu/scf/*"
    pass_when: |
      `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      && cmake --build build -j && ctest --test-dir build --output-on-failure` passes the
      whole suite in seconds; test_cpu_scf_method_driver_parity passes and every existing
      tests/cpu/scf/* test (including test_scf_solver against the retained reference) still
      passes. No regressions/ case is added in 03a (the literature-value e2e ships in 03b).
    status: pending
---

# Acceptance criteria — scf-driver-refactor-03a-convert

把 `MeanFieldMethod<D>` 改成 `CPUDriver<scf::Phase>` 驱动 02 stage 链、逐位对拍在树 `SCFSolver` 参考，参考实现留到 03b 才删。全部 `code`/`docs`/`runtime`，逐条可判。

- **ac-001..004** —— 四方法逐位一致（≤1e-12 Ha），unit tier、秒级、对在树参考直接对拍。
- **ac-005** —— 派生 `CPUDriver<scf::Phase>`、成员→stage-ctor 无 config 接线、`drive_scf_` 为私有成员、无新 `make_*`。
- **ac-006** —— run/launch 共存债务显式有界、双入口共用同一循环、代码内指向 04。
- **ac-007** —— `mp2`/`fci` 仍为成员且结果不变。
- **ac-008** —— ROHF 可扩展性设计评审 + per-method 文件只含构造器默认/resolve_occupation/name。
- **ac-009 / ac-010** —— 文件拆分/伞头/源登记/文档；整套快套件秒级全绿且旧参考及其测试不变。

无 `scientific`/`performance`/`ui_runtime`：逐位对拍由 ctest 单测在库内完成（归 `code`）；out-of-library 文献值 e2e（`scientific`/`runtime`）属 03b；无性能预算；无前端。
