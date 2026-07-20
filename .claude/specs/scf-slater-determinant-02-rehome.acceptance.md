---
slug: scf-slater-determinant-02-rehome
criteria:
  - id: ac-001
    summary: Three enums live verbatim in new public scf_enums.h
    type: code
    evaluator_hint: "inspect: include/atomiverse/cpu/scf/scf_enums.h"
    pass_when: |
      MeanFieldType{HF,KS}, SpinType{Restricted,Unrestricted}, InitGuess{Hcore,MINAO} are
      defined (enumerators unchanged) in include/atomiverse/cpu/scf/scf_enums.h and no longer
      in scf_result.h; the scf_types.h and gradient.h #error tombstones are untouched.
    status: verified
    last_checked: 2026-07-20
  - id: ac-002
    summary: SCFInitialGuess rehomed into scf_guess.h
    type: code
    evaluator_hint: "inspect: scf_guess.h"
    pass_when: |
      struct SCFInitialGuess is defined verbatim in include/atomiverse/cpu/scf/scf_guess.h
      (with the xtensor include added) and not in scf_result.h.
    status: verified
    last_checked: 2026-07-20
  - id: ac-003
    summary: GradientResult public visibility preserved in gradient_result.h
    type: code
    evaluator_hint: "inspect: gradient_result.h + internal/gradient.h + mean_field_method.h"
    pass_when: |
      struct GradientResult is defined verbatim in include/atomiverse/cpu/scf/gradient_result.h
      (a public, non-internal header); mean_field_method.h (forces() return, last_grad_) and
      internal/gradient.h compile against it; it is no longer in scf_result.h.
    status: verified
    last_checked: 2026-07-20
  - id: ac-004
    summary: IterationRecord deleted as redundant; zero survivors
    type: code
    evaluator_hint: "grep: IterationRecord|\\.trace\\b in cpu scf"
    pass_when: |
      grep for IterationRecord over include/ src/ tests/ examples/ docs/ returns zero matches;
      SCFResult has no trace member; the push_back formerly at
      src/cpu/scf/mean_field_method.cpp:450-451 is removed; no other SCFResult member changed.
    status: verified
    last_checked: 2026-07-20
  - id: ac-005
    summary: scf_result.h end state is SCFResult-only; exemptions untouched
    type: code
    evaluator_hint: "inspect: scf_result.h; grep: toy SCFResult"
    pass_when: |
      include/atomiverse/cpu/scf/scf_result.h defines exactly one type (SCFResult), with
      <vector> dropped; scf::toy::SCFResult in solve.h and the NDDO semiempirical types are
      byte-identical to before.
    status: verified
    last_checked: 2026-07-20
  - id: ac-006
    summary: Include surface rerouted; no includer needs scf_result.h for moved types
    type: code
    evaluator_hint: "inspect: the 8 direct includers"
    pass_when: |
      internal/hamiltonian_builder.h:16 includes scf_enums.h instead of scf_result.h;
      mean_field_method.h, internal/gradient.h, internal/integral_derivatives.h, and
      src/cpu/recorder/density_cube_writer.cpp include the new homes for the enum/guess/gradient
      types they use while retaining scf_result.h only where SCFResult itself is used;
      internal/mp2_pack.h, internal/fci_pack.h, and orbital_cube_writer.cpp are unchanged.
    status: verified
    last_checked: 2026-07-20
  - id: ac-007
    summary: New-home include-surface unit test passes in ctest
    type: code
    evaluator_hint: "target: test_scf_type_homes"
    pass_when: |
      tests/cpu/scf/test_scf_type_homes.cpp is registered in tests/CMakeLists.txt, includes only
      scf_enums.h / scf_guess.h / gradient_result.h (not scf_result.h), exercises all three enums
      and constructs SCFInitialGuess and GradientResult, and passes.
    status: verified
    last_checked: 2026-07-20
  - id: ac-008
    summary: Zero behavior change — full build and test suite green
    type: runtime
    pass_when: |
      cmake --build build succeeds and ctest --test-dir build --output-on-failure passes with no
      assertion edits in any pre-existing test (test_method_driver_parity, delta-SCF tests, SCF
      stages, cube-writer tests all green as-is).
    status: verified
    last_checked: 2026-07-20
---

# Acceptance criteria

ac-001/002/003 固定三处新家与逐字搬迁；ac-004 固定 `IterationRecord` 冗余删除（含 `SCFResult::trace` 及唯一写点）；ac-005 固定 `scf_result.h` 末态与两处点名豁免；ac-006 固定八处直接 include 面的重路由（含三处明确不动项）；ac-007 是 03-unbag 依赖的可独立编译性质的绑定门；ac-008 是行为零变化的运行时总门。
