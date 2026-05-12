---
slug: driver-infra-auto-injection
schema: v1
criteria:
  - id: ac-001
    summary: install_infra() hook fires before validate_components()
    type: code
    pass_when: "test_install_infra_hook test binary passes"
    test_path: "tests/cpu/test_install_infra_hook.cpp"

  - id: ac-002
    summary: set_* API produces same Plan as raw add()/on().add() — NVE, NPT, and custom-hook coexistence
    type: code
    pass_when: "MDPlanEquivalence.NVELJSetApiMatchesRawAddPath, MDPlanEquivalence.NPTSetApiMatchesRawAddPathWithExplicitVirialStress, MDPlanEquivalence.CustomHookViaOnAddCoexistsWithSetters all pass"
    test_path: "tests/cpu/md/test_md_plan_equivalence.cpp"

  - id: ac-003
    summary: Singular setters — idempotent on same value, throw on conflict
    type: code
    pass_when: "MDSetApi.{SetPairTwiceWithSameTypeIsIdempotent, SetNeighborTwiceThrowsWhenNameDiffers, SetIntegratorTwiceWithSameNameIdempotent, SetIntegratorRejectsDifferentTypeAfterFirst} pass; MinimizeSetApi.{SetOptimizerSingularSameNameIdempotent, SetOptimizerRejectsDifferentTypeAfterFirst} pass"
    test_path: "tests/cpu/md/test_md_set_api.cpp; tests/cpu/test_minimize_evaluate_set_api.cpp"

  - id: ac-004
    summary: Capability-dep auto-injection — NPT/Barostat reads()==Stress, LJ needs_neighbor_list(), with driver-level dedup by canonical-provider name
    type: code
    pass_when: "MDSetApi.{NPTReadsStress, BarostatReadsStress, LJNeedsNeighborList, NPTAutoInjectsVirialStress, ExplicitVirialStressSuppressesAutoInject, NoStressReadsLeavesVirialStressOut} all pass"
    test_path: "tests/cpu/md/test_md_set_api.cpp"

  - id: ac-005
    summary: Missing set_neighbor when needs_neighbor_list()==true throws at build() across MD / Minimize / Evaluate
    type: code
    pass_when: "MDSetApi.PairWithoutNeighborListThrowsAtBuild, MinimizeSetApi.PairWithoutNeighborListThrowsAtBuild, EvaluateSetApi.PairWithoutNeighborListThrowsAtBuild all pass; happy paths PairWithNeighborListBuildsClean (×3) also pass; EvaluateSetApi.EnergyOnlyEvaluatorBuildsWithoutPair confirms non-pair eval path stays open"
    test_path: "tests/cpu/md/test_md_set_api.cpp; tests/cpu/test_minimize_evaluate_set_api.cpp"

  - id: ac-006
    summary: Full check + test suite green
    type: code
    pass_when: "`cmake --build build -j$(nproc)` and `ctest --test-dir build --output-on-failure` both succeed"
    test_path: "ctest --test-dir build"

  - id: ac-007
    summary: CLAUDE.md Example block advertises set_* API and auto-injection rules
    type: docs
    pass_when: "CLAUDE.md '## Example' section uses set_pair / set_neighbor / set_integrator / set_thermo, NPT example explicitly notes 'VirialStress auto-injected via reads()==Stress', AIMD example uses set_integrator + add(SCFEvaluator) for the not-yet-categorical case"
    evaluator_hint: "manual review of CLAUDE.md"
---

# Acceptance: driver-infra-auto-injection

LAMMPS-style categorical setter API + capability-driven auto-injection
landed on three CPU drivers (`cpu::md::MD`, `cpu::minimizer::Minimize`,
`cpu::Evaluate`).

All `type: code` criteria are gated by `/mol:impl`'s test runs.
`type: docs` is verified by inspection at close-out.

T8 (test/example/tutorial migration to set_* surface) and T13 (MPI
compile-time gating CI matrix) are explicitly deferred to follow-up
specs and are NOT criteria of this spec — see the parent spec's
"Out of scope" section.
