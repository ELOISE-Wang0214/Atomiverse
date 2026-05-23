---
slug: msu-pe-cuda-components-02-cuda-categorical-api
criteria:
  - id: ac-001
    summary: cuda::md::MD exposes 11 categorical setters matching cpu::md::MD
    type: code
    pass_when: |
      `include/atomiverse/cuda/md/md.cuh` declares public member functions
      `set_pair`, `set_bond`, `set_angle`, `set_dihedral`, `set_neighbor`,
      `set_integrator`, `set_thermostat`, `set_barostat`, `set_observer`,
      `set_record`, `set_thermo`; every name also present in
      `include/atomiverse/cpu/md/md.h` (directly) or `include/atomiverse/cpu/cpu_driver.h`
      (inherited). A grep of the CPU set against the CUDA header returns
      11/11 matches.
    status: pending
  - id: ac-002
    summary: Singular setters reject mismatched re-call, accept same-name re-call
    type: runtime
    pass_when: |
      `tests/cuda/md/test_cuda_md_components.cpp::CudaMdCategoricalApi.SingularRejectsMismatch`
      and `CudaMdCategoricalApi.SingularIdempotentOnSame` pass; the
      throwing case raises `std::invalid_argument` whose `what()` contains
      both component names; the idempotent case results in exactly one
      pending binding for the slot.
    status: pending
  - id: ac-003
    summary: Plural setters append in insertion order
    type: runtime
    pass_when: |
      `CudaMdCategoricalApi.PluralAppends` registers three observers /
      recorders / thermostats via `set_observer` / `set_record` /
      `set_thermostat` and asserts the pending order matches insertion
      order across all three.
    status: pending
  - id: ac-004
    summary: VirialStress auto-injected when a Barostat-role component is present
    type: runtime
    pass_when: |
      `CudaMdCategoricalApi.VirialStressAutoInjectedOnBarostat` builds an
      MD with `set_barostat(cuda::md::Barostat{...})` (no explicit
      VirialStress) and asserts a component with `name()=="virial_stress"`
      is in the pending list after `build()`.
      `CudaMdCategoricalApi.VirialStressDedupedWhenExplicit` shows
      exactly one `virial_stress` when the user pre-registers one.
      `CudaMdCategoricalApi.NoAutoInjectWithoutBarostat` shows zero
      `virial_stress` when no barostat is set.
    status: pending
  - id: ac-005
    summary: Missing NeighborList with a pair evaluator throws at build()
    type: runtime
    pass_when: |
      `CudaMdCategoricalApi.MissingNeighborThrows` calls
      `sim.set_pair(cuda::LJ{...})` without `set_neighbor()`; `sim.build()`
      throws `std::runtime_error` whose `what()` names the evaluator and
      suggests `set_neighbor`.
    status: pending
  - id: ac-006
    summary: Categorical-API smoke MD run completes on CUDA without NaN
    type: runtime
    pass_when: |
      `CudaMdCategoricalApi.SmokeAllSetters` builds a 64-atom LJ system
      entirely via `set_pair` / `set_neighbor` / `set_integrator` /
      `set_thermo`, calls `build()` + `run(frame, 10)`, and asserts
      `std::isfinite(frame.scalars.energy)` and
      `std::isfinite(frame.scalars.kinetic)`.
    status: pending
  - id: ac-007
    summary: Raw add() / on(phase).add() paths still work — zero regression
    type: runtime
    pass_when: |
      Existing `tests/cuda/md/test_md.cu` MD.NVEEnergyConservation,
      MD.ThermoOutput, MDTaxonomy.IsADynamicsDriver,
      MDTaxonomy.ThrowsOnMissingIntegrator,
      MDTaxonomy.ThrowsOnMissingForceWriter continue to pass unchanged.
    status: pending
  - id: ac-008
    summary: New test file registered in CMakeLists.txt
    type: code
    pass_when: |
      `CMakeLists.txt` lists `tests/cuda/md/test_cuda_md_components.cpp`
      in the CUDA tests target so a fresh `cmake --build` followed by
      `ctest` runs the new `CudaMdCategoricalApi.*` cases.
    status: pending
  - id: ac-009
    summary: All new public API has doxygen with units where applicable
    type: docs
    pass_when: |
      Every new public method in `include/atomiverse/cuda/md/md.cuh`
      carries a `///` doxygen block with `@brief` and (where the parameter
      is a physical quantity) units in square brackets, matching CPU
      style in `include/atomiverse/cpu/md/md.h`.
    status: pending
  - id: ac-010
    summary: Full check + test suite green
    type: runtime
    pass_when: |
      The project's full check + test suite (per CLAUDE.md `mol_project`
      build target) exits 0 with all previously-green tests plus the new
      `CudaMdCategoricalApi.*` cases passing.
    status: pending
---

# Acceptance criteria

ac-001 / ac-008 / ac-009 are static (file / build-system / doc inspection). ac-002 through ac-007 and ac-010 are runtime gtest assertions. ac-007 is the explicit zero-regression contract — additive setter API must not break the existing raw `add()` / `on(phase).add()` paths exercised by `tests/cuda/md/test_md.cu`. ac-004 is the load-bearing semantic — `VirialStress` auto-injection must work via `role()==Barostat`, since `cuda::ComponentBase` has no `reads()` virtual yet (per `cuda/md/barostat.cuh` deferral note).
