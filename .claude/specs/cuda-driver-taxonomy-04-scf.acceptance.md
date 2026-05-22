---
slug: cuda-driver-taxonomy-04-scf
criteria:
  - id: ac-001
    summary: No cuda::scf::SCF driver class and no cuda::scf::SinglePoint
    type: code
    pass_when: |
      No type named `SCF` exists in namespace `ATV::cuda::scf` and no type
      named `SinglePoint` exists in `ATV::cuda::scf`; `include/atomiverse/cuda/scf/scf.cuh`
      and `src/cuda/scf/driver.cu` are deleted; a grep for `cuda::scf::SCF`
      as a driver class returns no production-code matches.
    status: pending
  - id: ac-002
    summary: cuda::scf::SCFSolver exists and is not a Driver subclass
    type: code
    pass_when: |
      `ATV::cuda::scf::SCFSolver` is declared in
      `include/atomiverse/cuda/scf/scf_solver.cuh`, and a
      `static_assert(!std::is_base_of_v<ATV::Driver<ATV::scf::Phase,
      ATV::cuda::device_tag>, ATV::cuda::scf::SCFSolver>)` (or equivalent
      type-level check) compiles, proving it is not a Driver subclass; it
      exposes a build/run lifecycle with no `scf::Phase` parameter.
    status: pending
  - id: ac-003
    summary: SCF step components live in cuda::scf::internal:: and are not Plan-registered
    type: code
    pass_when: |
      Hamiltonian, EigenSolver, DensityBuilder, RHFFock, RHFFockDevice,
      RHFDensity, ConvergenceCheck are declared in namespace
      `ATV::cuda::scf::internal` under `include/atomiverse/cuda/scf/internal/`,
      none inherits `ATV::cuda::ComponentBase`, none defines
      `default_bindings()`, and no facade-public production header includes
      a `cuda/scf/internal/` header.
    status: pending
  - id: ac-004
    summary: cuda::Evaluate is a one-pass driver inheriting cuda::EvaluateDriver
    type: code
    pass_when: |
      `ATV::cuda::Evaluate` is declared in `include/atomiverse/cuda/evaluate.cuh`
      and `static_assert(std::is_base_of_v<ATV::cuda::EvaluateDriver,
      ATV::cuda::Evaluate>)` compiles; `ATV::cuda::single_point` is a free
      function template in `include/atomiverse/cuda/single_point.cuh`, not a class.
    status: pending
  - id: ac-005
    summary: cuda::scf::SCFEvaluator is an Evaluator-role component with Energy|DensityMatrix|Orbitals
    type: code
    pass_when: |
      `ATV::cuda::scf::SCFEvaluator` inherits `ATV::cuda::ComponentBase`,
      `role()` returns `ComponentRole::Evaluator`, and `capabilities()`
      returns exactly `Capability::Energy | Capability::DensityMatrix |
      Capability::Orbitals` (no `Force`); it owns an `SCFSolver` member.
    status: pending
  - id: ac-006
    summary: Single-point CUDA RHF energy via cuda::Evaluate + SCFEvaluator reproduces the reference
    type: runtime
    evaluator_hint: requires CUDA toolkit + GPU
    pass_when: |
      A test running single-point H2/STO-3G RHF through
      `cuda::Evaluate + cuda::scf::SCFEvaluator` (or `cuda::single_point`)
      writes a frame energy equal to the textbook reference -1.1167 Ha
      within 0.01 Ha, matching the pre-migration `test_cuda_scf_rhf` value.
    status: pending
  - id: ac-007
    summary: The full CUDA SCF regression suite passes after migration
    type: runtime
    evaluator_hint: requires CUDA toolkit + GPU
    pass_when: |
      `ctest --test-dir build -R cuda_scf --output-on-failure` reports all
      of test_cuda_scf, test_cuda_scf_rhf, test_cuda_scf_device_integrals
      passing, with every energy matching its pre-migration value to the
      suite's existing tolerance.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003 / ac-004 / ac-005** — structural code checks that the SCF migration matches the chain's fixed decisions: the SCF driver class and any `SinglePoint` class are gone, `SCFSolver` is a non-`Driver` inner solver, step components are de-`Component`-ised into `cuda::scf::internal::`, and the new `cuda::Evaluate` driver / `cuda::single_point` helper / `cuda::scf::SCFEvaluator` component exist with the correct base classes, role, and capabilities. These are verifiable by reading headers and compiling type-level assertions; no GPU needed.
- **ac-006** — the binding numeric-preservation check for the new single-point path: the user-facing `cuda::Evaluate + cuda::scf::SCFEvaluator` composition produces the same RHF energy the old `cuda::scf::SCF` driver did. Requires a CUDA toolkit and GPU.
- **ac-007** — the global regression gate: the entire pre-existing CUDA SCF test suite still passes, proving the structural refactor changed no SCF numerics. Requires a CUDA toolkit and GPU.
