---
slug: cuda-driver-taxonomy-03-minimizer
criteria:
  - id: ac-001
    summary: cuda::minimizer::SD and CG driver classes no longer exist
    type: code
    pass_when: |
      grep across include/atomiverse/cuda/minimizer/ and src/cuda/minimizer/
      finds no `class SD` and no `class CG` driver definitions;
      include/atomiverse/cuda/minimizer/cg.cuh is deleted and no file
      includes it.
    status: pending
  - id: ac-002
    summary: SteepestDescent Optimizer component replaces SteepestDescentStep
    type: code
    pass_when: |
      include/atomiverse/cuda/minimizer/sd.cuh defines
      `class SteepestDescent` with `ComponentRole role() const override`
      returning `ComponentRole::Optimizer`; the identifier
      `SteepestDescentStep` no longer appears in the cuda/ tree.
    status: pending
  - id: ac-003
    summary: cuda::minimizer::Minimize exists and inherits cuda::OptimizeDriver
    type: code
    pass_when: |
      include/atomiverse/cuda/minimizer/minimize.cuh (new) defines
      `class Minimize : public OptimizeDriver` in namespace
      ATV::cuda::minimizer, with a constructor that binds Zero* via
      on(minimizer::Phase::force).add(...).
    status: pending
  - id: ac-004
    summary: Minimize + SteepestDescent relaxes a test geometry, energy decreasing
    type: runtime
    pass_when: |
      `ctest --test-dir build -R cuda_minimizer --output-on-failure`
      passes; test_cuda_minimizer_sd reports final energy strictly less
      than the post-single-evaluation energy, matching the
      pre-migration relaxed result. Requires a CUDA toolkit.
    status: pending
  - id: ac-005
    summary: build() throws when no Optimizer component is registered
    type: runtime
    pass_when: |
      a test constructs cuda::minimizer::Minimize, adds only a force
      evaluator (no SteepestDescent), calls build(), and observes an
      exception from OptimizeDriver validation; the test passes under
      `ctest --test-dir build -R cuda_minimizer`. Requires a CUDA toolkit.
    status: pending
  - id: ac-006
    summary: Doxygen present on SteepestDescent and Minimize per doc style
    type: code
    pass_when: |
      SteepestDescent and cuda::minimizer::Minimize and their public
      members carry @brief, and @param entries state units (alpha in fs);
      no stale references to the deleted SD / CG drivers in their Doxygen.
    status: pending
---

# Acceptance criteria

- **ac-001** — The algorithm-as-driver anti-pattern is gone: both the `SD`
  and `CG` driver classes are deleted, and `cg.cuh` (an empty shell) is
  removed from the tree with no dangling includes.
- **ac-002** — The steepest-descent step is a first-class `Optimizer`
  component, renamed `SteepestDescent` for parity with CPU `Optimizer`
  component naming, and tagged with `role()==ComponentRole::Optimizer`.
- **ac-003** — The single user-facing CUDA minimizer driver
  `cuda::minimizer::Minimize` exists, subclasses `cuda::OptimizeDriver`, and
  binds the phase-agnostic `Zero*` infrastructure at `minimizer::Phase::force`
  in its constructor.
- **ac-004** — Behavior preservation: the migrated `Minimize + SteepestDescent`
  shape relaxes the perturbed LJ lattice to a lower-energy minimum, matching
  the pre-migration trajectory. Requires a CUDA toolkit.
- **ac-005** — `cuda::OptimizeDriver`'s build-time validation rejects a
  `Minimize` with no `Optimizer` component. Requires a CUDA toolkit.
- **ac-006** — Public symbols carry project-standard Doxygen with units, and
  no documentation references the deleted driver classes.
