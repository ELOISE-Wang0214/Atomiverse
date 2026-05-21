---
slug: phase4-vibrations-01-numerical-hessian
criteria:
  - id: ac-001
    summary: Response category symbols present in source
    type: code
    pass_when: |
      grep finds `class ResponseDriver` in
      include/atomiverse/cpu/response_driver.h, `enum class Phase` inside
      `namespace vib` in include/atomiverse/facade/stage.h,
      `CategoryTraits<ATV::vib::Phase>` in include/atomiverse/cpu/cpu_driver.h,
      and `Hessian` (value 1u << 7) in the Capability enum in
      include/atomiverse/facade/component_role.h.
    status: pending
  - id: ac-002
    summary: vib domain headers and CMake registration present
    type: code
    pass_when: |
      include/atomiverse/cpu/vib/hessian.h, hessian_result.h,
      reference_geometry.h, displace_step.h, hessian_accumulator.h and
      symmetrize.h all exist; CMakeLists.txt lists src/cpu/vib/hessian.cpp in
      ATOMIVERSE_SRCS; cpu/vib/hessian.h contains no `#include` of any
      cpu/scf/ header.
    status: pending
  - id: ac-003
    summary: Project builds with the new Response category
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON &&
      cmake --build build -j` completes with no errors; the
      test_cpu_vib_numerical_hessian target is produced.
    status: pending
  - id: ac-004
    summary: Numerical Hessian matches in-test gradient-difference oracle
    type: runtime
    evaluator_hint: scientific
    pass_when: |
      In test_cpu_vib_numerical_hessian, the cpu::vib::Hessian-produced
      Hessian for H2O (RHF or PBE, STO-3G or 6-31G) and an oracle that
      central-differences scf_gradient directly with the same h=1.0e-3 A,
      same method and same convergence agree elementwise to <= 1e-6 eV/Angstrom^2.
    status: pending
  - id: ac-005
    summary: Returned Hessian is symmetric
    type: runtime
    pass_when: |
      In test_cpu_vib_numerical_hessian, the HessianResult.hessian for the
      H2O case satisfies max_{iα,jβ} |H[iα,jβ] - H[jβ,iα]| <= 1e-12.
    status: pending
  - id: ac-006
    summary: Hessian satisfies the translational sum rule near a stationary point
    type: runtime
    evaluator_hint: scientific
    pass_when: |
      In test_cpu_vib_numerical_hessian, at a near-stationary H2O geometry,
      for each axis α the translational-block row sum
      |Σ_j H[iα,jα]| <= 1e-4 eV/Angstrom^2 for every atom i.
    status: pending
  - id: ac-007
    summary: ResponseDriver validates required components
    type: runtime
    pass_when: |
      In test_cpu_vib_numerical_hessian, building a cpu::vib::Hessian with no
      Capability::Force evaluator throws, and building one with the
      HessianAccumulator removed (no Capability::Hessian component) throws;
      a correctly configured driver builds without throwing.
    status: pending
  - id: ac-008
    summary: Single-atom frame runs without crash
    type: runtime
    pass_when: |
      In test_cpu_vib_numerical_hessian, running cpu::vib::Hessian on an
      N=1 frame completes 6 displacement passes and returns a HessianResult
      whose hessian has shape [3,3].
    status: pending
  - id: ac-009
    summary: Full test suite stays green
    type: runtime
    pass_when: |
      `ctest --test-dir build --output-on-failure` reports 100% of tests
      passed, including test_cpu_vib_numerical_hessian.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 (code)** — Static presence of the fourth category, its phase enum, the `CategoryTraits` specialization, the `Capability::Hessian` bit, the new `cpu/vib/` headers, the CMake source registration, and the no-`cpu/scf/`-edge layering rule.
- **ac-003 (runtime)** — The project compiles with the new category and capability bit; `CategoryTraits<vib::Phase>` is required for `CPUDriver<vib::Phase>` to instantiate, so a clean build is a real gate.
- **ac-004 (runtime, scientific)** — The self-consistency proof: the driver's `6N`-displacement assembly equals an oracle that differences `scf_gradient` directly under identical conditions. This is the binding correctness bar for the finite-difference Hessian.
- **ac-005 (runtime)** — `Symmetrize` correctly projects onto the symmetric subspace.
- **ac-006 (runtime, scientific)** — A cheap physical sanity check (translational-invariance sum rule) standing in for the full 6-near-zero-eigenvalue analysis, which is deferred to link 2.
- **ac-007 (runtime)** — `ResponseDriver::validate_components()` enforces both required ingredients: a `Capability::Force` evaluator and a `Capability::Hessian` accumulator.
- **ac-008 (runtime)** — Degenerate `N = 1` frame exercises the `6N` loop bound without special-casing.
- **ac-009 (runtime)** — No regression from adding a new outer category and a new capability bit.
