---
slug: cuda-driver-taxonomy-02-md
criteria:
  - id: ac-001
    summary: cuda::md::MD inherits cuda::DynamicsDriver
    type: code
    pass_when: |
      include/atomiverse/cuda/md/md.cuh declares `class MD` with a
      base-specifier of `cuda::DynamicsDriver` (or `DynamicsDriver`
      resolved in the cuda namespace), not `Driver<Phase>` /
      `Driver<md::Phase, cuda::device_tag>`.
    status: pending
  - id: ac-002
    summary: cuda::md::MD no longer names the raw Driver template
    type: code
    pass_when: |
      No occurrence of `Driver<` remains in
      include/atomiverse/cuda/md/md.cuh (base clause or member
      signatures); the only driver type named is the category base.
    status: pending
  - id: ac-003
    summary: Zero* MD components are phase-agnostic
    type: code
    pass_when: |
      In include/atomiverse/cuda/md/init.cuh, none of ZeroForces,
      ZeroEnergy, ZeroKinetic override default_bindings() with a
      hardcoded md::Phase entry; each inherits the empty-vector
      ComponentBase default. cuda/md/md.cuh binds them via explicit
      on(md::Phase::...).add().
    status: pending
  - id: ac-004
    summary: Every CUDA MD-domain component declares role()/capabilities()
    type: code
    pass_when: |
      VelocityVerlet and Deform override role()==ComponentRole::Integrator;
      NoseHoover overrides role()==ComponentRole::Thermostat; Barostat
      overrides role()==ComponentRole::Barostat. InitVelocity,
      KineticEnergy and the Zero* components leave role()==None and
      capabilities()==None. Tags match the CPU MD equivalents under
      include/atomiverse/cpu/md/.
    status: pending
  - id: ac-005
    summary: MD build-time validation throws on a missing Integrator
    type: runtime
    pass_when: |
      A cuda::md::MD configured with a Force-capable evaluator
      (cuda::PairLJCut) but no Integrator-role component throws from
      build() (via cuda::DynamicsDriver::validate_components); the RED
      test asserting this in tests/cuda/md/test_md.cu passes. Requires
      a CUDA toolkit.
    status: pending
  - id: ac-006
    summary: cuda::md::MD is-a DynamicsDriver and Zero* binding test passes
    type: runtime
    pass_when: |
      The structural test in tests/cuda/md/test_md.cu — asserting
      std::is_base_of_v<cuda::DynamicsDriver, cuda::md::MD> and that
      ZeroForces/ZeroEnergy/ZeroKinetic return an empty
      default_bindings() — compiles and passes. Requires a CUDA toolkit.
    status: pending
  - id: ac-007
    summary: CUDA MD regression suite passes (NVE drift + thermostat)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R cuda_md --output-on-failure` reports
      all matched tests passing, including test_cuda_md (NVE energy
      conservation / drift bound) and test_cuda_md_thermostat
      (Nosé-Hoover target temperature). Requires a CUDA toolkit and a
      CUDA-capable device.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 (code)** — The reparent is the load-bearing change:
  `cuda::md::MD` must derive from `cuda::DynamicsDriver` and must not
  name the raw `Driver<...>` template anywhere in `md.cuh`. Together
  these confirm the inheritance chain is `MD → cuda::DynamicsDriver →
  cuda::CUDADriver<md::Phase>`.
- **ac-003 (code)** — Latent-bug fix #2: the Zero* components carry no
  hardcoded `md::Phase` binding, so 03-minimizer can reuse them against
  `minimizer::Phase`. `cuda::md::MD` re-supplies the `md::Phase`
  bindings explicitly in its constructor.
- **ac-004 (code)** — Component tagging mirrors the CPU MD components
  exactly. The force-evaluator bases (`cuda::Pair`/`Bond`/`Angle`/
  `Dihedral`) already tag `Evaluator` + `Energy|Force|Stress`, so this
  criterion concerns only the MD-domain components.
- **ac-005 (runtime)** — `cuda::DynamicsDriver::validate_components()`
  rejects a malformed MD setup at build time, proving the validation is
  wired and the `role()` tags are observed.
- **ac-006 (runtime)** — The compile-time / structural assertions in the
  RED test pass after the reparent + phase-agnostic Zero* land.
- **ac-007 (runtime)** — Behaviour-preservation gate: the unchanged
  physics produces unchanged NVE drift and thermostat behaviour.
