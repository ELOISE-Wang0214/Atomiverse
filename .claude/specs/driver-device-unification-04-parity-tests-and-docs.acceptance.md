---
slug: driver-device-unification-04-parity-tests-and-docs
criteria:
  - id: ac-001
    summary: Plan-equivalence parity test passes across CPU and CUDA
    type: runtime
    pass_when: |
      `ctest -R test_driver_device_parity --output-on-failure` exits 0. The test
      builds CPU and CUDA `md::MD` (and Evaluate/Minimize) from identical
      categorical setter sequences and asserts the resulting Plans have an
      identical (phase, order, role, capabilities) sequence.
    status: pending
  - id: ac-002
    summary: Validation parity across backends
    type: runtime
    pass_when: |
      Within the same parity test: missing-Integrator, missing-Force,
      missing-Evaluator(+capabilities), and missing-NeighborList configurations
      each throw the same category of exception on both CPU and CUDA — proving a
      single shared validate_components()/install_defaults().
    status: pending
  - id: ac-003
    summary: CUDA molrs.Frame overload round-trips via the category template
    type: runtime
    pass_when: |
      `ctest -R test_cuda_molrs_overload --output-on-failure` exits 0 —
      materialize a molrs.Frame, run one `cuda::Evaluate.run(molrs::FrameRef&)`,
      write back; the overload (now on the device-generic template) compiles,
      links, and round-trips for the CUDA instantiation.
    status: pending
  - id: ac-004
    summary: CLAUDE.md taxonomy diagram updated to the single-template stack
    type: docs
    pass_when: |
      `CLAUDE.md` ## Architecture § Driver taxonomy shows
      `Driver<P,DT>` → `BackendDriver<P,DT>` →
      `{Evaluate,Dynamics,Optimize,Response}Driver<DT>`, with no remaining
      `cpu::CPUDriver` / `cuda::CUDADriver` per-backend layer description.
    status: pending
  - id: ac-005
    summary: Full suite green; review gate clean
    type: runtime
    pass_when: |
      `ctest --test-dir build --output-on-failure` exits 0 with the new parity
      tests included. A `/mol:review` pass on the full chain diff reports 0 🚨
      and no unaddressed 🔴.
    status: pending
---
