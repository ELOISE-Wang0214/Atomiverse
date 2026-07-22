---
spec: backend-mirror-04-driver
created: 2026-07-22
criteria:
  - id: ac-001
    summary: CPUDriver and CUDADriver symbols removed
    type: code
    pass_when: |
      rg '\bCPUDriver\b|\bCUDADriver\b' over include/, src/, tests/, examples/,
      regressions/, python/ returns no matches in product sources.
    status: pending
  - id: ac-002
    summary: DriverNode is the erased nestable base
    type: code
    pass_when: |
      cpu::DriverNode (and cuda::DriverNode if cuda nesting exists) is the
      non-template nestable base; phase loop class is cpu::Driver / cuda::Driver.
    status: pending
  - id: ac-003
    summary: CUDA uses Criterion and add_criterion
    type: code
    pass_when: |
      No LoopHandler or add_handler in include/atomiverse/cuda/; cuda Driver exposes
      add_criterion; ForceConverged/StopOnNaN are Criterion-shaped.
    status: pending
  - id: ac-004
    summary: driver_names regression builds and runs
    type: runtime
    pass_when: |
      regressions/layout/driver_names/ builds and runs using new Driver/criterion
      names only (no CPUDriver/LoopHandler).
    status: pending
  - id: ac-005
    summary: full suite green after rename
    type: runtime
    pass_when: |
      Full check + ctest pass.
    status: pending
out_of_scope:
  - Backend dissolve; path renames; docs-only updates
---

# Acceptance — backend-mirror-04-driver

Done means public loop types are `Driver`/`DriverNode` and termination is Criterion on both backends, with tests rewritten.
