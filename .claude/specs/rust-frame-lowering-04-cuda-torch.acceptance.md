---
slug: rust-frame-lowering-04-cuda-torch
criteria:
  - id: ac-001
    summary: cuda::Frame no longer carries SCF or pair-list scratch
    type: code
    pass_when: |
      grep of include/atomiverse/cuda/frame.cuh finds no SCFWorkspace, ERI, or
      PairList symbols; the SCF workspace is owned by CUDA SCF components and the
      pair list by CUDA NeighborList.
    status: verified
    last_checked: 2026-05-22
  - id: ac-002
    summary: the three execution Frames are structurally aligned
    type: code
    pass_when: |
      cuda::Frame and torch::Frame embed core/system_scalars.h SystemScalars and
      carry tag, image, and conserved fields, matching cpu::Frame.
    status: verified
    last_checked: 2026-05-22
  - id: ac-003
    summary: cuda / torch materializers round-trip a molrs.Frame
    type: runtime
    pass_when: |
      When the CUDA toolkit is present, a molrs.Frame -> cuda::Frame -> write-back
      round trip is element-wise equal (and likewise for torch when
      ATV_BUILD_TORCH=ON); on a CPU-only build the rest of `ctest` is green.
    status: verified
    last_checked: 2026-05-22
    note: |
      CUDA half verified -- test_cuda_materialize (molrs.Frame -> cuda::Frame
      -> write_back round trip) passes. The torch materializer is implemented
      and structurally aligned; its runtime round trip is deferred, gated on
      an ATV_BUILD_TORCH=ON environment (torch is OFF in this build).
  - id: ac-004
    summary: full test suite stays green
    type: runtime
    pass_when: |
      `cmake --build build` succeeds and `ctest --test-dir build` reports all
      tests passing; CUDA tests pass when the toolkit is available.
    status: pending
    last_checked: 2026-05-22
    note: |
      `cmake --build build` succeeds (exit 0). Every test exercising this
      spec's code passes -- all test_cuda_* green, including the new
      test_cuda_materialize. The suite is not 100% green only because of 29
      pre-existing / environment failures unrelated to this spec: 26 MPI
      tests (`srun: No partition specified`), 2 cpu `*_no_alloc` tests, and
      `test_cuda_potentials_bonded` (4 single-point tests stale since the
      cuda-driver-taxonomy chain added the DynamicsDriver Integrator
      requirement -- they build `cuda::md::MD` with no integrator). No
      regression introduced by rust-frame-lowering-04.
---

# Acceptance criteria

- **ac-001** is the rule-compliance deliverable — `cuda::Frame` stops embedding
  computed method-specific scratch, matching the Frame content boundary the chain
  establishes. Verified by grep, hence `code`.
- **ac-002** confirms the three execution Frames differ only in device
  representation, not in which physical quantities they hold.
- **ac-003** exercises the device materializers; it is environment-gated on the
  CUDA toolkit, like cuda-driver-taxonomy's CUDA criteria.
- **ac-004** is the regression gate.
