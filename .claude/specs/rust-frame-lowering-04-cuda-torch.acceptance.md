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
    status: pending
  - id: ac-002
    summary: the three execution Frames are structurally aligned
    type: code
    pass_when: |
      cuda::Frame and torch::Frame embed core/system_scalars.h SystemScalars and
      carry tag, image, and conserved fields, matching cpu::Frame.
    status: pending
  - id: ac-003
    summary: cuda / torch materializers round-trip a molrs.Frame
    type: runtime
    pass_when: |
      When the CUDA toolkit is present, a molrs.Frame -> cuda::Frame -> write-back
      round trip is element-wise equal (and likewise for torch when
      ATV_BUILD_TORCH=ON); on a CPU-only build the rest of `ctest` is green.
    status: pending
  - id: ac-004
    summary: full test suite stays green
    type: runtime
    pass_when: |
      `cmake --build build` succeeds and `ctest --test-dir build` reports all
      tests passing; CUDA tests pass when the toolkit is available.
    status: pending
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
