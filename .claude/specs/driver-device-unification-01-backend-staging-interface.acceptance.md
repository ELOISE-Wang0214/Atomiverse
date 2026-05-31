---
slug: driver-device-unification-01-backend-staging-interface
criteria:
  - id: ac-001
    summary: cpu::Backend exposes stage_in / stage_out / drain delegating to existing ops
    type: code
    pass_when: |
      `include/atomiverse/cpu/backend.h` declares `void stage_in(cpu::Frame&)`,
      `void stage_out(cpu::Frame&)`, `void drain()`. `src/cpu/backend.cpp`
      implements them as: stage_in → set_frame(frame); stage_out → no-op
      (zero-copy); drain → drain_analysis(). The old public names remain (this
      link is additive). Project builds via `cmake --build build -j$(nproc)`.
    status: pending
  - id: ac-002
    summary: cuda::Backend exposes stage_in / stage_out / drain delegating to existing ops
    type: code
    pass_when: |
      `include/atomiverse/cuda/backend.cuh` declares `stage_in(cpu::Frame&)`,
      `stage_out(cpu::Frame&)`, `drain()`. `src/cuda/driver.cu` implements them
      as: stage_in → materialize(frame); stage_out → retrieve(frame); drain →
      sync(). Old public names remain. CUDA target builds clean.
    status: pending
  - id: ac-003
    summary: cpu::ComponentBase declares provides_neighbor_list(); NeighborList overrides true
    type: code
    pass_when: |
      `include/atomiverse/cpu/component.h` declares
      `virtual bool provides_neighbor_list() const` defaulting to `false`;
      `cpu::md::NeighborList` (`include/atomiverse/cpu/md/neighbor_list.h`)
      overrides it to return `true`. Surface matches `cuda::ComponentBase`.
    status: pending
  - id: ac-004
    summary: New backend-staging test passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_backend_staging --output-on-failure`
      exits 0. The test asserts (a) a CPU build→stage_in→step→stage_out→drain
      run on a small LJ frame produces the same final positions/energy as the
      legacy set_frame path, and (b) `cpu::md::NeighborList{}.provides_neighbor_list()`
      is true while a bare component stub is false.
    status: pending
  - id: ac-005
    summary: Full existing suite stays green (zero behavior change)
    type: runtime
    pass_when: |
      `ctest --test-dir build --output-on-failure` exits 0 — `test_md`,
      `cuda_md*`, minimizer, evaluate, vib/hessian unchanged. This link adds
      methods only; no driver call site is modified and no public name removed.
    status: pending
---
