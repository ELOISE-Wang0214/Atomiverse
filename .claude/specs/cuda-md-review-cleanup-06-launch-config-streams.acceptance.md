---
slug: cuda-md-review-cleanup-06-launch-config-streams
criteria:
  - id: ac-001
    summary: __launch_bounds__ annotated on all hot kernels
    type: code
    pass_when: |
      `grep -nE '__launch_bounds__' src/cuda/potentials/pair_lj_cut.cu
      src/cuda/potentials/bond_harmonic.cu src/cuda/potentials/angle_harmonic.cu
      src/cuda/potentials/dihedral_multi_harmonic.cu src/cuda/md/kinetic_energy.cu
      src/cuda/observer/stress.cu` returns at least one hit per file. Default
      is `(256, 4)`; `dihedral_multi_harmonic_k` uses `(128, 4)`.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: dihedral_multi_harmonic_k register usage within budget
    type: code
    pass_when: |
      `nvcc --ptxas-options=-v -res-usage` build output (or
      `cuobjdump --dump-resource-usage`) reports the
      `dihedral_multi_harmonic_k` kernel uses ≤ 128 registers per thread
      with spill ≤ 128 bytes. Captured in a build-log artifact or
      docstring comment.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: cuda::Context exposes force_streams_ pool + events
    type: code
    pass_when: |
      `include/atomiverse/cuda/context.cuh` `cuda::Context` declares
      `std::array<cudaStream_t, 4> force_streams_`, `cudaEvent_t
      force_start_event_`, and `std::array<cudaEvent_t, 4> force_done_events_`
      (or equivalent). Created in `Context::build()`, destroyed in
      `Context::destroy()`.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: Backend::dispatch_force_phase routes evaluators to dedicated streams
    type: code
    pass_when: |
      `src/cuda/driver.cu` Backend::dispatch_force_phase (or equivalent)
      uses `cudaEventRecord(force_start_event_, main_stream)` then
      `cudaStreamWaitEvent(force_streams_[i], force_start_event_)` before
      each force-evaluator's `launch()`. End-of-phase events fan-in to
      `main_stream` for the post_force phase.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: LJ uses __drcp_rn (or __frcp_rn for float) instead of FP64 division
    type: code
    pass_when: |
      `grep -nE '__drcp_rn|__frcp_rn' src/cuda/potentials/pair_lj_cut.cu`
      returns at least 1 hit. `grep -n 'work_t\(1\) / r\|1.0 / r'`
      inside the LJ kernel body returns 0 hits.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: export_atom_pairs_kernel stores r² (no sqrt)
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/export_pairs.cu` `export_atom_pairs_kernel`
      writes `r2` (squared distance) to the pair buffer; `grep -n 'sqrt'`
      inside the kernel body returns 0 hits. LJ kernel consumes r² directly.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: Correctness preserved
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials' --output-on-failure`
      exits 0. Force / energy / stress / kinetic-energy values within
      existing tolerances.
    status: pending
    last_checked: ""
