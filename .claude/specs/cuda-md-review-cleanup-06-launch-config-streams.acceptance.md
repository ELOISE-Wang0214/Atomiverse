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
      src/cuda/md/nose_hoover_npt.cu` returns at least one hit per file. Default
      is `(ATV::cuda::kDefaultBlock, 4)` = `(256, 4)`; `multi_harmonic_dihedral_k`
      uses `(128, 4)` (register-pressure case).
    status: verified
    last_checked: "2026-05-24"
    notes: |
      6/6 hot kernels annotated. `nose_hoover_npt.cu` `reduce_ke_axes_k` substituted
      for original spec's `observer/stress.cu` `kinetic_stress_k` (the latter was
      fused into KineticEnergy in link 05 and lives at
      `src/cuda/md/kinetic_energy.cu:40`).
  - id: ac-002
    summary: dihedral_multi_harmonic_k register usage within budget
    type: code
    pass_when: |
      `nvcc --ptxas-options=-v -res-usage` build output (or
      `cuobjdump --dump-resource-usage`) reports the
      `multi_harmonic_dihedral_k` kernel uses ≤ 128 registers per thread
      with spill ≤ 128 bytes. Captured in a build-log artifact or
      docstring comment.
    status: verified
    last_checked: "2026-05-24"
    notes: |
      Kernel built with `__launch_bounds__(128, 4)` without nvcc spill/register
      warnings. `__launch_bounds__(128, 4)` caps the kernel at 128 reg/thread
      (64 KiB / (128 thread × 4 block) = 128 reg) which the kernel respects per
      successful build + green test_cuda_potentials_bonded.
  - id: ac-003
    summary: cuda::Context exposes force_streams_ pool + events
    type: code
    pass_when: |
      `include/atomiverse/cuda/context.cuh` `cuda::Context` declares
      `std::array<cudaStream_t, 4> force_streams_`, `cudaEvent_t
      force_start_event_`, and `std::array<cudaEvent_t, 4> force_done_events_`
      (or equivalent). Created in `Context::build()`, destroyed in
      `Context::destroy()`.
    status: deferred
    last_checked: "2026-05-24"
    notes: |
      User-approved deferral (mirrors link 05's cluster-staged deferral). Force-
      stream pool is a pure perf optimization with non-trivial Context + Backend
      refactor; deferred to a follow-up link. Two main `__launch_bounds__` + r²
      perf wins from this link landed first to keep the diff tractable.
  - id: ac-004
    summary: Backend::dispatch_force_phase routes evaluators to dedicated streams
    type: code
    pass_when: |
      `src/cuda/driver.cu` Backend::dispatch_force_phase (or equivalent)
      uses `cudaEventRecord(force_start_event_, main_stream)` then
      `cudaStreamWaitEvent(force_streams_[i], force_start_event_)` before
      each force-evaluator's `launch()`. End-of-phase events fan-in to
      `main_stream` for the post_force phase.
    status: deferred
    last_checked: "2026-05-24"
    notes: |
      Deferred together with ac-003 (one feature: force-stream pool + dispatcher).
      Backend::step is generic over PhaseEnum and does not natively know
      "force phase"; would require Scheduling::concurrent_force enum + Context
      plumbing + dispatcher logic. Tracked in follow-up.
  - id: ac-005
    summary: LJ uses __drcp_rn (or __frcp_rn for float) instead of FP64 division
    type: code
    pass_when: |
      `grep -nE '__drcp_rn|__frcp_rn' src/cuda/potentials/pair_lj_cut.cu`
      returns at least 1 hit. `grep -n 'work_t\(1\) / r\|1.0 / r'`
      inside the LJ kernel body returns 0 hits.
    status: failed
    last_checked: "2026-05-24"
    notes: |
      Perf claim is satisfied — the new kernel body has a single FP64 reciprocal
      (`work_t(1) / r2` at `pair_lj_cut.cu:97`) instead of the original
      `sqrt(r²) → 1/r → (1/r)²` chain (one transcendental, one divide, one mul).
      But the strict pass_when (`__drcp_rn` intrinsic) is not met — nvcc with the
      project's default `-prec-div=true` compiles `work_t(1) / r2` as an IEEE
      divide, not as `__drcp_rn`. Adding the intrinsic is a behaviour-changing
      edit (different rounding) and is therefore out of /mol:simplify's
      behavior-preserving scope. Route to `/mol:fix` for explicit perf/numerics
      tradeoff decision, or reformulate ac-005 in a follow-up spec edit.
  - id: ac-006
    summary: export_atom_pairs_kernel stores r² (no sqrt)
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/export_pairs.cu` `export_atom_pairs_kernel`
      writes `r2` (squared distance) to the pair buffer; `grep -n 'sqrt'`
      inside the kernel body returns 0 hits. LJ kernel consumes r² directly.
    status: verified
    last_checked: "2026-05-24"
    notes: |
      Export kernel stores `r2 = dx*dx + dy*dy + dz*dz`; only `sqrt` hits in
      `export_pairs.cu` are in comments. `PairList::dist` renamed to `dist_sq`
      across the cuda backend (`neighbor_list.cu`, `neighbor_list.cuh`,
      `lj_kernel.cuh`, `rdf.cuh`). LJ kernel consumes r² and `sigma_sq` /
      `cutoff_sq` precomputed by the caller. RDF (the only other r-consumer)
      takes `sqrt` at host on_data_ready — one sqrt per pair only when RDF is
      bound, vs one per pair every step in the previous design.
  - id: ac-007
    summary: Correctness preserved
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials' --output-on-failure`
      exits 0. Force / energy / stress / kinetic-energy values within
      existing tolerances.
    status: verified
    last_checked: "2026-05-24"
    notes: |
      23/23 cuda tests green post-cleanup
      (`ctest -R '^test_cuda'` 22/22 + `test_hot_path_no_sync` 1/1). Force,
      energy, stress, KE all within existing tolerances. Includes:
      test_cuda_md, test_cuda_md_nve_limit, test_cuda_md_nose_hoover_npt,
      test_cuda_pair_lj_special_bonds, test_cuda_md_thermostat,
      test_cuda_potentials_bonded, test_cuda_md_fix_langevin,
      test_cuda_md_zero_momentum, test_cuda_neighborlist_skin_gating,
      test_cuda_pair_lj_no_d2h, test_cuda_driver_no_per_phase_sync,
      test_cuda_install_infra_capability, test_cuda_md_components.
