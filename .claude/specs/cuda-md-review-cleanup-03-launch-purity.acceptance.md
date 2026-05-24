---
slug: cuda-md-review-cleanup-03-launch-purity
criteria:
  - id: ac-001
    summary: PairLJCut::compute() free of per-step cudaStreamSynchronize and D2H
    type: code
    pass_when: |
      `tests/cuda/test_pair_lj_no_d2h.cu` static-grep gate passes:
      `PairLJCut::compute()` body (with the gated first-launch
      `pair_capacity_host_ == 0` block stripped) contains zero
      `cudaStreamSynchronize(` and zero `cudaMemcpyDeviceToHost`
      occurrences. The first-launch read of `pair_capacity_host_`
      runs once per simulation (build-time-equivalent) and the kernel
      already early-exits via `if (idx < n_pairs)` line 63. Verified
      2026-05-24.
    status: verified
    last_checked: 2026-05-24
  - id: ac-002
    summary: NoseHooverNPT consolidates reductions into a single 8-double pinned buffer + defer-stream + event-sync
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu` declares `gather_npt_scalars_k`
      (mirror stress[0..2] + energy → scratch[4..7]) paired with
      `reduce_ke_axes_k` (fills scratch[0..3]); both kernels write
      into `State::d_npt_scratch_` (8 doubles, gated by
      `kNptScratchSize` constexpr). `refresh_npt_scalars_` runs both
      kernels on defer_stream, async-copies the 8-double scratch to
      `h_npt_scratch_` (pinned), records `npt_scratch_event_`, and
      cudaEventSynchronize-waits on the host thread (no
      cudaStreamSynchronize on main_stream). Old per-helper D2H
      patterns at the original lines 250, 268, 399, 507-512, 542-545
      are removed; readers consume `h_npt_scratch_[0..7]` instead.
    status: verified
    last_checked: 2026-05-24
  - id: ac-003
    summary: NPT box-diagonal D2H gated by device event (no main-stream sync)
    type: code
    pass_when: |
      `propagate_x_and_box_` calls `refresh_box_diag_(ctx)` instead
      of a synchronous `cudaMemcpy` triple. The helper runs
      `gather_box_diag_k` on defer_stream, async-copies 3 doubles to
      `h_box_diag_pinned_`, records `box_diag_event_`, and
      cudaEventSynchronize-waits. Picks up cooperative Deform-stream
      box mutations without draining main_stream. `kBoxDiagSize`
      constexpr documents the layout.
    status: verified
    last_checked: 2026-05-24
  - id: ac-004
    summary: NHC chain is a single fused kernel launch per nhc_atoms_half_ call
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu` `nhc_atoms_half_` and
      `nhc_baro_half_` launch exactly one kernel (`scale_vel_uniform_k`)
      per call; the NHC chain is integrated entirely on host (eta_t,
      eta_t_dot host-side loops, std::exp host-side). No mtchain-many
      interleaved launches — the spec's original phrasing matched a
      pre-existing rewrite. No further fusion needed.
    status: verified
    last_checked: 2026-05-24
  - id: ac-005
    summary: Backend::step() free of per-phase cudaStreamSynchronize
    type: code
    pass_when: |
      `tests/cuda/test_driver_no_per_phase_sync.cu` static-grep gate
      passes: `Backend::step()` body contains 0 `cudaStreamSynchronize(`
      calls. Intra-stream FIFO + AsyncBridge's event-based
      cross-stream ordering replace the previous double-stream sync.
      `Backend::sync()` and `retrieve()` may still sync at end-of-run.
      Verified 2026-05-24.
    status: verified
    last_checked: 2026-05-24
  - id: ac-006
    summary: AsyncBridge uses cudaLaunchHostFunc + drain-on-destroy
    type: code
    pass_when: |
      `tests/cuda/test_driver_no_per_phase_sync.cu` second gate:
      `async_bridge.cu` contains ≥1 `cudaLaunchHostFunc(` call (in
      `AsyncBridge::launch`) and at most 1 `cudaStreamSynchronize(`
      (the destroy-time drain to flush pending worker-thread
      callbacks before pinned buffers free). The previous mid-step
      `cudaStreamSynchronize(defer_stream)` is removed. Callbacks
      receive a heap-allocated `CallbackPayload` capturing the
      LaunchContext scalar fields + frame/runtime pointers; freed
      after the user hook returns. Verified 2026-05-24.
    status: verified
    last_checked: 2026-05-24
  - id: ac-007
    summary: ZeroMomentum runs without D2H+sync (on-device v_cm)
    type: code
    pass_when: |
      `src/cuda/md/init.cu` `ZeroMomentum::launch()` no longer
      contains `cudaStreamSynchronize` or D2H. A new
      `remove_com_from_scratch` kernel reads the 4-double device
      scratch (Σmv + Σm), computes v_cm = Σmv / Σm per thread, and
      subtracts directly from `frame.atoms.v[x,y,z]`. No host
      barrier on the per-step path.
    status: verified
    last_checked: 2026-05-24
  - id: ac-008
    summary: extended_hamiltonian_ computed without per-step D2H/H2D
    type: code
    pass_when: |
      `extended_hamiltonian_` reads E_pot from `h_npt_scratch_[7]`
      (populated by the most-recent `refresh_npt_scalars_` call in
      end_of_step P8 / P9) and KE from `KE_global_` (host-side,
      analytically maintained). The previous 2 cudaMemcpy(D2H) reads
      of energy + kinetic are removed. End-of-step H2D writes of
      `kinetic` and `conserved` go through the pinned
      `h_kinetic_conserved_pinned_[2]` buffer (source must outlive
      the queued async copy — a stack-local double would race).
    status: verified
    last_checked: 2026-05-24
  - id: ac-009
    summary: Extended hot-path test passes with zero wrapper-sync over a 100-step NPT loop
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_hot_path_no_sync --output-on-failure`
      exits 0. The new `ZeroSyncOverNPTLoop` test drives a 100-step
      `cuda::md::MD` loop with `NoseHooverNPT` + `NeighborList` +
      `PairLJCut` + `InitVelocity` + a `NoOpAsyncRecorder`, and
      asserts `cuda::detail::sync_count() == 0` (wrapper counts only
      `cuda::stream_synchronize` calls; component-level
      `cudaStreamSynchronize` is excluded from the gate by design
      per the link-07 sync-guard rule). Bonded evaluators
      (HarmonicBond, HarmonicAngle, MultiHarmonicDihedral) were
      pre-existing components in the framework but exercising them
      requires constructing topology — covered by
      `test_cuda_potentials_bonded` as a separate gate. Verified
      2026-05-24: ZeroSyncOver100Steps + ZeroSyncOverNPTLoop +
      LifecycleRoundTrip all pass.
    status: verified
    last_checked: 2026-05-24
  - id: ac-010
    summary: Driver per-phase sync gate test passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_cuda_driver_no_per_phase_sync
      --output-on-failure` exits 0. Both static-grep tests pass:
      `Backend::step()` body has 0 cudaStreamSynchronize calls;
      `async_bridge.cu` has 1 cudaStreamSynchronize (destroy drain)
      and ≥1 cudaLaunchHostFunc. Verified 2026-05-24.
    status: verified
    last_checked: 2026-05-24
  - id: ac-011
    summary: NPT correctness preserved (existing test passes)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_cuda_md_nose_hoover_npt
      --output-on-failure` exits 0. The defer-stream + event-sync
      consolidation does not destabilize NHC feedback (per-call
      cudaEventSynchronize keeps the host thread in step with
      device state — equivalent host wait semantics to the old
      cudaStreamSynchronize but without draining the main compute
      stream). Verified 2026-05-24: full 23/23 cuda regression
      passed including the existing NPT test.
    status: verified
    last_checked: 2026-05-24
