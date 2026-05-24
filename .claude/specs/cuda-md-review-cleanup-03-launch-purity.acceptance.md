---
slug: cuda-md-review-cleanup-03-launch-purity
criteria:
  - id: ac-001
    summary: PairLJCut::compute() free of cudaStreamSynchronize and D2H
    type: code
    pass_when: |
      `grep -nE 'cudaStreamSynchronize|cudaMemcpy.*HostToDevice|cudaMemcpy.*DeviceToHost'
      src/cuda/potentials/pair_lj_cut.cu` returns 0 lines. The grid is sized
      via `pair_capacity` (constant at build) and the kernel's existing
      `if (idx >= *n_pairs) return;` early-exit handles the over-launch.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: NoseHooverNPT consolidates reductions into a single 8-double pinned buffer + event polling
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu` declares a kernel
      `compute_npt_reductions_k` that writes `{ke_total, ke_x, ke_y, ke_z,
      stress_xx, stress_yy, stress_zz, energy}` (or a documented superset);
      readback is via `cudaMemcpyAsync(pinned, dev, 8 * sizeof(double),
      cudaMemcpyDeviceToHost, defer_stream)` followed by `cudaEventQuery`
      polling (no `cudaStreamSynchronize` on the main stream); old per-component
      D2H copies at lines 250, 268, 399, 507-512, 542-545 are removed.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: NPT box-diagonal D2H gated by device version counter
    type: code
    pass_when: |
      `cuda::Box` (or `NoseHooverNPT::State`) carries a device int
      `box_version_counter_`; `propagate_x_and_box_` increments it after
      writing new diagonal values. Host side caches last-read counter and
      skips the D2H if unchanged. `git grep cudaMemcpy.*box`
      `src/cuda/md/nose_hoover_npt.cu` shows the copy guarded by a
      version-mismatch branch.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: NHC chain fused into single scale_vel_chain_k kernel
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu` `nhc_atoms_half_` and `nhc_baro_half_`
      launch a single `scale_vel_chain_k` kernel per call, passing the
      host-computed `eta_t_dot[mtchain]` array via kernel argument or
      `cudaMemcpyToSymbolAsync`. The previous mtchain-many interleaved
      launches are removed.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: Backend::step() free of per-phase cudaStreamSynchronize
    type: code
    pass_when: |
      `src/cuda/driver.cu:71-72` no longer contains
      `cudaStreamSynchronize(context_.stream)` or
      `cudaStreamSynchronize(context_.defer_stream)` inside the phase loop.
      `Backend::run()` may sync at the end of the run. `grep -nE
      'cudaStreamSynchronize' src/cuda/driver.cu` returns at most 1 hit
      (the end-of-run sync).
    status: pending
    last_checked: ""
  - id: ac-006
    summary: AsyncBridge uses cudaLaunchHostFunc instead of stream sync
    type: code
    pass_when: |
      `src/cuda/async_bridge.cu:39` region: `cudaStreamSynchronize(defer_stream)`
      before `on_data_ready` is replaced with `cudaLaunchHostFunc(defer_stream,
      &dispatch_callback, this)` (or equivalent). `grep -n cudaStreamSynchronize
      src/cuda/async_bridge.cu` returns 0 hits.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: ZeroMomentum runs without D2H+sync, OR documented integrity exception in NOTES
    type: code
    pass_when: |
      EITHER (a) `src/cuda/md/init.cu:233` no longer issues
      `cudaStreamSynchronize` or D2H; v_cm is computed entirely on device
      via a single-thread kernel reading device-resident `{Σmv, Σm}`;
      `remove_com` reads the device scratch directly. OR (b)
      `.claude/NOTES.md` gains an "integrity components" classification row
      explicitly listing `ZeroMomentum` as a permitted host-barrier exception
      with rationale (per `CLAUDE.md ## Architecture § CUDA execution model`
      classification table).
    status: pending
    last_checked: ""
  - id: ac-008
    summary: extended_hamiltonian_ computed on device (no per-step D2H/H2D)
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu:540-545, 598-602` removes the 2 D2H +
      2 H2D copies for `H`/`KE_global_`. The conserved Hamiltonian is
      computed by a device kernel; host readout happens only on demand
      (e.g. via the 8-double pinned buffer in ac-002).
    status: pending
    last_checked: ""
  - id: ac-009
    summary: Extended hot-path test passes with zero syncs over a 100-step NPT loop
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_hot_path_no_sync --output-on-failure`
      exits 0. The test drives a 100-step `cuda::md::MD` loop with
      `NoseHooverNPT` + `NeighborList` + `PairLJCut` + `HarmonicBond` +
      `HarmonicAngle` + `MultiHarmonicDihedral` + `VirialStress`, and asserts
      `cuda::detail::sync_count() == 0` at the end of every step.
    status: pending
    last_checked: ""
  - id: ac-010
    summary: Driver per-phase sync gate test passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_driver_no_per_phase_sync --output-on-failure`
      exits 0. A 50-step `cuda::md::MD` NVE loop with `NeighborList` +
      `PairLJCut` reports zero `cuda::stream_synchronize` calls during the
      hot path.
    status: pending
    last_checked: ""
  - id: ac-011
    summary: NPT correctness preserved (existing test passes)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_cuda_nose_hoover_npt --output-on-failure`
      exits 0. The conserved-Hamiltonian drift and T/P targets remain within
      the existing test's tolerances; the 1-step lag in pinned-buffer reads
      does not destabilize the NHC feedback loop.
    status: pending
    last_checked: ""
