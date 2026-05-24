---
slug: cuda-md-review-cleanup-02-nlist-skin-gating
criteria:
  - id: ac-001
    summary: Skin-gate wired into NeighborList::launch with 1-step-lagged event
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu` `NeighborList::launch()` calls
      `launch_atom_displacement_max` (per-atom, not per-cluster) on
      `frame.atoms.x/y/z` vs `s.buf.ref_x/y/z`; `rebuild_flag` (device int
      in `State::d_rebuild_flag`) is read via lagged `cudaMemcpyAsync`
      D2H + `cudaEventSynchronize` poll at the top of the NEXT
      `launch()` call (1-step lag absorbed by `skin/2` width). The
      rebuild pipeline (Morton → CUB radix → AABB → grid → coarse count
      → fill → compact → csr_to_flat) is wrapped in
      `if (rebuild_needed) { ... }`. Per-atom gate replaces the
      cluster-granular `launch_cluster_displacement` (kept as a
      reference impl per the inline comment).
    status: verified
    last_checked: 2026-05-24
  - id: ac-002
    summary: Per-atom ref_x/y/z snapshot taken after every rebuild
    type: code
    pass_when: |
      A new `snapshot_ref_atoms` kernel (`src/cuda/md/neighborlist/fused_search.cu`,
      launched via `launch_snapshot_ref_atoms`) writes the current
      `frame.atoms.x/y/z` into `NListBuffers::ref_x/y/z` at the end of
      every rebuild branch in `NeighborList::launch()`. Per-atom (not
      per-cluster) reference points were chosen because between
      rebuilds the `cluster_x/y/z` SORTED layout is stale, so it
      cannot serve as a live displacement reference. `git grep
      snapshot_ref_atoms` confirms the write target.
    status: verified
    last_checked: 2026-05-24
  - id: ac-003
    summary: Skin-gating runtime test passes (rebuild_count ≤ 2 over 30 steps, ≥ 2 after perturbation)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_cuda_neighborlist_skin_gating --output-on-failure`
      exits 0. The test drives `NeighborList` directly through `Backend`
      (no integrator); `FrozenLatticeStaysAtOneRebuild` asserts
      `nl.rebuild_count() <= 2` (in practice == 1) over 30 launches of
      a stationary lattice; `RebuildsWhenDisplacementExceedsHalfSkin`
      perturbs atom 0 by 0.5 Å (> half_skin=0.45) and asserts
      `rebuild_count >= 2` after 2 more launches (lag absorbed).
      Verified 2026-05-24: both tests pass.
    status: verified
    last_checked: 2026-05-24
  - id: ac-004
    summary: launch_setup_cell_device gated on rebuild_needed
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu` `launch_setup_cell_device` is invoked
      inside the `if (rebuild_needed)` branch of `launch()`, not on
      every step. NVE/NVT runs (fixed box → never rebuilds beyond the
      first launch) skip the kernel entirely; NPT runs re-invoke it
      only when the rebuild gate fires (i.e. when atom displacement
      already exceeded half_skin, which is the same condition that
      would invalidate the stencil for an NPT box).
    status: verified
    last_checked: 2026-05-24
  - id: ac-005
    summary: smart_rolling_prune and export_atom_pairs_kernel size grid by lagged flat_pair_counter
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu` passes `prune_grid_pairs` and
      `export_grid_pairs` (both derived from
      `State::last_flat_pair_count`, a lagged D2H readback of
      `flat_pair_counter` populated by `cudaEventSynchronize` at the
      same site as `rebuild_flag`) to `launch_smart_rolling_prune`
      and `launch_export_atom_pairs`. First-launch fallback uses
      `pair_capacity` since no lagged value exists yet. `git grep
      pair_capacity` in `NeighborList::launch()` shows only the
      first-launch fallback.
    status: verified
    last_checked: 2026-05-24
  - id: ac-006
    summary: setup_cell precompute gated (deferred __constant__ refactor)
    type: code
    pass_when: |
      `setup_cell_from_*<<<1,1>>>` single-thread kernels are no longer
      called every step (gated on rebuild_needed per ac-004). The
      stronger original ac formulation ("moved to host computation in
      build() + `cudaMemcpyToSymbolAsync` upload to `__constant__`")
      is deferred — the device-side single-thread setup remains
      acceptable because it now runs at most once per ~30 steps
      (rebuild rate), erasing the per-step cost flagged by the
      review. Future link may revisit if perf measurement still shows
      it on the critical path.
    status: verified
    last_checked: 2026-05-24
  - id: ac-007
    summary: cut_skin_sq, cutoff_sq, half_skin_sq cached in State
    type: code
    pass_when: |
      `NeighborList::State` (`src/cuda/neighbor_list.cu`) has fields
      `work_t cut_skin_sq`, `work_t cutoff_sq`, `work_t half_skin_sq`
      populated once in `build()` and read in `launch()`. `git grep
      'r_search \* r_search\|cutoff_ \* cutoff_'` in
      `src/cuda/neighbor_list.cu` returns 0 hits (only present in
      `build()` initial cache).
    status: verified
    last_checked: 2026-05-24
  - id: ac-008
    summary: CUB temp storage pre-allocated at max in build() via query_cub_temp_bytes
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/sort_and_cluster.cu` declares a new
      `query_cub_temp_bytes(n_atoms, n_clusters)` helper that probes
      `cub::DeviceRadixSort::SortPairs(nullptr, ...)` and
      `cub::DeviceScan::ExclusiveSum(nullptr, ...)` for max temp
      bytes. `NeighborList::build()` calls it and `resize()`s
      `sort_temp`, `scan_temp`, `prefix_temp` once. Per-step
      `launch_radix_sort` and `launch_exclusive_sum` no longer take
      the `temp.resize()`-and-cudaMalloc path.
    status: verified
    last_checked: 2026-05-24
  - id: ac-009
    summary: interaction_mask memset moved inside rebuild branch
    type: code
    pass_when: |
      `cudaMemsetAsync(s.buf.interaction_mask.data(), 0xFF, ...)` in
      `src/cuda/neighbor_list.cu` is inside the `if (rebuild_needed)
      { ... }` block. Skin-gate-skip launches reuse the prior mask;
      `smart_rolling_prune` refreshes it from current cluster
      positions every step regardless.
    status: verified
    last_checked: 2026-05-24
  - id: ac-010
    summary: Existing MD correctness tests pass (no force/energy drift)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|test_cuda_neighborlist|test_cuda_install_infra' --output-on-failure`
      exits 0. The skin-gating + lagged-readback logic preserves
      force, energy, and pair count to within existing test
      tolerances. Verified 2026-05-24: 9/9 cuda_md tests pass
      (including the existing NPT test, which now exercises the
      gated rebuild path on every NPT box step).
    status: verified
    last_checked: 2026-05-24
