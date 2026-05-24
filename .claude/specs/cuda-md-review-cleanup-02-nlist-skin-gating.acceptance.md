---
slug: cuda-md-review-cleanup-02-nlist-skin-gating
criteria:
  - id: ac-001
    summary: launch_cluster_displacement wired into NeighborList::launch
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu` `NeighborList::launch()` calls
      `launch_cluster_displacement` near the top; `rebuild_flag` (device int)
      is read via lagged `cudaMemcpyAsync` D2H + `cudaEventQuery` poll;
      the rest of the rebuild pipeline (Morton sort, CUB radix, AABB, grid,
      coarse count, fill, compact, csr_to_flat, prune, export) is wrapped in
      `if (rebuild_needed) { ... }`.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: ref_cluster_{x,y,z} populated after rebuild
    type: code
    pass_when: |
      After every full rebuild path in `NeighborList::launch()`, a new
      `update_ref_clusters` kernel (or inline equivalent) writes the current
      cluster centers into `State::buf.ref_cluster_x/y/z`. `git grep`
      confirms `ref_cluster_x` appears as a write target in
      `src/cuda/md/neighborlist/fused_search.cu` or `src/cuda/neighbor_list.cu`.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: Skin-gating runtime test passes (≤2 radix-sort invocations over 30 steps)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_neighborlist_skin_gating --output-on-failure`
      exits 0. The test runs 30 MD steps with atoms held stationary (zero force,
      or fixed-position dump) at skin=0.3·cutoff, instruments the CUB radix-sort
      kernel invocation count (via a sigil counter or stream-trace), and
      asserts `radix_sort_invocations <= 2`.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: launch_setup_cell_device not called per step for fixed box
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu:190` region: `launch_setup_cell_device<<<1,1>>>`
      is guarded by a host-side box-change check (e.g. `if
      (box_changed_since_last_call_) { ... }`). For NVE/NVT runs the kernel is
      not invoked after the first call.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: smart_rolling_prune and export_atom_pairs_kernel size grid by actual pair count, not capacity
    type: code
    pass_when: |
      Neither `src/cuda/md/neighborlist/fused_search.cu:233-235` nor
      `src/cuda/neighbor_list.cu:312` launches with `pair_capacity` as the
      grid-sizing argument. Instead they use a lagged-read of
      `flat_pair_counter` (or `n_pairs`) from `State`. `git grep pair_capacity`
      inside these two files returns no launch-config-related hits.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: setup_cell precompute moved to build() + __constant__ upload
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/setup.cu:178,193` `setup_cell_from_*<<<1,1>>>`
      single-thread kernels are not called per step. Host computation in
      `NeighborList::build()` produces stencil + dims, and
      `cudaMemcpyToSymbolAsync` uploads them to a `__constant__` symbol.
      Re-upload occurs only on box change.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: cut_skin_sq and cutoff_sq cached in State
    type: code
    pass_when: |
      `NeighborList::State` (`src/cuda/md/neighborlist/nlist_impl.cuh`) has
      fields `Real cut_skin_sq_` and `Real cutoff_sq_` populated in `build()`
      and read in `launch()`. `git grep 'r_search * r_search\|cutoff_ \* cutoff_'`
      in `src/cuda/neighbor_list.cu` returns 0 hits.
    status: pending
    last_checked: ""
  - id: ac-008
    summary: CUB temp storage pre-allocated at max in build()
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/sort_and_cluster.cu:112-120` no longer issues
      `cudaMalloc` via `temp.resize()` after `build()` completes. Build queries
      CUB API for the max temp size and resizes once. Verified by inserting a
      runtime assertion that `temp.size()` does not change between `launch()`
      calls.
    status: pending
    last_checked: ""
  - id: ac-009
    summary: interaction_mask memset moved off per-step path
    type: code
    pass_when: |
      `src/cuda/neighbor_list.cu:267-269` `cudaMemsetAsync(interaction_mask,
      0xFF, ms2 * sizeof(uint64_t))` is either (a) moved to `build()`
      (one-time) or (b) sized to actually-written portion only. `git grep
      'cudaMemsetAsync.*interaction_mask'` in `src/cuda/neighbor_list.cu`
      shows the call no longer uses `pair_capacity * 8` as the byte count for
      every step.
    status: pending
    last_checked: ""
  - id: ac-010
    summary: Existing MD correctness tests pass (no force/energy drift)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|test_md\.cu$' --output-on-failure`
      exits 0. The skin-gating logic preserves force and energy values to
      within existing tolerances.
    status: pending
    last_checked: ""
