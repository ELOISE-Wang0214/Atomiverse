---
slug: cuda-md-review-cleanup-05-reductions-and-atomics
criteria:
  - id: ac-001
    summary: cuda::block_reduce_sum + block_reduce_sum_n in launch.cuh
    type: code
    pass_when: |
      `include/atomiverse/cuda/launch.cuh` declares
      `template <typename T> __device__ __forceinline__ T
      block_reduce_sum(T value, T* warp_scratch)` and
      `template <int N, typename T> __device__ __forceinline__ void
      block_reduce_sum_n(const T (&values)[N], T* warp_scratch, T (&out)[N])`
      backed by `cooperative_groups::reduce` (compiles to
      `__reduce_add_sync` on sm_80+, shfl_xor fallback otherwise).
      Also exposes `cuda::kWarpSize = 32` and
      `cuda::warps_per_block(int block = kDefaultBlock)` so launch
      sites size the per-warp shared scratch without magic numbers.
      Initially landed as a separate `cg_reduce.cuh`; merged into
      `launch.cuh` per user feedback to avoid `cg` / coarse-grain
      naming ambiguity and to keep kernel-author helpers centralised.
    status: verified
    last_checked: 2026-05-24
  - id: ac-002
    summary: Hot force kernels use warp-shuffle reduction (not 6 sequential shared-mem)
    type: code
    pass_when: |
      `src/cuda/potentials/pair_lj_cut.cu`, `bond_harmonic.cu`,
      `angle_harmonic.cu`, `dihedral_multi_harmonic.cu` each replace
      the legacy 7-pass `extern __shared__ work_t sdata[]` tree
      reduction with a single `cuda::block_reduce_sum_n<7>` call
      packing energy + 6 stress components. `grep -c '__syncthreads'`
      drops by ≥ 4 lines per file vs the pre-link-05 baseline. Shared
      memory per block drops from `blockDim.x * sizeof(work_t)` to
      `warps_per_block * sizeof(work_t)` (~8× less).
    status: verified
    last_checked: 2026-05-24
  - id: ac-003
    summary: KineticEnergy and kinetic_stress_k fused into single 7-value kernel
    type: code
    pass_when: |
      `src/cuda/md/kinetic_energy.cu` declares
      `kinetic_energy_stress` (renamed from `kinetic_energy`) that
      atomic-adds 7 partials per block: 1 ke scalar +
      6 stress_kinetic components (with kMvv2e folded in).
      `src/cuda/observer/stress.cu` no longer contains a
      `kinetic_stress_k` kernel definition or launch; its `launch()`
      is narrowed to a single `cudaMemsetAsync` clearing
      `stress_kinetic[6]` (so KineticEnergy's atomic-add lands in a
      known-zero buffer). VirialStress::State no longer allocates
      `d_stress_tmp`. The previous separate launch
      (`kinetic_stress_k`) is eliminated.
    status: verified
    last_checked: 2026-05-24
  - id: ac-004
    summary: kinetic_stress_k shared usage reduced to <1 KiB (warp-shuffle path)
    type: code
    pass_when: |
      The fused `kinetic_energy_stress` kernel declares `__shared__
      work_t warp_scratch[]` sized to `warps_per_block(block) *
      sizeof(work_t)` = ≤ 64 bytes for the default 256-thread block.
      The previous `6 * blockDim * sizeof(work_t)` ≈ 12 KiB
      allocation pattern is gone.
    status: verified
    last_checked: 2026-05-24
  - id: ac-005
    summary: reduce_ke_axes_k uses warp-shuffle (not 6 KiB shared)
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu` `reduce_ke_axes_k` uses
      `cuda::block_reduce_sum_n<3>` for the 3-axis KE reduction.
      Shared memory sized via `warps_per_block(block) * sizeof(double)`
      = ≤ 64 bytes (vs the pre-link-05 `3 * blockDim * sizeof(double)`
      ≈ 6 KiB).
    status: verified
    last_checked: 2026-05-24
  - id: ac-006
    summary: 'if (s_energy[0] != 0) atomicAdd' pattern removed
    type: code
    pass_when: |
      `grep -rn 'if.*s_energy\[0\].*!= .*atomicAdd' src/cuda/potentials/
      src/cuda/observer/ src/cuda/md/` returns 0 hits. The new
      warp-shuffle reductions emit unconditional atomicAdds —
      `atomicAdd(_, 0)` is free on modern GPUs and the conditional
      branch wasn't pulling its weight.
    status: verified
    last_checked: 2026-05-24
  - id: ac-007
    summary: LJ pair uses cluster-staged force write (≤2 atomicAdds per atom per block)
    type: code
    pass_when: |
      DEFERRED to a follow-up link. The cluster-staged restructure
      requires reorganising `pair_lj_cut` to consume cluster-pair
      input (flat_ci/flat_cj groupings) and stage partial forces in
      shared memory before a single atomicAdd per atom per block.
      This is a non-mechanical data-structure rewrite with
      regression risk; user-approved deferral to keep link 05
      shipping the 8 reduction migrations + KE/stress fusion that
      did land. Backlog: track under
      `cuda-md-review-cleanup-09-bench-and-regression` umbrella or
      a new follow-up link.
    status: deferred
    last_checked: 2026-05-24
  - id: ac-008
    summary: Bonded kernels reduce per-topology atomicAdd count
    type: code
    pass_when: |
      DEFERRED to a follow-up link. The per-atom owner-list reverse
      index requires building a new topology data structure (per-atom
      owner lists in BondHarmonic/AngleHarmonic/MultiHarmonicDihedral
      build()) and rewriting the kernels to iterate by destination
      atom instead of by topology entry. Non-mechanical change;
      user-approved deferral. Backlog as ac-007.
    status: deferred
    last_checked: 2026-05-24
  - id: ac-009
    summary: Correctness preserved (force/energy/stress match within tolerances)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials|cuda_observer'
      --output-on-failure` exits 0. Per-step force, energy, stress,
      and kinetic-energy values match pre-fix to within existing test
      tolerances on both sparse (random-LJ) and dense (PE topology)
      systems. Verified 2026-05-24: 23/23 cuda tests pass including
      existing NPT, install-infra-capability, neighborlist-skin-gating,
      pair-LJ-special-bonds, potentials-bonded, hot-path-no-sync.
    status: verified
    last_checked: 2026-05-24
