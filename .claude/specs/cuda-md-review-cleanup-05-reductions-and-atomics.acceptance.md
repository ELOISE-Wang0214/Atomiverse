---
slug: cuda-md-review-cleanup-05-reductions-and-atomics
criteria:
  - id: ac-001
    summary: cg_reduce.cuh provides warp-shuffle reduction helper
    type: code
    pass_when: |
      `include/atomiverse/cuda/cg_reduce.cuh` declares a templated function
      that performs warp-level sum reduction over N values using
      `cooperative_groups::reduce` or `__reduce_add_sync(0xffffffff, val)`.
      Compiles standalone via `cmake --build build --target atomiverse_cuda`.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: Hot force kernels use warp-shuffle reduction (not 6 sequential shared-mem)
    type: code
    pass_when: |
      `src/cuda/potentials/pair_lj_cut.cu`, `bond_harmonic.cu`,
      `angle_harmonic.cu`, `dihedral_multi_harmonic.cu` each invoke either
      `__reduce_add_sync` or a function from `cuda/cg_reduce.cuh` in their
      stress-reduction block. The previous pattern of 6 `__syncthreads()`
      calls in a row is removed: `grep -c '__syncthreads()' <each file>` is
      reduced by at least 4 lines per file compared to baseline.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: KineticEnergy and kinetic_stress_k fused into single 7-value kernel
    type: code
    pass_when: |
      A single device kernel (call it `kinetic_energy_and_stress_k` or
      equivalent) in `src/cuda/md/kinetic_energy.cu` (or `observer/stress.cu`)
      writes 7 reduction values `{ke, stress_xx, stress_yy, stress_zz,
      stress_xy, stress_xz, stress_yz}`. The previous separate launch of
      `kinetic_stress_k` from `cuda::VirialStress::launch()` is removed.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: kinetic_stress_k shared usage reduced to <1 KiB (warp-shuffle path)
    type: code
    pass_when: |
      `src/cuda/observer/stress.cu` `kinetic_stress_k` (or the new fused
      kernel) declares `__shared__` storage totaling < 1024 bytes for the
      reduction. The previous `6 * blockDim * sizeof(work_t)` allocation
      pattern is gone.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: reduce_ke_axes_k uses warp-shuffle (not 6 KiB shared)
    type: code
    pass_when: |
      `src/cuda/md/nose_hoover_npt.cu:129-167` `reduce_ke_axes_k` uses
      warp-shuffle (or cooperative_groups warp reduce) for the 3-axis
      reduction. Shared memory usage < 1 KiB.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: 'if (s_energy[0] != 0) atomicAdd' pattern removed
    type: code
    pass_when: |
      `grep -rn 'if.*s_energy\[0\].*!= .*atomicAdd' src/cuda/potentials/ src/cuda/observer/ src/cuda/md/`
      returns 0 hits.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: LJ pair uses cluster-staged force write (≤2 atomicAdds per atom per block)
    type: code
    pass_when: |
      `src/cuda/potentials/pair_lj_cut.cu` force kernel declares
      `__shared__ work_t fi_local[CLUSTER_SIZE][3]` (or similar staging
      buffer); the per-pair `atomicAdd` count inside the inner pair loop
      is 0; final `atomicAdd` to global `force` happens once per atom per
      block. `grep -c atomicAdd src/cuda/potentials/pair_lj_cut.cu` shows
      the count dropped from 6 (pre-fix) to ≤ 2 inside the kernel body.
    status: pending
    last_checked: ""
  - id: ac-008
    summary: Bonded kernels reduce per-topology atomicAdd count
    type: code
    pass_when: |
      EITHER `bond_harmonic.cu`, `angle_harmonic.cu` introduce a per-atom
      owner-list reverse index and accumulate in register (with at most
      one atomicAdd per dest atom per topology batch), OR they use
      shared-mem partial accumulation per block. `grep -c atomicAdd`
      shows the per-topology atomic count reduced by at least 3× from
      baseline.
    status: pending
    last_checked: ""
  - id: ac-009
    summary: Correctness preserved (force/energy/stress match within tolerances)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials|cuda_observer'
      --output-on-failure` exits 0. Per-step force, energy, and stress
      values match pre-fix to within existing test tolerances on both
      sparse (random-LJ) and dense (PE topology) systems.
    status: pending
    last_checked: ""
