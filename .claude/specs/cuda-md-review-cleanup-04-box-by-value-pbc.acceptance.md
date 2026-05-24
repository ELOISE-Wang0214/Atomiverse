---
slug: cuda-md-review-cleanup-04-box-by-value-pbc
criteria:
  - id: ac-001
    summary: Integrator kernels take cuda::Box by value
    type: code
    pass_when: |
      The kernel signatures of `velocity_verlet_initial_k` /
      `velocity_verlet_final_k`, `nve_limit_*_k`, and the NPT per-atom kernels
      in `nose_hoover_npt.cu` accept a `cuda::Box` parameter by value (not
      raw `Real* box`). `grep -nE 'Real.*box\[0\]|Real.*box\[4\]|Real.*box\[8\]'
      src/cuda/md/{velocity_verlet,nve_limit,nose_hoover_npt}.cu` returns 0
      occurrences inside the kernel bodies.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: Integrator kernels do not compute 1/L per thread
    type: code
    pass_when: |
      `grep -nE 'work_t\(1\) / .*box\[|1\.0 / .*L' src/cuda/md/{velocity_verlet,nve_limit,nose_hoover_npt}.cu`
      returns 0 hits inside `__global__` kernel bodies. All inverse-length
      values are read from `box.inv_l[]`.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: Bonded kernels use cuda::Box::min_image (no round(dx/bx) pattern)
    type: code
    pass_when: |
      `grep -nE 'round\(d[xyz] ?/' src/cuda/potentials/{bond_harmonic,angle_harmonic,dihedral_multi_harmonic}.cu`
      returns 0 hits. The minimum-image computation is delegated to
      `cuda::Box::min_image(dx, dy, dz)` (or equivalent ortho fast path).
    status: pending
    last_checked: ""
  - id: ac-004
    summary: smart_rolling_prune / export_atom_pairs_kernel / fill_pairs template-specialized on Ortho
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/fused_search.cu`, `src/cuda/neighbor_list.cu`,
      and `src/cuda/md/neighborlist/pair_search.cu` each contain at least one
      `template <bool Ortho>` declaration over the relevant kernels. The
      launcher selects the ortho-specialized instantiation when
      `cuda::Box::ortho` is true.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: cuda::Box::min_image and wrap_axis inline in __device__ context
    type: code
    pass_when: |
      `include/atomiverse/cuda/box.cuh` `min_image` and `wrap_axis` carry
      `__forceinline__` (or are confirmed by SASS inspection to inline at
      every call site). No `cuda::Box::min_image` call generates a function
      call in `cuobjdump --dump-sass` output for the hot pair / bond / angle
      kernels.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: Correctness preserved (force/energy match within existing tolerances)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials|cuda_neighbor'
      --output-on-failure` exits 0. Per-step forces and energies match the
      pre-fix values to within the existing test tolerances on both ortho
      and triclinic boxes.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: Triclinic-box test exists and passes
    type: runtime
    pass_when: |
      A test exists (existing or new — e.g. `test_cuda_neighborlist_triclinic`)
      that exercises the `Ortho=false` template path. `ctest --test-dir build
      -R triclinic --output-on-failure` exits 0.
    status: pending
    last_checked: ""
