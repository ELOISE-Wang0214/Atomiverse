---
slug: cuda-md-review-cleanup-07-kernel-header-and-microfixes
criteria:
  - id: ac-001
    summary: lj_kernel.cuh header matches pair_lj_cut.cu definition (single source of truth)
    type: code
    pass_when: |
      `include/atomiverse/cuda/potentials/lj_kernel.cuh` declares
      `cuda::kernels::lj_force` with the SAME argument count and types as
      the definition at `src/cuda/potentials/pair_lj_cut.cu`. The torch
      backend at `src/torch/potentials/lj.cu:96-105` invokes the kernel with
      the matching signature and compiles via
      `cmake --build build -DATV_BUILD_TORCH=ON` (when torch is enabled).
    status: pending
    last_checked: ""
  - id: ac-002
    summary: harmonic_bond_kernel.cuh signature aligned similarly
    type: code
    pass_when: |
      `include/atomiverse/cuda/potentials/harmonic_bond_kernel.cuh` declares
      `cuda::kernels::harmonic_bond_force` with the same argument count and
      types as the definition in `src/cuda/potentials/bond_harmonic.cu`.
      Torch backend call site matches.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: box_muller uses sincos (single PTX op)
    type: code
    pass_when: |
      `grep -nE 'sincos|sincospi' src/cuda/md/fix_langevin.cu` returns at
      least 1 hit. `grep -nE 'cos\(2.*PI|sin\(2.*PI' src/cuda/md/fix_langevin.cu`
      returns 0 hits as separate calls.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: berendsen_scale_k receives μ as kernel arg (no per-block cbrt)
    type: code
    pass_when: |
      `grep -n cbrt src/cuda/md/barostat.cu` returns 0 hits inside any
      `__global__` kernel body. `cbrt` may appear in host `launch()` code
      computing μ before kernel invocation.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: dihedral dv_dcos coefficients precomputed in build()
    type: code
    pass_when: |
      `cuda::MultiHarmonicDihedral::State` (in
      `include/atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh`)
      holds a device array `dv_dcos_coeffs_` populated in `build()` with
      `{C[0], C[1], 2*C[2], 3*C[3], 4*C[4]}`. The dihedral kernel reads from
      this array instead of computing `2*Cn[2]` etc. inline.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: inv_cell cached and re-uploaded only on box change
    type: code
    pass_when: |
      `cuda::NeighborList::State` carries an `inv_cell_cache_[9]`; the host
      side re-computes only when box volume / shape changes and uploads via
      `cudaMemcpyToSymbolAsync` to `__constant__`. The per-step
      `setup_cell` kernel (already moved off the hot path in link 02) no
      longer recomputes inv_cell device-side.
    status: pending
    last_checked: ""
  - id: ac-007
    summary: deform.cu box-diagonal write factored into separate <<<1,1>>> kernel
    type: code
    pass_when: |
      `src/cuda/md/deform.cu:113` no longer contains `if (i==0) { box[0] =
      ...; }` inside the per-atom propagate kernel. A separate
      `update_box_diagonal_k<<<1,1>>>` (or equivalent) kernel handles the
      box-diagonal write, scheduled after `propagate_x_k` on the same stream.
    status: pending
    last_checked: ""
  - id: ac-008
    summary: All existing CUDA MD + torch tests pass
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_potentials|cuda_fix_langevin'
      --output-on-failure` exits 0. If `ATV_BUILD_TORCH=ON`, torch backend
      tests also pass. Microfix correctness preserved.
    status: pending
    last_checked: ""
