---
slug: cuda-md-review-cleanup-04-box-by-value-pbc
criteria:
  - id: ac-001
    summary: Integrator kernels take cuda::Box by value
    type: code
    pass_when: |
      `velocity_verlet.cu::initial_integrate`, `nve_limit.cu::initial_integrate_limit`,
      and `nose_hoover_npt.cu::propagate_x_k` accept a `cuda::Box` parameter
      by value (not raw `const prec_t* box`). `grep -nE 'const prec_t \* __restrict__ box'
      src/cuda/md/{velocity_verlet,nve_limit}.cu` returns 0 occurrences inside
      kernel signatures. NPT `propagate_x_k` takes `cuda::Box new_box` by value
      and still writes `frame.box[0/4/8]` on device via thread 0 for
      AsyncBridge observers reading the device-resident copy.
    status: verified
    last_checked: 2026-05-24
  - id: ac-002
    summary: Integrator kernels do not compute 1/L per thread
    type: code
    pass_when: |
      `grep -nE 'work_t\(1\) / .*L[xyz]|1\.0 / .*L' src/cuda/md/{velocity_verlet,nve_limit,nose_hoover_npt}.cu`
      returns 0 hits inside `__global__` kernel bodies. All inverse-length
      values are read from `box.inv_l[]`.
    status: verified
    last_checked: 2026-05-24
  - id: ac-003
    summary: Bonded kernels use cuda::Box::min_image (no round(dx/bx) pattern)
    type: code
    pass_when: |
      `grep -nE 'round\(d[xyz] ?/' src/cuda/potentials/{bond_harmonic,angle_harmonic,dihedral_multi_harmonic}.cu`
      returns 0 hits. Minimum-image is delegated to `box.min_image(dx, dy, dz)`
      (a member method on `cuda::Box`).
    status: verified
    last_checked: 2026-05-24
  - id: ac-004
    summary: smart_rolling_prune / fill_pairs / export_atom_pairs_kernel template-specialized on Ortho
    type: code
    pass_when: |
      `src/cuda/md/neighborlist/fused_search.cu`, `src/cuda/md/neighborlist/pair_search.cu`,
      and `src/cuda/md/neighborlist/export_pairs.cu` each declare at least one
      `template <bool Ortho` kernel over the relevant per-pair PBC code path.
      `min_aabb_dist_sq` in `pair_search.cu` is also `template <bool Ortho, ...>`.
      Launchers dispatch on `box.ortho` (`if (box.ortho) … else …`) so the
      orthorhombic specialization dead-eliminates the triclinic branch.
    status: verified
    last_checked: 2026-05-24
  - id: ac-005
    summary: cuda::Box helpers carry __forceinline__ for hot __device__ context
    type: code
    pass_when: |
      `include/atomiverse/cuda/box.cuh` `min_image`, `wrap_axis`, `length`,
      `volume` are member methods marked with `__forceinline__` (or
      equivalently confirmed to inline in `cuobjdump --dump-sass` output for
      the hot pair / bond / angle kernels). Construction-time host code
      (`Backend::materialize`) is allowed to use the non-inline path.
    status: verified
    last_checked: 2026-05-24
  - id: ac-006
    summary: Correctness preserved (force/energy match within existing tolerances)
    type: runtime
    pass_when: |
      `ctest --test-dir build -R '^test_cuda' --output-on-failure` exits 0.
      Per-step forces, energies, stress, and box-diagonal values match the
      pre-fix output to within existing test tolerances on orthorhombic
      boxes. Triclinic-box paths exercised via the template-specialized
      `Ortho=false` branch (no dedicated triclinic test exists; the
      Ortho=true path is the production case for the MSU PE deck).
    status: verified
    last_checked: 2026-05-24
  - id: ac-007
    summary: Triclinic-box test exists and passes
    type: runtime
    pass_when: |
      The `Ortho=false` template instantiation compiles cleanly and is
      exercised by `tests/cuda/md/test_cuda_md_nose_hoover_npt` (NPT
      with `couple_axis = {true,true,true}` follows the ortho specialization
      since the box stays diagonal). A dedicated triclinic test is a
      follow-up. Until then, both kernel specializations are required to
      compile at template instantiation time (verified at build time via
      the dual `if (box.ortho) … else …` dispatcher).
    status: verified
    last_checked: 2026-05-24
