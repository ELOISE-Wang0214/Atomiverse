---
slug: cuda-md-review-cleanup-08-placement-and-conventions
criteria:
  - id: ac-001
    summary: src/cuda/md/neighborlist/ relocated to src/cuda/neighborlist/
    type: code
    pass_when: |
      `find src/cuda/md/neighborlist -type f` returns empty.
      `find src/cuda/neighborlist -type f` returns the previous file set
      (setup.cu, grid.cu, sort_and_cluster.cu, pair_search.cu,
      fused_search.cu, export_pairs.cu, plus the *.cuh internals). All
      includes updated; `cmake --build build` is clean.
    status: pending
    last_checked: ""
  - id: ac-002
    summary: cuda::VirialStress relocated to cuda/md/virial_stress.{cuh,cu}
    type: code
    pass_when: |
      `find include/atomiverse/cuda/observer -name stress.cuh` returns empty.
      `find src/cuda/observer -name stress.cu` returns empty.
      `include/atomiverse/cuda/md/virial_stress.cuh` and
      `src/cuda/md/virial_stress.cu` exist with the relocated content. All
      includes updated; `cmake --build build` is clean.
    status: pending
    last_checked: ""
  - id: ac-003
    summary: angle_harmonic.cuh docstring fixed to [rad]
    type: code
    pass_when: |
      `include/atomiverse/cuda/potentials/angle_harmonic.cuh` `@param theta0`
      line (around line 18) reads `[rad]` and does NOT contain the legacy
      `[degrees] (converted to radians internally)` text. The new wording
      should reference the API contract: `theta0_ already in radians`.
    status: pending
    last_checked: ""
  - id: ac-004
    summary: Backward-compat alias headers removed
    type: code
    pass_when: |
      `ls include/atomiverse/cuda/potentials/lj.cuh
       include/atomiverse/cuda/potentials/harmonic_bond.cuh
       include/atomiverse/cuda/potentials/harmonic_angle.cuh
       include/atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh`
      all return "No such file" errors. Repo-wide
      `grep -rn 'include.*cuda/potentials/\(lj\|harmonic_bond\|harmonic_angle\|multi_harmonic_dihedral\)\.cuh'`
      returns 0 hits.
    status: pending
    last_checked: ""
  - id: ac-005
    summary: Multi-GPU halo-exchange CUDA gap documented in NOTES
    type: code
    pass_when: |
      `.claude/NOTES.md` contains an entry titled (or referencing)
      "Multi-GPU halo exchange CUDA backend gap" that documents the
      divergence from CPU `#ifdef ATV_BUILD_MPI` auto-injection. Marked
      as "out of scope for current chain" with cross-ref to a future
      multi-GPU spec.
    status: pending
    last_checked: ""
  - id: ac-006
    summary: .claude/architecture.md refreshed (no stale SCFWorkspace claim)
    type: code
    pass_when: |
      `grep -n 'SCFWorkspace' .claude/architecture.md` does NOT match any
      line claiming `cuda::Frame` embeds it. The architecture blueprint
      is regenerated post-`/mol:map` and reflects the current
      `cuda::Frame` (post-`rust-frame-lowering-04-cuda-torch` state).
    status: pending
    last_checked: ""
  - id: ac-007
    summary: Full repo build + test pass after relocations
    type: runtime
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON` configures
      cleanly; `cmake --build build -j$(nproc)` builds with zero errors;
      `ctest --test-dir build --output-on-failure` exits 0.
    status: pending
    last_checked: ""
