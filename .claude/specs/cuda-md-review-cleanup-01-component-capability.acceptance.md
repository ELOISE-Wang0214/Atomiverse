---
slug: cuda-md-review-cleanup-01-component-capability
criteria:
  - id: ac-001
    summary: cuda::ComponentBase declares reads(), needs_neighbor_list(), provides_neighbor_list() virtuals
    type: code
    pass_when: |
      `include/atomiverse/cuda/component.cuh` declares
      `virtual Capability reads() const`, `virtual bool needs_neighbor_list() const`,
      and `virtual bool provides_neighbor_list() const` on `cuda::ComponentBase`,
      with default returns `Capability::None`, `false`, `false`. Project compiles
      via `cmake --build build --target atomiverse` with zero warnings related
      to virtual override misuse.
    status: verified
    last_checked: 2026-05-24
  - id: ac-002
    summary: cuda::VirialStress declares Observer + Stress
    type: code
    pass_when: |
      `include/atomiverse/cuda/observer/stress.cuh` `cuda::VirialStress` overrides
      `role()` returning `ComponentRole::Observer` and `capabilities()` returning
      a value with `Capability::Stress` bit set.
    status: verified
    last_checked: 2026-05-24
  - id: ac-003
    summary: cuda::Angle and cuda::Dihedral declare Stress capability
    type: code
    pass_when: |
      `cuda::Angle::capabilities()` (in `include/atomiverse/cuda/potentials/angle.cuh`)
      and `cuda::Dihedral::capabilities()` (in `include/atomiverse/cuda/potentials/dihedral.cuh`)
      both return values with `Capability::Stress` bit set. Cluster 1 also
      adds `Stress` to `cuda::md::Barostat::reads()` for symmetry with
      `cpu::md::Barostat::reads()`.
    status: verified
    last_checked: 2026-05-24
  - id: ac-004
    summary: cuda::NoseHooverNPT declares reads()==Stress
    type: code
    pass_when: |
      `cuda::NoseHooverNPT::reads()` returns a value with `Capability::Stress` bit set
      (`include/atomiverse/cuda/md/nose_hoover_npt.cuh`). `cuda::md::Barostat`
      receives the same declaration.
    status: verified
    last_checked: 2026-05-24
  - id: ac-005
    summary: cuda::Pair (or all pair subclasses) declare needs_neighbor_list()==true
    type: code
    pass_when: |
      `cuda::PairLJCut::needs_neighbor_list()` returns `true` via the
      `cuda::Pair` base override (`final` — subclasses inherit).
    status: verified
    last_checked: 2026-05-24
  - id: ac-006
    summary: install_infra_ no longer uses name-string proxies or dynamic_cast<cuda::Pair*>
    type: code
    pass_when: |
      `grep -rn 'name() == "virial_stress"\|name() == "neighbor_list"\|dynamic_cast<const cuda::Pair\*>' include/atomiverse/cuda/md/`
      returns 0 lines (verified: exit 1). `cuda::md::MD::install_infra_`
      consults `reads()`, `needs_neighbor_list()`, `capabilities()`, and
      `provides_neighbor_list()` instead. The symmetric
      `provides_neighbor_list()` virtual was added to `ComponentBase` +
      `NeighborList` to retire the last `name() == "neighbor_list"` proxy.
    status: verified
    last_checked: 2026-05-24
  - id: ac-007
    summary: AsyncBridge-derived recorders/observers declare role()
    type: code
    pass_when: |
      `cuda::Thermo` (Recorder), `cuda::Recorder` (Recorder),
      `cuda::DensityRecorder` (Recorder), `cuda::RDF` (Observer),
      `cuda::MullikenObs` (Observer) each override `role()` (matching
      their CPU sibling). None inherit the default `None`.
    status: verified
    last_checked: 2026-05-24
  - id: ac-008
    summary: New capability-based auto-injection test passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_install_infra_capability --output-on-failure`
      exits 0. The test builds a `cuda::md::MD` with `PairLJCut` +
      `NoseHooverNPT`, calls `build()`, and asserts the resulting Plan
      contains: (a) a `VirialStress` component triggered by
      `NoseHooverNPT::reads()==Stress`, (b) end-to-end NPT run completes
      via `PairLJCut::needs_neighbor_list()` triggering the
      missing-NeighborList throw path. Verified 2026-05-24: 6 tests
      passed (smoke virtuals + end-to-end injection + missing-NL throw).
    status: verified
    last_checked: 2026-05-24
  - id: ac-009
    summary: Existing CUDA MD regression suite passes
    type: runtime
    pass_when: |
      `ctest --test-dir build -R 'cuda_md|cuda_nose_hoover|cuda_md_components' --output-on-failure`
      exits 0. No behavior change is expected; this gate guards against
      accidental regressions in the install_infra_ rewrite. Verified
      2026-05-24: 8/8 cuda_md tests passed; broader 19/19 cuda tests
      passed; full ctest suite pending background completion.
    status: verified
    last_checked: 2026-05-24
