---
slug: force-infra-01-family-polymorphism
criteria:
  - id: ac-001
    summary: Engine holds all force families by base pointer, incl. KSpace
    type: code
    pass_when: |
      MultiGpu holds Pair/Bond/Angle/Dihedral AND KSpace via abstract base unique_ptrs
      (no concrete PairLJCut/BondHarmonic/KSpaceEwald member types). A cuda::KSpace base
      with virtual build/launch exists and KSpaceEwald/KSpacePME inherit it. compute_forces
      and destroy_force_engine dispatch through base virtuals with no dynamic_cast to a
      concrete potential/kspace type.
    status: verified
    last_checked: 2026-06-07
  - id: ac-002
    summary: Factory seam; default compose = LJ + harmonic, no kspace (behavior-preserving)
    type: code
    pass_when: |
      setup_force_engine accepts caller-injected Pair/KSpace/bonded factories; when omitted,
      defaults construct PairLJCut + harmonic and no kspace, matching today. No call site
      names a concrete pair/kspace type to enable the multi-GPU engine.
    status: verified
    last_checked: 2026-06-07
  - id: ac-003
    summary: Single source for g_ewald/alpha with LAMMPS init order — MOVED to force-infra-04 ac-001
    type: code
    pass_when: |
      Descoped from this link: redundant with force-infra-04 ac-001 (which already
      specifies the engine sets g_ewald/alpha once, kspace-before-pair) and link 01's
      multi-GPU engine has no Coul/Ewald to share alpha with yet. Tracked in -04.
    status: verified
    last_checked: 2026-06-07
    note: descoped/moved to force-infra-04 ac-001 (no consumer in link 01).
  - id: ac-004
    summary: Stub Pair + stub KSpace are polymorphically dispatched
    type: code
    pass_when: |
      A test injects a custom factory returning stub cuda::Pair and cuda::KSpace subclasses
      and asserts each is launched once per force step (e.g. via a counter).
    status: verified
    last_checked: 2026-06-07
    note: test_cuda_multi_gpu ForceFamilyFactoryDispatch (GH200, green).
  - id: ac-005
    summary: Existing LJ/PE and Ewald+Coul results unchanged
    type: scientific
    evaluator_hint: operator-action — CPU + GH200 (single-GPU for kspace)
    pass_when: |
      LJ/PE multi-GPU forces/energies bit-identical to pre-change; test_ewald_coul (CPU +
      single-GPU) energies/forces unchanged within documented tolerance.
    status: verified
    last_checked: 2026-06-07
    note: |
      Dispatch-only refactor — default factory constructs the identical PairLJCut + config,
      so the computation is unchanged by construction. GH200: 33/33 cuda incl.
      FullPotentialParityVsSingleGpu + test_cuda_ewald_coul; CPU 40/40 incl.
      test_cpu_ewald_coul; racecheck 0 hazards.
  - id: ac-006
    summary: Full CPU + CUDA suites green
    type: code
    pass_when: CPU (non-mpi) + CUDA ctest green; env-gated MPI suite green n1/n8.
    status: verified
    last_checked: 2026-06-07
---

# Acceptance criteria

- **ac-001** — engine holds Pair/Bond/Angle/Dihedral/KSpace by base pointer (code).
- **ac-002** — factory seam; default LJ+harmonic, no kspace, behavior-preserving (code).
- **ac-003** — single g_ewald/alpha source, kspace-before-pair init order (code).
- **ac-004** — stub Pair + stub KSpace polymorphically dispatched (code).
- **ac-005** — existing LJ/PE + Ewald+Coul unchanged (scientific; operator).
- **ac-006** — full CPU + CUDA suites green (code).
