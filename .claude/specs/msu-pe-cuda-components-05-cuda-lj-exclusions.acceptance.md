---
slug: msu-pe-cuda-components-05-cuda-lj-exclusions
criteria:
  - id: ac-001
    summary: cuda::PairLJCut exposes set_special_bonds(s12,s13,s14) returning self
    type: code
    pass_when: |
      Compiling
      `cuda::PairLJCut(1.0, 3.0, 10.0).set_special_bonds(0.0, 0.0, 0.0)`
      succeeds and returns a reference/value of type cuda::PairLJCut;
      include/atomiverse/cuda/potentials/pair_lj_cut.cuh declares the
      member with the signature `PairLJCut& set_special_bonds(double,
      double, double)`.
    status: pending
  - id: ac-002
    summary: set_special_bonds rejects mixed-scaling combinations
    type: runtime
    pass_when: |
      A gtest case calling `set_special_bonds(0.0, 1.0, 0.0)` causes the
      call to throw `std::invalid_argument`; `set_special_bonds(0.5,
      0.5, 0.5)` and `set_special_bonds(1.0, 1.0, 1.0)` and
      `set_special_bonds(0.0, 0.0, 0.0)` do NOT throw.
    status: pending
  - id: ac-003
    summary: 4-atom chain with (0,0,0) yields zero CUDA LJ energy
    type: scientific
    pass_when: |
      Test fixture in tests/cuda/potentials/test_pair_lj_special_bonds.cu
      builds a 4-atom collinear chain with bonds (0-1,1-2,2-3), angles
      and the (0-1-2-3) dihedral; after one compute() with
      set_special_bonds(0,0,0) the device-side
      frame.scalars->energy == 0.0 exactly (or |E| < 1e-12 eV after D2H
      copy).
    status: pending
  - id: ac-004
    summary: CUDA vs CPU energy parity on 4-atom chain across (0,0,0) and (1,1,1)
    type: scientific
    pass_when: |
      For the same 4-atom chain geometry, running cpu::PairLJCut and
      cuda::PairLJCut with identical (epsilon, sigma, cutoff) and the
      same set_special_bonds value, |E_cuda - E_cpu| < 1e-5 eV
      (absolute) AND |E_cuda - E_cpu| / max(|E_cpu|, 1e-10) < 1e-4
      (relative) for both (0.0,0.0,0.0) and (1.0,1.0,1.0) configurations.
    status: pending
  - id: ac-005
    summary: empty exclusion list takes nullptr fast path with no regression
    type: scientific
    pass_when: |
      A 256-atom no-bond LJ test case running cuda::PairLJCut WITHOUT
      calling set_special_bonds produces total energy and per-atom forces
      bit-identical to a reference run captured before this change is
      merged (recorded as tests/cuda/potentials/refdata/lj_256_nobond.npz
      or equivalent), confirmed by gtest EXPECT_FLOAT_EQ on energy and
      max |Δf| == 0.
    status: pending
  - id: ac-006
    summary: build/destroy lifecycle frees device CSR without leak
    type: runtime
    pass_when: |
      Running the test under `compute-sanitizer --tool memcheck` on a
      `build → compute → destroy → build → compute → destroy` cycle
      reports 0 leaks attributable to the d_excl_offsets_ / d_excl_partners_
      allocations; after destroy(), both pointers equal nullptr (gtest
      EXPECT_EQ).
    status: pending
  - id: ac-007
    summary: lj_force kernel signature accepts exclusion CSR
    type: code
    pass_when: |
      include/atomiverse/cuda/potentials/lj_kernel.cuh declares
      `__global__ void lj_force(...)` with parameters including
      `const int* __restrict__ excl_offsets, const int* __restrict__
      excl_partners, cuda::work_t s_excl` (or 3 separate s12/s13/s14
      args); the kernel compiles cleanly with these added args.
    status: pending
  - id: ac-008
    summary: docstrings document units and supported s-combinations
    type: docs
    pass_when: |
      include/atomiverse/cuda/potentials/pair_lj_cut.cuh contains
      doxygen on set_special_bonds that states (a) s12/s13/s14 are
      dimensionless, (b) only all-equal or (0,0,0) combinations are
      supported in this revision, (c) the partner list is built at
      build() time from frame.bonds/angles/dihedrals.
    status: pending
  - id: ac-009
    summary: full check + test suite passes on gpugh node
    type: runtime
    pass_when: |
      `cmake --build build --target check` plus
      `ctest --test-dir build -L cuda --output-on-failure` exit 0 on the
      Dardel gpugh partition (GH200); test_pair_lj_special_bonds appears
      in the passing test list with no SKIP.
    status: pending
---

# Acceptance criteria

ac-001 — ac-002 lock the public API surface and input validation contract.

ac-003 is the headline physics test: a fully bonded 4-atom chain under
(0,0,0) special bonds must produce zero LJ energy. ac-004 is the
cross-backend parity check that closes Phase G.8 — CPU and CUDA must
agree across both the exclusion-on and exclusion-off configurations.

ac-005 protects users without bonded topology: the new code path must
not regress unbonded LJ runs.

ac-006 ensures lifecycle cleanliness — repeated build/destroy must not
leak device memory.

ac-007 — ac-008 are static / documentation gates verifiable without
running CUDA. ac-009 is the full-suite gate the chain ends on.
