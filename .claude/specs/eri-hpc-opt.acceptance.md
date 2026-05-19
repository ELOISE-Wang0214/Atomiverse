---
slug: eri-hpc-opt
criteria:
  - id: ac-001
    summary: Build succeeds with -march=native and xsimd 11.2.0 wired
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_TESTS=ON -DATV_BUILD_BENCHMARKS=ON
      -DATV_BUILD_LIBCINT=ON -DXTENSOR_USE_XSIMD=OFF` followed by `cmake --build build`
      completes with exit 0; xsimd 11.2.0 is fetched; no `XTENSOR_USE_XSIMD` compile
      definition is set globally on the `atomiverse` target.
    status: verified
    last_checked: 2026-05-19
  - id: ac-002
    summary: Ltot-dispatched eri_block + xsimd ss|ss path pass correctness gates
    type: code
    pass_when: |
      `ctest --test-dir build -R "BlockG|NativeIntegrals" --output-on-failure`
      passes with all BlockG cases to 1e-9 and NativeVsLibcint to 1e-7.
    status: verified
    last_checked: 2026-05-19
  - id: ac-003
    summary: Aspirin STO-3G Fock build beats PySCF direct
    type: performance
    pass_when: |
      `./build/benchmarks/bench_native --benchmark_filter=BM_GBuild_Block_Aspirin
      --benchmark_repetitions=5` reports median wall time ≤ 950 ms (PySCF baseline),
      with Phase 3 target ~500 ms.
    status: verified
    last_checked: 2026-05-19
  - id: ac-004
    summary: cc-pVDZ ERI matches libcint within 1e-7 (Phase 4)
    type: scientific
    pass_when: |
      A new cc-pVDZ test in `tests/scf/test_native_integrals.cpp` driving the
      BLAS angular contraction path agrees with libcint reference to max abs
      error ≤ 1e-7 over all (dd|dd), (dp|dd), (dd|pp), (dd|ss) shell quartets
      in the test molecule.
    status: pending
  - id: ac-005
    summary: cc-pVDZ Fock build benchmark recorded vs PySCF (Phase 4)
    type: performance
    pass_when: |
      `BM_GBuild_Block_Aspirin_ccPVDZ` exists in `bench_native` and emits a
      median wall time; the benchmark output lists PySCF cc-pVDZ direct-Fock
      baseline alongside, so caller can verify ≤ 1.0× ratio.
    status: pending
  - id: ac-006
    summary: HPC knowledge files document the canonical ERI patterns (Phase 5)
    type: docs
    pass_when: |
      `.claude/agents/av-hpc.md`, `.claude/skills/av-hpc/SKILL.md`, and
      `.claude/skills/av-perf/SKILL.md` each contain a section mentioning:
      xsimd 11.2.0 pin, `-march=native` language-guarded CMake pattern,
      Ltot-based eri_block dispatch, and the Boys SIMD pattern with safe-T select.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** are Phase 1-3 regression gates: any future change must keep build green, correctness within tolerance, and the perf win locked in.
- **ac-004 / ac-005** are Phase 4 gates (BLAS angular contraction for `l ≥ 2`, cc-pVDZ enablement).
- **ac-006** is the Phase 5 knowledge-propagation gate so the HPC agent and skills encode the patterns canonically.
