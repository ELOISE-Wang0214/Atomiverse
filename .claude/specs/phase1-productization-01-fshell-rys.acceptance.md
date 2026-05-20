---
slug: phase1-productization-01-fshell-rys
criteria:
  - id: ac-001
    summary: kRysMaxL==3 and kRysMaxRoots==7 in rys_eri.h
    type: code
    pass_when: |
      grep in include/atomiverse/cpu/scf/rys_eri.h shows
      `constexpr int kRysMaxL = 3;` and
      `constexpr int kRysMaxRoots = 7;`
    status: pending
  - id: ac-002
    summary: kMaxL==3 in cart.h
    type: code
    pass_when: |
      grep in include/atomiverse/cpu/scf/cart.h shows
      `inline constexpr int kMaxL = 3;`
    status: pending
  - id: ac-003
    summary: Stale STO-3G nroots<=5 comments removed
    type: code
    pass_when: |
      No comment in rys_eri.h or src/cpu/scf/rys_roots.cpp still
      claims the engine supports "only" l in {0,1} or nroots <= 5;
      grep for "kRysMaxRoots = 5" and "nroots <= 5" returns no
      live (non-historical) assertion of the old limit.
    status: pending
  - id: ac-004
    summary: rys_roots() dispatches nroots 6 and 7 without aborting
    type: runtime
    pass_when: |
      A unit test in tests/cpu/scf/test_rys_roots.cpp calls
      rys_roots() with nroots=6 and nroots=7 across small/
      intermediate/large x and the returned roots and weights
      match an independent reference to <= 1e-10; the test
      passes under ctest.
    status: pending
  - id: ac-005
    summary: f-quartet ERIs are finite and non-zero
    type: runtime
    pass_when: |
      tests/cpu/scf/test_rys_eri_fshell.cpp computes
      rys_eri_block for (f,s|s,s), (f,p|s,s), (f,d|s,s) and
      (f,f|f,f); every block entry is finite and at least one
      entry per quartet is non-zero; the test passes under ctest.
    status: pending
  - id: ac-006
    summary: f-quartet ERI derivatives pass FD self-consistency
    type: runtime
    pass_when: |
      For each of (f,s|s,s), (f,p|s,s), (f,d|s,s), (f,f|f,f),
      rys_eri_block_deriv agrees with central finite differences
      of rys_eri_block (h=1e-4) to <= 1e-6 relative (or 1e-7
      absolute) on every component; the test passes under ctest.
    status: pending
  - id: ac-007
    summary: f-quartet ERI matches libcint reference
    type: runtime
    evaluator_hint: gated on ATV_HAS_LIBCINT and l=3 harness support
    pass_when: |
      With -DATV_BUILD_LIBCINT=ON, at least one f-quartet ERI
      from rys_eri_block matches build_eri_libcint_from_pack
      (l=3 g_l branch added) to <= 1e-9; test passes under ctest.
      If the harness cannot be made l=3-correct in this link the
      criterion is satisfied instead by a code check that the
      l=3 g_l branch is present in tests/libcint_reference.cuh.
    status: pending
  - id: ac-008
    summary: def2-TZVP build_integrals succeeds end-to-end
    type: runtime
    pass_when: |
      A test builds integrals for a small molecule (CH4 or H2O)
      with the def2-TZVP basis; build_integrals returns an
      IntegralPack whose ERI tensor is fully finite and non-zero,
      exercising f-shells; the test passes under ctest.
    status: pending
  - id: ac-009
    summary: New public symbols carry doxygen tags
    type: code
    pass_when: |
      grep of src/cpu/scf/rys_roots_iterative.h shows every public
      function/constant declaration is preceded by a doxygen block
      containing @brief, and @param/@return where the signature has
      parameters / a return value.
    status: pending
  - id: ac-010
    summary: Full ctest suite stays green, no l<=2 regression
    type: runtime
    pass_when: |
      `cmake --build build` succeeds and `ctest --test-dir build`
      reports 100% pass, including all pre-existing s/p/d ERI,
      gradient, and SCF tests.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** — the constant bumps and comment hygiene. Pure grep checks; `kMaxNCart` follows `kMaxL` automatically and needs no separate criterion.
- **ac-004** — the load-bearing new-algorithm gate: the iterative nroots-6/7 root finder must produce correct roots/weights, not merely not crash. RED test exists before the port lands.
- **ac-005 / ac-006** — f-shell ERI correctness. ac-005 is the GREEN counterpart to "write failing f-quartet tests"; ac-006 is the FD self-consistency domain check that exercises the +1 derivative shift and 7-root budget — a deterministic ctest comparison, hence `runtime`.
- **ac-007** — absolute correctness against libcint. Typed `runtime` with a documented fallback to a `code` check because the harness's l=3 support is itself part of this link's scope.
- **ac-008** — the end-user-visible outcome: a real f-shell basis (def2-TZVP) runs through `build_integrals`.
- **ac-009** — doxygen contract for the only genuinely new public header; verified by grep, hence `code`.
- **ac-010** — the no-regression guarantee; the `kMaxL` bump must leave every l≤2 result bit-identical.
