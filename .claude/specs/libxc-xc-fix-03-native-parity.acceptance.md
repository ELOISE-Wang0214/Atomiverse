---
slug: libxc-xc-fix-03-native-parity
criteria:
  - id: ac-001
    summary: Native evaluate_xc matches libxc to 1e-12 on every enumerated point
    type: scientific
    pass_when: |
      With -DATV_USE_LIBXC=ON -DATV_BUILD_RHF=ON -DATV_BUILD_TESTS=ON, the
      binary build/tests/test_xc_vs_libxc runs to completion under ctest and
      reports all assertions passed. For each (functional, point) in the
      battery {LSDA, PBE} x {(0.5,0.5), (0.8,0.2), (1.0,0.0),
      (1e-6,1e-8), (100,50)} with PBE sigma sweep
      {0.0, 0.01, 0.1, 1.0}, every output channel (exc, vrho_alpha,
      vrho_beta, and for PBE vsigma_aa, vsigma_ab, vsigma_bb) satisfies
      |native - libxc| <= 1e-12 in Hartree atomic units, where libxc is
      invoked with functional IDs XC_LDA_X + XC_LDA_C_VWN (=7, VWN5) and
      XC_GGA_X_PBE + XC_GGA_C_PBE (=101, =130). Scientist-agent review at
      /mol:impl finalize confirms (a) units are Hartree/Bohr end-to-end,
      (b) cited references (VWN 1980 DOI 10.1139/p80-159, PBE 1996 DOI
      10.1103/PhysRevLett.77.3865, libxc DOI 10.1016/j.softx.2017.11.002)
      match the implemented formulas, and (c) the vrho/vsigma convention
      adopted in evaluate_xc is libxc's d(rho*exc)/d(rho_sigma) form.
    status: pending
  - id: ac-002
    summary: test_xc_vs_libxc target builds and discharges cleanly under ctest
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      -DATV_USE_LIBXC=ON` configures without error; `cmake --build build
      --target test_xc_vs_libxc -j$(nproc)` compiles with no warnings
      escalated to errors; `ctest --test-dir build --output-on-failure -R
      '^(test_xc_vs_libxc|test_xc_functional)$'` exits 0 with both tests
      reported as Passed. When libxc is *not* available
      (-DATV_USE_LIBXC=OFF), the same configure/build/ctest sequence
      still succeeds: test_xc_vs_libxc compiles to a trivial main and
      either is skipped at registration or returns success; no other
      target or test regresses.
    status: pending
---

# Acceptance criteria

## ac-001 — Native vs. libxc parity (scientific)

This is the binding scientific gate of the libxc-xc-fix chain. Every output channel of `evaluate_xc()` must reproduce libxc's value to within 1e-12 Hartree at every tabulated `(ρ_α, ρ_β, σ_αα, σ_αβ, σ_ββ)`. The tolerance is tight on purpose: 1e-12 is *numerical* equivalence in double precision, leaving no slack for "almost right" parameter tables or sign errors. The scientist agent verifies that implemented constants and formulas trace back to the cited primary literature and that the chosen libxc functional IDs are exactly the ones referenced (VWN5 = 7, not VWN3 = 5).

## ac-002 — Build and CTest discharge (code)

Guards against three regressions: (a) the new translation unit fails to compile under `-DATV_USE_LIBXC=ON`, (b) the libxc-guarded code leaks symbols when `-DATV_USE_LIBXC=OFF`, (c) `test_xc_functional`'s updated expectations diverge from the corrected native code.
