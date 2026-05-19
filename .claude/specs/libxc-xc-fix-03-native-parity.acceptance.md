---
slug: libxc-xc-fix-03-native-parity
criteria:
  - id: ac-001
    summary: Native evaluate_xc matches the libxc oracle to 1e-12 on every battery point
    type: scientific
    pass_when: |
      With `-DATV_BUILD_XC_ORACLE=ON -DATV_BUILD_RHF=ON -DATV_BUILD_TESTS=ON`,
      the binary `build/tests/test_cpu_scf_xc_vs_libxc` runs to completion
      under ctest and reports all assertions passed. For each (functional,
      point) in the battery {LSDA, PBE} x {(0.5,0.5), (0.8,0.2),
      (1.0,0.0), (1e-6,1e-8), (100,50)} with PBE sigma sweep
      {0.0, 0.01, 0.1, 1.0}, every output channel of
      `ATV::cpu::scf::evaluate_xc(...)` (exc, vrho_alpha, vrho_beta, and
      for PBE vsigma_aa, vsigma_ab, vsigma_bb) satisfies
      `|native - oracle| <= 1e-12` in Hartree atomic units, where `oracle`
      is `ATV::test::xc_oracle::evaluate_lsda(rho_a, rho_b)` or
      `ATV::test::xc_oracle::evaluate_pbe(rho_a, rho_b, sigma_aa,
      sigma_ab, sigma_bb)` from sub-spec 02. Scientist-agent review at
      /mol:impl finalize confirms (a) units are Hartree/Bohr end-to-end,
      (b) cited references (VWN 1980 DOI 10.1139/p80-159, PBE 1996 DOI
      10.1103/PhysRevLett.77.3865, libxc DOI 10.1016/j.softx.2017.11.002)
      match the implemented formulas, (c) the vrho/vsigma convention
      adopted in evaluate_xc is libxc's d(rho*exc)/d(rho_sigma) form,
      and (d) every libxc-citation comment in xc_functional.cpp points
      at a real primitive in libxc/maple/ or libxc/src/.
    status: pending
  - id: ac-002
    summary: Parity target builds and discharges cleanly under both oracle gates
    type: code
    pass_when: |
      With the oracle enabled: `cmake -B build -DATV_BUILD_TESTS=ON
      -DATV_BUILD_RHF=ON -DATV_BUILD_XC_ORACLE=ON` configures without
      error; `cmake --build build --target test_cpu_scf_xc_vs_libxc
      -j$(nproc)` compiles cleanly; `ctest --test-dir build
      --output-on-failure -R '^(test_cpu_scf_xc_vs_libxc|test_cpu_scf_xc_functional)$'`
      exits 0. With the oracle disabled: `cmake -B build_off
      -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_XC_ORACLE=OFF`
      and `cmake --build build_off -j$(nproc)` succeed; `ctest --test-dir
      build_off --output-on-failure` exits 0; the parity test target is
      not registered; `nm libatomiverse` shows zero `xc_func` symbols;
      `grep -R "ATV_HAS_LIBXC\|ATV_USE_LIBXC" include/ src/` returns
      zero hits.
    status: pending
---

# Acceptance criteria

## ac-001 — Native vs. oracle parity (scientific)

The binding scientific gate of the chain. Every output channel of `evaluate_xc()` must reproduce the libxc-derived oracle's value to within 1e-12 Hartree. The scientist agent additionally verifies the in-source libxc citations are real — proving honest reimplementation, not opaque copy-paste.

## ac-002 — Build and CTest discharge (code)

Guards four regressions: (a) parity test fails to compile/link against the oracle when `-DATV_BUILD_XC_ORACLE=ON`; (b) oracle leaks into the production build when OFF; (c) `test_cpu_scf_xc_functional`'s updated expectations diverge from the corrected native code; (d) stray `ATV_HAS_LIBXC` / `ATV_USE_LIBXC` macros survive in production, silently re-creating a libxc dependency.
