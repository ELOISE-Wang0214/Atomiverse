---
slug: libxc-xc-fix-02-libxc-wrapper
criteria:
  - id: ac-001
    summary: xc_libxc_oracle static library builds and references xc_func_init
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      -DATV_BUILD_XC_ORACLE=ON` configures without error and
      `cmake --build build --target xc_libxc_oracle -j` exits 0.
      Running `nm $(find build -name 'libxc_libxc_oracle.a' -o -name
      'xc_libxc_oracle.cpp.o' | head -1)` shows at least one reference
      to a symbol matching `xc_func_init`, confirming the oracle
      actually calls libxc.
    status: pending

  - id: ac-002
    summary: libatomiverse contains zero libxc symbols regardless of oracle option
    type: code
    pass_when: |
      For both configurations
      `ATV_BUILD_XC_ORACLE=OFF` and `ATV_BUILD_XC_ORACLE=ON`, the command
      `nm $(find build -name 'libatomiverse.*' | head -1) 2>/dev/null
      | grep -c 'xc_func\|xc_oracle::evaluate_lsda\|xc_oracle::evaluate_pbe'`
      returns `0`. The oracle never leaks into the production library.
    status: pending

  - id: ac-003
    summary: Oracle reproduces Slater exchange textbook value to 1e-10 Ha
    type: scientific
    pass_when: |
      In the `test_xc_libxc_oracle` binary, the LDA-exchange-only
      sub-test evaluates libxc `XC_LDA_X` alone at
      `rho_alpha = rho_beta = 0.5` (zeta = 0) in `XC_POLARIZED` mode
      and asserts `std::abs(exc - (-0.7385587663820224)) < 1e-10`
      Hartree, matching the closed form `-(3/4)(3/pi)^(1/3)`. The
      full `ctest -R test_xc_libxc_oracle` exit code is 0.
    status: pending

  - id: ac-004
    summary: Default-OFF build produces no oracle artifacts
    type: code
    pass_when: |
      With `cmake -B build_off -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON`
      (i.e. `ATV_BUILD_XC_ORACLE` unset, default OFF) followed by
      `cmake --build build_off -j`,
      `find . -path ./build_off -prune -o -name 'xc_libxc_oracle*' -print`
      finds the source files under `tests/cpu/scf/oracle/`, while
      `find build_off -name 'xc_libxc_oracle*'` returns nothing — no
      object file, no static library, no test binary is emitted.
    status: pending
---

# Acceptance criteria

- **ac-001** asserts the oracle actually links against libxc when its build option is enabled; the `nm` check prevents a silent-empty-TU regression.
- **ac-002** is the load-bearing architectural invariant: `libatomiverse` never contains libxc symbols. This is what makes the oracle deletable later.
- **ac-003** is the numerical anchor: `−(3/4)(3/π)^{1/3}` is the canonical cross-library reference for LDA exchange.
- **ac-004** enforces the default-OFF lifecycle property. Together with ac-002, retirement is `rm -rf tests/cpu/scf/oracle/` plus one CMake-option removal — no production surgery.
