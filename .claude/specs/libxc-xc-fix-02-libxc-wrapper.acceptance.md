---
slug: libxc-xc-fix-02-libxc-wrapper
criteria:
  - id: ac-001
    summary: Build succeeds with ATV_USE_LIBXC=ON and links libxc symbols
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      -DATV_USE_LIBXC=ON && cmake --build build -j` completes with exit
      code 0, the resulting `libatomiverse` exports
      `ATV::cpu::scf::evaluate_xc_libxc`, and `xc_libxc.cpp` references
      `xc_func_init` (verified via `nm`/`objdump` on the object file).
    status: pending

  - id: ac-002
    summary: Build with ATV_USE_LIBXC=OFF still succeeds and emits no libxc symbols
    type: code
    pass_when: |
      `cmake -B build_off -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      -DATV_USE_LIBXC=OFF && cmake --build build_off -j` succeeds, and
      `nm build_off/.../xc_libxc.cpp.o` shows no `xc_func_init` or
      `evaluate_xc_libxc` symbol (the TU is empty without ATV_HAS_LIBXC).
    status: pending

  - id: ac-003
    summary: test_xc_libxc passes — finite XCOutput for LSDA and PBE at ρ=1
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_xc_libxc --output-on-failure`
      passes (exit 0), including the happy-path LSDA, happy-path PBE,
      spin-polarized smoke, and density-cutoff sub-tests defined in
      *Testing strategy*.
    status: pending

  - id: ac-004
    summary: libxc XC_LDA_X reproduces textbook Slater exchange at ρ=1 to 1e-10 Ha
    type: scientific
    pass_when: |
      In `test_xc_libxc`, the LDA-exchange-only sub-test evaluates libxc
      `XC_LDA_X` at ρ_α = ρ_β = 0.5 (ζ = 0) in `XC_POLARIZED` mode and
      asserts `|exc - (-0.7385587663820224)| < 1e-10` Hartree, matching
      the closed-form `-(3/4)(3/π)^{1/3}` independent of VWN.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** establish the build contract: the wrapper compiles and links with libxc when requested, and disappears cleanly when not.
- **ac-003** is the runtime smoke contract: every code path in `evaluate_xc_libxc` (LSDA branch, PBE branch including `vsigma_*`, spin-polarized init, cutoff) is exercised at least once.
- **ac-004** is the only numerical anchor this sub-spec verifies without depending on the native side. The Slater-exchange value `−(3/4)(3/π)^{1/3}` is the canonical cross-library reference for LDA exchange and confirms libxc is driven with the correct convention (`XC_POLARIZED`, density indices, `exc`-per-electron units).
