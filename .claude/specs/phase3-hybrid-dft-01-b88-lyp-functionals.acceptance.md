---
slug: phase3-hybrid-dft-01-b88-lyp-functionals
criteria:
  - id: ac-001
    summary: VWN-RPA correlation matches libxc XC_LDA_C_VWN_RPA (id 8) to 1e-12
    type: runtime
    pass_when: |
      Built with -DATV_BUILD_RHF=ON -DATV_BUILD_XC_ORACLE=ON, the new
      VWN-RPA parity case in tests/cpu/scf/test_xc_vs_libxc.cpp passes:
      every channel (exc, vrho_alpha, vrho_beta) of vwn_rpa_correlation
      matches the libxc id-8 oracle within the combined abs/rel 1e-12
      bound over the full kPoints density battery including
      spin-polarized points.
    status: pending
  - id: ac-002
    summary: B88 exchange matches libxc XC_GGA_X_B88 (id 106) to 1e-12
    type: runtime
    pass_when: |
      The new B88 parity case in tests/cpu/scf/test_xc_vs_libxc.cpp
      passes: exc, vrho_alpha, vrho_beta, vsigma_aa, vsigma_ab, vsigma_bb
      of b88_exchange_spin match the libxc id-106 oracle within the
      combined abs/rel 1e-12 bound over the kPoints x kSigmaValues
      sweep including spin-polarized points.
    status: pending
  - id: ac-003
    summary: LYP correlation matches libxc XC_GGA_C_LYP (id 131) to 1e-12
    type: runtime
    pass_when: |
      The new LYP parity case in tests/cpu/scf/test_xc_vs_libxc.cpp
      passes: exc, vrho_alpha, vrho_beta, vsigma_aa, vsigma_ab, vsigma_bb
      of lyp_correlation match the libxc id-131 oracle within the
      combined abs/rel 1e-12 bound over the kPoints x kSigmaValues
      sweep including spin-polarized points.
    status: pending
  - id: ac-004
    summary: Full ctest suite stays green (no LSDA/PBE regression)
    type: runtime
    pass_when: |
      ctest --test-dir build --output-on-failure passes with zero
      failures, including the pre-existing XCvsLibxc.LsdaBattery,
      XCvsLibxc.PbeBattery and XCvsLibxc.DensityCutoff cases.
    status: pending
  - id: ac-005
    summary: Three new functionals present in native XC source
    type: code
    pass_when: |
      vwn_rpa_correlation is defined in src/cpu/scf/xc_functional.cpp;
      b88_exchange_spin and lyp_correlation are defined in the new file
      src/cpu/scf/xc_gga.cpp; all three are declared in the new private
      header include/atomiverse/cpu/scf/xc_internal.h; no symbol is
      named "vwn3".
    status: pending
  - id: ac-006
    summary: All three functionals carry Doxygen tags
    type: code
    pass_when: |
      grep of include/atomiverse/cpu/scf/xc_internal.h and the
      definitions of vwn_rpa_correlation / b88_exchange_spin /
      lyp_correlation shows each is preceded by a doxygen block with
      @brief and @param/@return where the signature has parameters /
      a return value.
    status: pending
---

# Acceptance criteria

- **ac-001 — VWN-RPA parity.** `vwn_rpa_correlation` is a re-parametrization of the existing VWN Pade machinery; correctness is verified solely by exact 1e-12 agreement with libxc functional id 8 on all LDA output channels.
- **ac-002 — B88 parity.** `b88_exchange_spin` per-spin GGA exchange must match libxc id 106 on energy density and all `vrho`/`vsigma` channels across the density and gradient sweep, including the empty-spin-channel zero limit.
- **ac-003 — LYP parity.** `lyp_correlation` spin-polarized GGA correlation must match libxc id 131 on all six output channels across the full sweep.
- **ac-004 — No regression.** The file split (moving shared VWN helpers behind `xc_internal.h`, adding `xc_gga.cpp`) must leave Slater / VWN5 / PW92 / PBE behavior bit-identical; the whole suite stays green.
- **ac-005 — Source presence.** Structural check that the three functionals exist in the native (non-oracle) engine, with the file split applied and the naming trap avoided.
- **ac-006 — Documentation.** Every new engine symbol is Doxygen-documented with units; verified by grep, hence `code`.
