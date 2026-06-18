---
spec: native-xc-libxc-parity
title: Acceptance criteria — Native XC library with libxc/CP2K bit-parity
---

# Acceptance — Native XC (libxc & CP2K parity)

Oracle: libxc.so.15 + headers in the Spack view of the CP2K build at
`/nobackup/proj/disk/teoroo/personal/jicli594/work/cp2k/spack/spack/opt/spack/view`.
CP2K reference for PADE: local `cp2k` (`~/.local/bin/cp2k`, v2026.1).

---

ac-001: libxc oracle harness + sweep generator
  type: structural
  status: pending
  test_file: tests/cpu/xc/test_cpu_xc_oracle_sweep.cpp
  description: >
    A test-only oracle (ATV_BUILD_XC_ORACLE) calls xc_lda_exc_vxc / xc_gga_exc_vxc.
    A generator sweeps rho in [1e-10,1e3], zeta in [0,1] (polarized + unpolarized),
    reduced gradient s in [0,5] with edge cases (rho->0, fully polarized, sigma=0,
    large sigma). Harness runs and reports max relative error per quantity.

ac-002: LDA building blocks match libxc to 1e-12
  type: scientific
  status: pending
  test_file: tests/cpu/xc/test_cpu_xc_lda.cpp
  description: >
    LDA_X (Slater), LDA_C_PW (PW92), LDA_C_VWN (VWN5) native kernels match the libxc
    oracle with max relative error <= 1e-12 on exc, vrho (both spin channels) across
    the full sweep; any exceedance is a documented bounded last-ulp libm difference.

ac-003: PADE (CP2K-native) reproduces CP2K's XC energy
  type: scientific
  status: pending
  test_file: tests/cpu/xc/test_cpu_xc_pade.cpp
  description: >
    The GTH Pade LDA kernel, clean-roomed from the Goedecker-Teter-Hutter formula,
    reproduces CP2K's E_xc for tests/QS/regtest-gpw-1/Ar.inp and H2.inp (CP2K run
    locally) consistent with matching the published total energies to <= 1e-6 Ha.

ac-004: PBE matches libxc bit-for-bit
  type: scientific
  status: pending
  test_file: tests/cpu/xc/test_cpu_xc_pbe.cpp
  description: >
    GGA_X_PBE + GGA_C_PBE (spin-polarized, with zeta=0 unpolarized special case)
    match the libxc oracle to <= 1e-12 max relative error on exc, vrho, vsigma over
    the full sweep including high-gradient and low-density edge cases.

ac-005: Shared cpu/xc module; scf and gpw consumers unchanged
  type: structural
  status: pending
  test_file: tests/cpu/xc/test_cpu_xc_module_api.cpp
  description: >
    XC lives in cpu/xc/ with the existing evaluate_xc / evaluate_xc_batch signature.
    cpu::scf (hamiltonian_builder) and cpu::gpw (collocate/integrate path) both call
    it without source changes beyond the include path; build green with ATV_BUILD_RHF
    and ATV_BUILD_GPW.

ac-006: Drives the GPW tutorial parity gates
  type: scientific
  status: pending
  test_file: tests/cpu/gpw/test_cpu_gpw_m3_h2o_tutorial.cpp
  description: >
    With the native PBE kernel, GPW reproduces the getting-started h2o.inp total
    energy to <= 1e-5 Ha (the runtime-libxc dependency in gpw-gaussian-plane-wave D7
    is removed); with native PADE, GPW M1 matches Ar.inp/H2.inp to <= 1e-6 Ha.

ac-007: Vendored-source provenance + relicense recorded
  type: structural
  status: pending
  test_file: tests/cpu/xc/LICENSE_PROVENANCE.md
  description: >
    Vendored libxc kernels retain upstream MPL-2.0 headers and a provenance note
    (libxc version 15, file list); the PADE port records its CP2K xc_pade.F origin.
    Atomiverse's LICENSE/README are updated to GPL-2.0 to reflect the combined work.
