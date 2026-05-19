---
slug: libxc-xc-fix-05-uks-validation
criteria:
  - id: ac-001
    summary: UKS-LSDA matches PySCF on all 5 validation systems within tolerance budget
    type: scientific
    pass_when: |
      In tests/cpu/scf/test_uks_solve.cpp, the new UKS-vs-PySCF reference
      assertions pass for every system under -DATV_USE_LIBXC=ON build:
        |E_atv - E_ref| < 1e-6 Ha for H atom (doublet)
        |E_atv - E_ref| < 1e-5 Ha for He atom (singlet)
        |E_atv - E_ref| < 1e-4 Ha for Li atom (doublet)
        |E_atv - E_ref| < 1e-5 Ha for H2 at R=0.74 Å
        |E_atv - E_ref| < 1e-4 Ha at every R in the 10-point H2 PES scan
      Reference energies are baked-in constexpr values from PySCF
      dft.UKS xc="lda,vwn", basis sto-3g, atom_grid=(50,110), prune=None,
      conv_tol=1e-12, with the version and invocation documented in the
      test file header.
    status: pending
  - id: ac-002
    summary: Native and libxc E_xc agree to 1e-12 Ha on converged UKS densities
    type: scientific
    pass_when: |
      Under -DATV_USE_LIBXC=ON, the #ifdef ATV_HAS_LIBXC-gated parity test
      in tests/cpu/scf/test_uks_solve.cpp passes for every system in
      {H, He, Li, H2@0.74, H2 scan points}: re-evaluating E_xc from the
      same converged (P_alpha, P_beta) on the same grid via
      cpu::scf::evaluate_xc and cpu::scf::evaluate_xc_libxc yields
      |E_xc_native - E_xc_libxc| < 1e-12 Ha.
    status: pending
  - id: ac-003
    summary: Full ctest suite passes under both -DATV_USE_LIBXC=ON and OFF
    type: runtime
    pass_when: |
      Both invocations exit 0 with zero failed tests:
        cmake -B build_on  -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_USE_LIBXC=ON
        cmake --build build_on -j$(nproc)
        ctest --test-dir build_on --output-on-failure
        cmake -B build_off -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_USE_LIBXC=OFF
        cmake --build build_off -j$(nproc)
        ctest --test-dir build_off --output-on-failure
      The OFF build skips the ac-002 parity block via #ifdef and still
      passes every other test in the suite (including the new UKS
      reference assertions from ac-001, which do not require libxc).
    status: pending
---

# Acceptance criteria

**ac-001 — UKS-LSDA vs PySCF reference (scientific).** Closes the original spec's success criterion #3 ("He atom energy matches PySCF to 1e-5 Ha") and extends it to the full five-system budget. PySCF reference values are hardcoded `constexpr` in the test file; the file header documents the exact PySCF version, basis, grid spec, and `dft.UKS` invocation used to generate them so the numbers stay reproducible.

**ac-002 — Native↔libxc E_xc parity (scientific).** Encodes Iron Rule #2 ("libxc as XC ground truth") at the energy-integral level: not just point-wise (spec 03's bar) but end-to-end on real converged UKS densities. The 1e-12 Ha bar proves the native XC fix is bit-identical to the reference wrapper after integration over the molecular grid.

**ac-003 — Build-matrix regression (runtime).** Closes the chain by asserting libxc remains genuinely optional. The OFF build must keep every non-parity test passing.
