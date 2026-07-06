---
slug: hf-koopmans
criteria:
  - id: ac-001
    summary: Full 5×2 −ε_HOMO sweep matches pinned PySCF values to 1e-8 Ha (E2E tier)
    type: scientific
    pass_when: |
      With -DATV_RUN_E2E_TESTS=ON, test_cpu_scf_koopmans_parity passes: for all
      10 (system, basis) pairs {He,H2,H2O,NH3,CH4} × {STO-3G, 6-31G} at the
      pinned geometries, |(−ε_HOMO)_engine − ref| < 1e-8 Ha against the pinned
      table (STO-3G: 0.87603551, 0.57797481, 0.39123677, 0.35252132,
      0.51931445; 6-31G: 0.91412663, 0.59539276, 0.50136813, 0.41503216,
      0.54404220).
    status: pending
  - id: ac-002
    summary: RHF total energies anchor the same solves to 1e-8 Ha (E2E tier)
    type: scientific
    pass_when: |
      test_cpu_scf_koopmans_parity also asserts, for all 10 pairs,
      |E_RHF_engine − ref| < 1e-8 Ha against the pinned totals (STO-3G:
      −2.8077839575, −1.1166843871, −74.9630231385, −55.4540871348,
      −39.7268336369; 6-31G: −2.8551604262, −1.1267339671, −75.9839744727,
      −56.1610211797, −40.1805127401), and asserts the conversion constant
      27.211386245988 eV/Ha in one place.
    status: pending
  - id: ac-003
    summary: ≥3 anchor systems pass inside the default ctest gate
    type: code
    pass_when: |
      A plain default-config ctest run (no ATV_RUN_E2E_TESTS) executes
      test_cpu_scf_koopmans and it passes with at least the four anchors
      He/H2/H2O @ STO-3G and He @ 6-31G at 1e-8 Ha, including the
      lumo_eV()==homo_eV() no-virtuals assertion for He/STO-3G.
    status: pending
  - id: ac-004
    summary: Example 11_koopmans builds and prints the Koopmans table (Ha and eV)
    type: runtime
    pass_when: |
      With -DATV_BUILD_EXAMPLES=ON -DATV_BUILD_RHF=ON, target 11_koopmans
      compiles and, when run, prints one row per (system, basis) pair
      containing −ε_HOMO in Ha, the Koopmans IP in eV, and an experimental IP
      column labelled adiabatic or vertical; exit code 0.
    status: pending
  - id: ac-005
    summary: Example demonstrates OrbitalInfoObserver, gated behind ATV_BUILD_RHF
    type: code
    pass_when: |
      examples/11_koopmans.cpp constructs cpu::scf::OrbitalInfoObserver from
      SCFEvaluator::result_handle() and reads homo_eV(); its add_executable
      registration in examples/CMakeLists.txt sits inside the
      if(ATV_BUILD_RHF) block; examples/README.md has a row 11.
    status: pending
  - id: ac-006
    summary: Tutorial numbers are engine-derived and regression-pinned
    type: docs
    pass_when: |
      Every numeric −ε_HOMO / IP value printed in the new Koopmans section of
      docs/tutorial/ch02-rhf.md appears verbatim (same digits, Ha, or as the
      documented eV conversion thereof) in the pinned constexpr reference
      table asserted by tests/cpu/scf/test_koopmans_parity.cpp, and the two
      files carry mutual cross-reference comments.
    status: pending
  - id: ac-007
    summary: Tutorial section covers the required Jensen-style content
    type: docs
    pass_when: |
      The ch02 Koopmans section contains: the frozen-orbital derivation
      ending in ε_i; the relaxation-vs-correlation error-cancellation
      argument; the adiabatic-vs-vertical IP distinction applied to the NIST
      values (incl. NH3 planarization and CH4 Jahn–Teller caveats); the H2O
      10.6→13.6 eV STO-3G→6-31G sensitivity; a limitations paragraph; and the
      Koopmans 1934 reference with DOI 10.1016/S0031-8914(34)90011-2.
    status: pending
  - id: ac-008
    summary: Default ctest gate respected — parity test disabled by default
    type: runtime
    pass_when: |
      In a default configure (ATV_RUN_E2E_TESTS=OFF), test_cpu_scf_koopmans_parity
      is compiled but listed as Disabled by ctest (ATV_E2E_TESTS tail
      mechanism), and the added default-gate wall time from
      test_cpu_scf_koopmans is under ~2 s.
    status: pending
---

# Acceptance criteria

ac-001/ac-002 是科学契约：教程表格里的每个数都由引擎重算并对 PySCF 2.13.1（cart=True, conv_tol=1e-12）钉死到 1e-8 Ha。ac-003/ac-008 保护 ~16 s 默认闸门分层。ac-004/ac-005 保证示例既能跑又确实示范了现有 `OrbitalInfoObserver` 而非绕过它。ac-006 是"教程数字即回归资产"的绑定条款；ac-007 钉住 Jensen 式教学要素与绝热/垂直辨析不被删减。
