---
slug: phase4-vibrations-03-thermochemistry
criteria:
  - id: ac-001
    summary: ZPE matches PySCF within 0.001 eV on H2O/NH3/CH4
    type: runtime
    pass_when: |
      test_cpu_vib_thermochem_vs_pyscf passes its ZPE assertions: for H2O,
      NH3, and CH4 the Thermochem-computed ZPE differs from the checked-in
      PySCF thermo() reference (pyscf_thermo_references.json) by <= 0.001 eV.
    status: pending
  - id: ac-002
    summary: Thermal Gibbs correction matches PySCF within 0.5 kJ/mol
    type: runtime
    pass_when: |
      test_cpu_vib_thermochem_vs_pyscf passes its Gibbs assertions: the
      thermal correction G_tot - E0 at T=298.15 K, P=1 atm differs from the
      PySCF thermo() reference by <= 0.5 kJ/mol on H2O, NH3, CH4 (and
      benzene or formaldehyde if included).
    status: pending
  - id: ac-003
    summary: Per-component entropies match PySCF on H2O
    type: runtime
    pass_when: |
      test_cpu_vib_thermochem_vs_pyscf passes its component check: for H2O
      the Thermochem-computed S_trans, S_rot, and S_vib each match the
      corresponding PySCF thermo() per-component values within the tight
      per-component tolerance set in the test.
    status: pending
  - id: ac-004
    summary: FrequencyReport writes a parseable TSV with count and ZPE
    type: runtime
    pass_when: |
      test_cpu_vib_output_components passes: after a Thermochem run the
      FrequencyReport output file parses as TSV, has one data row per
      frequency, and its trailer reports the correct frequency count and
      the ZPE value.
    status: pending
  - id: ac-005
    summary: NormalModeWriter writes a valid loadable XYZ animation
    type: runtime
    pass_when: |
      test_cpu_vib_output_components passes: a NormalModeWriter output XYZ
      file loads successfully, each frame has the input atom count, and the
      atom symbols match the input molecule.
    status: pending
  - id: ac-006
    summary: Component roles are correctly tagged
    type: code
    pass_when: |
      Thermochem::role() returns ComponentRole::Observer and both
      FrequencyReport::role() and NormalModeWriter::role() return
      ComponentRole::Recorder, as asserted in test_cpu_vib_output_components.
    status: pending
  - id: ac-007
    summary: Full ctest suite stays green
    type: runtime
    pass_when: |
      cmake --build with -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
      succeeds and ctest --test-dir build reports 100% tests passed.
    status: pending
---

# Acceptance criteria

- **ac-001** — Zero-point energy is the simplest electronic-energy-independent
  thermochemistry sub-check (`E_0K - E0 == ZPE`). Matching PySCF to 0.001 eV
  (~0.1 kJ/mol) on the three small molecules confirms the harmonic vibrational
  sum and the frequency-to-theta conversion are correct.
- **ac-002** — The thermal Gibbs correction `G_tot - E0` exercises every
  partition function (translational, rotational, vibrational) and the total
  aggregation. The 0.5 kJ/mol bar is the published chain target; sigma is
  pinned per molecule in the fixtures so rotational `1/sigma` is exact.
- **ac-003** — Component-resolved entropy check isolates which partition
  function is wrong if a total fails: `S_trans` validates Sackur-Tetrode,
  `S_rot` validates the inertia-tensor / rotor-type path, `S_vib` validates
  the harmonic-oscillator sum.
- **ac-004 / ac-005** — Output-component round trips: the TSV must be
  machine-parseable with the frequency count and ZPE present; the XYZ
  animation must load with atom count and symbols intact.
- **ac-006** — Structural-role tags gate correct placement in the driver
  loop: Observer computes, Recorder writes I/O, never both.
- **ac-007** — No regression in the existing CPU test suite.
