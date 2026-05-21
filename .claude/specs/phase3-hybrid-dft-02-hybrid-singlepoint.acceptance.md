---
slug: phase3-hybrid-dft-02-hybrid-singlepoint
criteria:
  - id: ac-001
    summary: B3LYP single-point energy matches PySCF to 1e-7 Ha on 5 molecules
    type: runtime
    pass_when: |
      test_hybrid_singlepoint runs B3LYP/6-31G* on H2O, CH4, NH3,
      benzene, glycine and every total energy is within 1e-7 Ha of the
      checked-in PySCF RKS(xc='b3lyp') reference in the *.json files;
      test passes under ctest.
    status: pending
  - id: ac-002
    summary: PBE0 single-point energy matches PySCF to 1e-7 Ha on 5 molecules
    type: runtime
    pass_when: |
      test_hybrid_singlepoint runs PBE0/6-31G* on H2O, CH4, NH3,
      benzene, glycine and every total energy is within 1e-7 Ha of the
      checked-in PySCF RKS(xc='pbe0') reference in the *.json files;
      test passes under ctest.
    status: pending
  - id: ac-003
    summary: HOMO/LUMO orbital energies match PySCF to 1e-5 eV
    type: runtime
    pass_when: |
      test_hybrid_singlepoint compares B3LYP and PBE0 HOMO and LUMO
      orbital energies on the 5-molecule set against the PySCF
      references; every value is within 1e-5 eV; test passes.
    status: pending
  - id: ac-004
    summary: alpha_hf=0 hybrid path reproduces pure-PBE bit-for-bit
    type: runtime
    pass_when: |
      a test case driving the hybrid Fock+energy path with alpha_hf
      forced to 0 yields a total energy bitwise identical to the
      existing pure-PBE SCF result on the same system; test passes.
    status: pending
  - id: ac-005
    summary: energy decomposition E_total == E_1e+E_J+alpha_hf*E_K+E_xc+E_nuc
    type: runtime
    pass_when: |
      test_hybrid_singlepoint asserts the hybrid SCFResult satisfies
      E_total == one_electron_energy + coulomb_energy + exchange_energy
      + xc_energy + nuclear_repulsion to within 1e-10 Ha; test passes.
    status: pending
  - id: ac-006
    summary: b3lyp_config / pbe0_config set correct functional and alpha_hf
    type: code
    pass_when: |
      b3lyp_config and pbe0_config exist in scf_solve.h; the SCFConfig
      they return has xc_type B3LYP / PBE0 respectively, and
      hybrid_exchange_fraction returns 0.20 for B3LYP and 0.25 for
      PBE0.
    status: pending
  - id: ac-007
    summary: full ctest suite stays green
    type: runtime
    pass_when: |
      ctest --test-dir build --output-on-failure passes with no
      regressions in existing HF / PBE / LSDA SCF tests.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** — Core correctness: hybrid total energies vs genuine
  checked-in PySCF reference data (with `pyscf_version` provenance) at
  6-31G* on the five-molecule set, tolerance 1e-7 Ha.
- **ac-003** — Orbital spectrum correctness: HOMO/LUMO energies vs PySCF,
  tolerance 1e-5 eV, confirming exact exchange is mixed correctly into the
  KS eigenproblem.
- **ac-004** — Reduction guarantee: the hybrid code path with `alpha_hf=0`
  must collapse exactly onto the pre-existing pure-PBE path (bit-for-bit),
  proving no `(1-alpha_hf)` rescaling artifact and no double counting.
- **ac-005** — Internal consistency: the energy decomposition closes
  exactly, proving `alpha_hf` is applied at one consistent site in both
  the Fock build and the energy accumulator.
- **ac-006** — User-facing entry points exist with the correct
  functional-specific `alpha_hf` constants.
- **ac-007** — No regression in the existing SCF test surface.
