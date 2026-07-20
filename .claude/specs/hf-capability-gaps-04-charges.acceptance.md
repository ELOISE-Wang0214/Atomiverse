---
slug: hf-capability-gaps-04-charges
criteria:
  - id: ac-001
    summary: Löwdin observer consumes MeanFieldOutputs/SlaterDeterminant, no debt types in API
    type: code
    pass_when: |
      include/atomiverse/cpu/observer/lowdin.h declares lowdin_charge/lowdin_spin
      taking scf::MeanFieldOutputs; grep finds no `IntegralPack` and no `SCFResult`
      token in lowdin.h or esp_charges.h public headers; all public symbols carry
      @brief/@param/@return with units [e].
    status: pending
  - id: ac-002
    summary: Löwdin unit tests pass without any SCF solve
    type: runtime
    pass_when: |
      tests/cpu/observer/test_lowdin.cpp is registered in ctest and green:
      single-s analytic case q = Z − P00 exact to 1e-12, symmetric two-atom
      synthetic density gives equal charges, sum(q) conserved; no run() of any
      mean-field method appears in the test file; densities come from a
      hand-built SlaterDeterminant.
    status: pending
  - id: ac-003
    summary: ESP constrained fit recovers synthetic charges on MK and CHELPG grids
    type: runtime
    pass_when: |
      tests/cpu/observer/test_esp_charges.cpp green: point charges
      (+0.4, −0.8, +0.4) at the H2O geometry recovered on both grid schemes,
      total-charge constraint satisfied exactly, symmetric H equal to 1e-8,
      buried-atom ill-conditioned case returns finite minimum-norm solution,
      unsupported-element vdW lookup throws.
    status: pending
  - id: ac-004
    summary: H2O charges regression matches PySCF-pinned values
    type: scientific
    pass_when: |
      `regressions` target runs regressions/scf/h2o_charges_pyscf in ≤60 s:
      Löwdin charges agree with PySCF-pinned values to ≤1e-6 e, CHELPG and MK
      each to ≤0.01 e, the two H charges equal; REFERENCE.md cites Löwdin 1950
      (10.1063/1.1747632), BMK 1990, Breneman-Wiberg 1990 and the pinned
      PySCF command/version.
    status: pending
  - id: ac-005
    summary: N2 ionization regression pins Koopmans and member-only ΔSCF
    type: scientific
    pass_when: |
      regressions/scf/n2_ionization green in ≤60 s: −ε_HOMO matches the
      PySCF-pinned value at the same basis to ≤1e-3 eV (≈17.3 eV HF band vs
      15.6 eV exp per Szabo-Ostlund §3.3, cited in REFERENCE.md); ΔSCF IP is
      computed via two converged runs configured only through the existing
      `charge` / `UHF::n_alpha` / `UHF::n_beta` members — grep of the diff
      shows no new *Config type and no factory function.
    status: pending
  - id: ac-006
    summary: OrbitalInfoObserver reports Koopmans IP column
    type: runtime
    pass_when: |
      tests/cpu/scf/test_orbital_info_observer.cpp green with the extension:
      ip_koopmans_eV() == −homo_eV() on a hand-built SlaterDeterminant and the
      TSV header gains an ip_koopmans column [eV].
    status: pending
  - id: ac-007
    summary: Tutorial ch02 Koopmans/ΔSCF and ch20 population/ESP section land
    type: docs
    pass_when: |
      docs/tutorial/ch02-rhf.md contains a hands-on Koopmans + ΔSCF subsection
      including an explicit warning that Koopmans EA is basis-set-artifactual
      (use ΔSCF); docs/tutorial/ch20-hf-properties-environment.md contains the
      "布居分析与 ESP 电荷" section with chapter-end PySCF comparison tables
      (scalars in tables, no bar charts) and a molplot/molvis ESP-grid figure.
    status: pending
  - id: ac-008
    summary: Full build and unit suite stay green in seconds
    type: runtime
    pass_when: |
      cmake --build build succeeds with -DATV_BUILD_HF=ON and
      `ctest --test-dir build --output-on-failure` passes with the new tests
      registered; no e2e/SCF solve added under tests/.
    status: pending
---

# Acceptance criteria

每条为二元判定;`scientific` 两条由 regressions 目标执行并对照 REFERENCE.md
中的引文与 PySCF 锚定值,`runtime` 各条由 ctest/构建判定,`docs` 由章节内容
存在性判定。SlaterDeterminant 对齐(scf-slater-determinant-03-unbag)由
ac-001/ac-002/ac-006 钉住。
