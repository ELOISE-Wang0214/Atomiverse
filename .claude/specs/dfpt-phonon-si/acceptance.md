---
spec: dfpt-phonon-si
created: 2026-07-25
criteria:
  - id: ac-001
    summary: Ewald dynamical matrix matches FD of Ewald forces
    type: runtime
    pass_when: |
      Unit test tests/cpu/planewave/test_dfpt.cpp: analytic
      ewald_dynamical_matrix vs central FD of EwaldSummation::forces_bohr on a
      2-atom periodic cell agrees elementwise to <= 5e-6 Ha/Bohr^2 (FD
      truncation-limited), and pair-accumulation row sums vanish to 1e-12.
    status: verified
    verified_by: test_cpu_pw_dfpt (EwaldDynamicalMatrix.*), 2026-07-25
  - id: ac-002
    summary: Sternheimer CG solves the free-electron response analytically
    type: runtime
    pass_when: |
      Unit test: with V_eff = 0 and no KB projectors, the projected CG
      solution equals Delta-psi(G) = -rhs(G)/(T(G) - eps) on the conduction
      space to <= 1e-8 per coefficient.
    status: verified
    verified_by: test_cpu_pw_dfpt (SternheimerCg.FreeElectronAnalytic), 2026-07-25
  - id: ac-003
    summary: DFPT force constant matches the discrete-energy second difference
    type: runtime
    pass_when: |
      Regression reg_dfpt-phonon-si: DFPT C[(0x),(0x)] agrees with the
      OFFLINE central second difference of the SCF total energy (identical
      LSDA Hamiltonian, ecut, k-mesh; h = 0.03 A, pinned in REFERENCE.md)
      within 0.30 eV/A^2. Revised from FD-of-analytic-forces: measurement
      showed pw_forces_analytic deviates ~2% from the exact discrete
      gradient (documented), so energy differences are the honest oracle.
      Measured agreement: 17.2277 vs 17.2393 eV/A^2 (0.07%).
    status: verified
    verified_by: run_reg_dfpt-phonon-si PASS, 2026-07-25
  - id: ac-004
    summary: Si Gamma spectrum — acoustic ~0, optical triplet in cited window
    type: runtime
    pass_when: |
      Regression: 3 lowest |omega| <= 30 cm^-1 with NO acoustic-sum-rule
      imposition (measured <= 4.4); 3 optical modes inside [500, 650] cm^-1
      around the literature DFPT-LDA Si TO(Gamma) 517 cm^-1 (Giannozzi et
      al., PRB 43, 7231 (1991)); triplet spread <= 10 cm^-1 (measured 2.7);
      row sums <= 0.02 eV/A^2. Values hard-coded; no live CP2K/QE.
    status: verified
    verified_by: run_reg_dfpt-phonon-si PASS, 2026-07-25
  - id: ac-005
    summary: Public API is members/ctor only; internals private
    type: code
    pass_when: |
      Public surface is cpu::pw::DfptPhonon members + compute(method, frame)
      returning DfptPhononResult; no public *Config bag, no make_*/build_*
      factory; Sternheimer/response machinery lives in
      planewave/internal/dfpt.{h,cpp} only.
    status: verified
    verified_by: agent-auto (surface review, 2026-07-25)
  - id: ac-006
    summary: Guard rails throw on unsupported regimes
    type: runtime
    pass_when: |
      Unit test: DfptPhonon::compute throws std::invalid_argument for
      xc != LSDA, empty upf_path, smearing != none, spin_polarized, and
      invalid mix_beta — all before any SCF or file access.
    status: verified
    verified_by: test_cpu_pw_dfpt (DfptPhonon.Throws*), 2026-07-25
  - id: ac-007
    summary: REFERENCE.md documents the offline pin and perf posture
    type: docs
    pass_when: |
      regressions/response/dfpt-phonon-si/REFERENCE.md cites the literature
      anchor + CP2K availability, input echo (cell, ecut, k-mesh, UPF, XC),
      measured ATV numbers and energy-FD pins, the FD-of-forces caveat, and
      perf_exempt=true with reason (no CP2K binary on the pin host).
    status: verified
    verified_by: agent-auto (REFERENCE.md review, 2026-07-25)
  - id: ac-008
    summary: pw_forces internal SCF honors method.xc (found defect fixed)
    type: code
    pass_when: |
      The promoted shared internal::run_scf_fixed_occ passes method.xc to
      XCBuilder::step and pw_forces_analytic_impl consumes it (previously
      its file-local run_scf silently evaluated the default PBE); the
      per-species v(G) form factor replaces the ill-defined S(G) division.
      Default-PBE callers unchanged: test_cpu_pw_forces / test_cpu_pw_stress
      still pass.
    status: verified
    verified_by: test_cpu_pw_forces + test_cpu_pw_stress green after fix, 2026-07-25
---

# Acceptance — dfpt-phonon-si

All criteria verified by the scoped gate: build `atomiverse`,
`test_cpu_pw_dfpt`, `test_cpu_pw_forces`, `test_cpu_pw_stress`,
`reg_dfpt-phonon-si`; ctest on those names; regression via
`--target run_reg_dfpt-phonon-si` (PASS, ≈ 25 s ≤ 60 s, hard-coded pins).
