---
spec: phonon-dispersion-si
created: 2026-07-25
criteria:
  - id: ac-001
    summary: Public API is members/ctor only; no Config/factory
    type: code
    pass_when: |
      Public surface is cpu::vib::PhononDispersion members + compute(...)
      returning PhononDispersionResult; no public *Config bag; no make_*/build_*
      factory; Fourier + Hermitian diagonalization live in
      vibration/phonon_dispersion.{h,cpp}.
    status: verified
    verified_by: agent-auto (surface review, 2026-07-25)
  - id: ac-002
    summary: Unit tests cover D(q) wiring and diamond path
    type: runtime
    pass_when: |
      tests/cpu/vibration/test_phonon_dispersion.cpp: set_diamond_path length,
      compute rejects bad shapes/factors, Gamma acoustic near zero on a
      synthetic translationally-invariant IFC, frequencies finite and sorted.
    status: verified
    verified_by: test_cpu_vib_phonon_dispersion PASS (0.4 s), 2026-07-25
  - id: ac-003
    summary: SW Si q-path regression — Gamma / X / L physics windows
    type: runtime
    pass_when: |
      regressions/response/phonon-dispersion-si: diamond Si conventional
      cubic 8-atom + SW, 2×2×2 supercell IFCs via Phonon; PhononDispersion
      Gamma/X/L path; Gamma |acoustic| <= 30 cm^-1 (measured 0.00);
      Gamma optical max in [500, 900] cm^-1 (measured 594.82);
      X/L max omega in [300, 900] (measured 562.56 / 559.27); hard-coded
      pins only (no live CP2K).
    status: verified
    verified_by: run_reg_phonon-dispersion-si PASS, wall 0.07 s, 2026-07-25
  - id: ac-004
    summary: REFERENCE.md offline pin + perf posture
    type: docs
    pass_when: |
      REFERENCE.md cites method, system (conventional cubic SW twin), offline
      pins (Gamma/X/L windows + measured values), CP2K availability, and
      perf_exempt with reason.
    status: verified
    verified_by: agent-auto (REFERENCE.md review, 2026-07-25)
  - id: ac-005
    summary: Dual gate wall budget
    type: runtime
    pass_when: |
      reg_phonon-dispersion-si finishes in <= 60 s (measured 0.07 s); scoped
      build of atomiverse + test_cpu_vib_phonon_dispersion +
      reg_phonon-dispersion-si green; ctest -R test_cpu_vib_phonon_dispersion
      green; run_reg_phonon-dispersion-si PASS.
    status: verified
    verified_by: dual gate 2026-07-25
---

# Acceptance — phonon-dispersion-si

All criteria verified by the scoped gate: build `atomiverse`,
`test_cpu_vib_phonon_dispersion`, `reg_phonon-dispersion-si`; ctest on the
unit name; regression via direct binary / `run_reg_phonon-dispersion-si`
(PASS, 0.07 s ≤ 60 s, hard-coded pins). `perf_exempt=true` (no CP2K binary).
