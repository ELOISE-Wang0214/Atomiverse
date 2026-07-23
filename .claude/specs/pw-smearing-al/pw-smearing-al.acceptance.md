---
spec: pw-smearing-al
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public PW smearing API is members-only (no Config bag)
    type: code
    pass_when: |
      cpu::pw::PW exposes smearing, degauss [Ry], mp_order, extra_bands as
      public members; PWResult exposes e_fermi, n_elec_occ, e_entropy.
      No public PWSmearingConfig / make_smearing / smearing_config.
    status: met
  - id: ac-002
    summary: Unit tests for occupation formulas and Fermi search
    type: runtime
    pass_when: |
      tests/cpu/planewave/test_pw_smearing (test_cpu_pw_smearing) exits 0.
      FD/Gaussian/MP-1 f(x) match closed form to 1e-12; Fermi search recovers
      Ne within 1e-8 on a synthetic spectrum.
    status: met
  - id: ac-003
    summary: Al FCC FD-smeared PW regression vs offline pin
    type: runtime
    pass_when: |
      regressions/pw/pw-smearing-al exits 0: SCF converged; |n_elec_occ - 3|
      < 1e-5; e_tot [Ry] within REFERENCE.md tol of offline pin; e_fermi finite
      and between band edges; occupations not all exactly 0 or 2.
    status: met
  - id: ac-004
    summary: Offline CP2K twin documented
    type: docs
    pass_when: |
      REFERENCE.md cites CP2K version, SMEAR method, geometry, and records
      pin numbers + perf note (within_10x or perf_exempt).
    status: met
  - id: ac-005
    summary: Default fixed occupations still work
    type: runtime
    pass_when: |
      Existing PW unit tests (test_cpu_pw_method / si) still pass with
      default Smearing::none.
    status: met
out_of_scope:
  - Spin-polarized metals
  - CUDA PW smearing
  - GPW SMEAR productization
  - Cold / Marzari–Vanderbilt smearing
---

# Acceptance — pw-smearing-al

Done means metallic FCC Al runs plane-wave SCF with Fermi–Dirac (and unit-tested
MP/Gaussian) fractional occupations, public API is concrete members only, and
the e2e regression passes against the offline pin in REFERENCE.md.
