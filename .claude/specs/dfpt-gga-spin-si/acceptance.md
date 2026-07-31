---
spec: dfpt-gga-spin-si
created: 2026-07-31
criteria:
  - id: ac-001
    summary: Joint spin x GGA DFPT C00 matches offline energy-FD pin
    type: runtime
    pass_when: |
      reg_dfpt-gga-spin-si runs DfptPhonon::compute with spin_polarized=true
      and XCType::PBE on the diamond-Si twin; |C[(0x),(0x)] - FD pin| <= 0.30
      eV/A^2 where the pin is the offline central second difference of the
      spin-polarized PBE SCF total energy (h = 0.03 A, tol 1e-9), hard-coded
      in REFERENCE.md.
    status: verified
  - id: ac-002
    summary: Physics window gates (acoustic, optical triplet, ASR, symmetry)
    type: runtime
    pass_when: |
      Acoustic |omega| <= 30 cm^-1 with no ASR imposition; optical triplet in
      [450, 650] cm^-1 with spread <= 10; DFPT row-sum <= 0.05 eV/A^2;
      pre-symmetrization asymmetry <= 1e-5 eV/A^2.
    status: verified
  - id: ac-003
    summary: Regression is an out-of-library consumer, <= 60 s
    type: code
    pass_when: |
      regressions/response/dfpt-gga-spin-si/reg_dfpt-gga-spin-si.cpp uses
      public headers + link atomiverse only, registered via
      atv_add_regression, never in ctest; wall <= 60 s on the login node.
    status: verified
  - id: ac-004
    summary: REFERENCE.md + oracle note document a re-runnable pin
    type: docs
    pass_when: |
      REFERENCE.md cites method, twin settings, measured pin values with date
      and build type; tools/cp2k_oracle/inputs/dfpt-gga-spin-si.note records
      the offline recipe and the CP2K 2026.1 run (binary, input, output line,
      wall) that anchors the physics window and the fair ratio.
    status: verified
  - id: ac-005
    summary: Queue + roadmap close the row
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md moves dfpt-gga-spin-si to done;
      docs/roadmap/dft.md marks joint spin x GGA DFPT shipped.
    status: verified
out_of_scope:
  - Finite-q GGA DFPT (held sibling dfpt-gga-finite-q-si)
  - Hybrid XC response; M != 0 cells; CUDA
  - Live CP2K/QE at ctest time
---

# Acceptance — `dfpt-gga-spin-si`

**Done** means the joint collinear-spin x pure-GGA (PBE) Γ DFPT phonon path of
`cpu::pw::DfptPhonon` is gated by an offline dual pin (same-Hamiltonian energy
FD + physics window) in one <= 60 s regression, with the queue row closed.
