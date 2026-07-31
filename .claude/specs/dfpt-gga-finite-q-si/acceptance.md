---
spec: dfpt-gga-finite-q-si
created: 2026-07-31
criteria:
  - id: ac-001
    summary: Finite-q PBE DFPT C00 matches offline supercell energy-FD pin
    type: runtime
    pass_when: |
      reg_dfpt-gga-finite-q-si runs DfptPhonon::compute with q_frac =
      (0.5, 0, 0) and XCType::PBE on the diamond-Si twin; the returned
      supercell force constant satisfies |C[(0x),(0x)] - FD pin| <= 0.30
      eV/A^2 where the pin is the offline central second difference of the
      restricted PBE PW SCF total energy on the 2x1x1 supercell
      (h = 0.03 A, tol 1e-9), hard-coded in REFERENCE.md.
    status: verified
  - id: ac-002
    summary: Physics window gates on the X-branch spectrum
    type: runtime
    pass_when: |
      Six finite modes at X; |omega| windows anchored by the CP2K 2026.1
      supercell VIBRATIONAL_ANALYSIS X-branch values (widened for the
      documented minimal-twin bias, final numbers in REFERENCE.md);
      n_imag bounded as in the LSDA parent; force-constant asymmetry and
      row-sum diagnostics at the family level.
    status: verified
  - id: ac-003
    summary: Regression is an out-of-library consumer, <= 60 s
    type: code
    pass_when: |
      regressions/response/dfpt-gga-finite-q-si/reg_dfpt-gga-finite-q-si.cpp
      uses public headers + link atomiverse only, registered via
      atv_add_regression, never in ctest; wall <= 60 s on the login node.
    status: verified
  - id: ac-004
    summary: REFERENCE.md + oracle note document a re-runnable pin
    type: docs
    pass_when: |
      REFERENCE.md cites method, twin settings, measured pin values with date
      and build type; tools/cp2k_oracle/inputs/dfpt-gga-finite-q-si.note
      records the offline recipe and the CP2K 2026.1 supercell run (binary,
      input, output line, wall) that anchors the physics window and the fair
      ratio (or a stated perf_exempt reason).
    status: verified
  - id: ac-005
    summary: Queue + roadmap close the row
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md moves dfpt-gga-finite-q-si to done;
      docs/roadmap/dft.md marks finite-q GGA DFPT shipped; tutorial method
      page indexed.
    status: verified
out_of_scope:
  - Joint spin x GGA x finite-q (census-held behind this row)
  - Dedicated incommensurate-q GGA pin (separate census decision)
  - Hybrid XC response; CUDA; live CP2K/QE at ctest time
---

# Acceptance — `dfpt-gga-finite-q-si`

**Done** means the finite-q (commensurate X) pure-GGA (PBE) DFPT phonon path
of `cpu::pw::DfptPhonon` is gated by an offline dual pin (same-Hamiltonian
supercell energy FD + CP2K-anchored physics window) in one <= 60 s
regression, with the queue row closed.
