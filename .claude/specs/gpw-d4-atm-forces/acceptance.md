---
spec: gpw-d4-atm-forces
created: 2026-07-31
criteria:
  - id: ac-001
    summary: dftd4() ATM leg is analytic (no FD, no grad_h)
    type: code
    pass_when: |
      src/cpu/scf/dftd4.cpp computes the s9 gradient in closed form (triple
      geometry derivatives + neutral-C6 CN chain); the grad_h parameter is
      removed from the public dftd4() signature and header docs updated.
    status: passed
  - id: ac-002
    summary: Periodic s9=1 forces match the dftd4 3.7.0 pins
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_dftd4.cpp gates max|F + g_pin| <= 5e-8
      Ha/Bohr (direct; above the oracle's own ~2.4e-8 periodic-ATM
      non-conservativeness, below the 4.05e-7 ATM increment), the ATM
      increment (s9=1 minus s9=0) vs the oracle pin delta to 5e-8 with
      non-vacuity > 1e-8, and F_x(O0) vs the FD of the oracle's printed
      displaced energies to 1e-9 Ha/Bohr.
    status: passed
  - id: ac-003
    summary: Isolated s9=1 forces match dftd4 3.7.0 --grad pin
    type: runtime
    pass_when: |
      Isolated water-dimer s9=1 forces match -gradient from
      pin_isolated_atm_grad.json to 1e-9 Ha/Bohr (oracle conservative when
      isolated).
    status: passed
  - id: ac-004
    summary: Analytic-vs-FD dual gate on the full s9=1 functional
    type: runtime
    pass_when: |
      test_cpu_scf_dispersion_gradient D4AtmForcesMatchFiniteDifference
      passes with the analytic path at <= 1e-7 Ha/Bohr (FD-step limited).
    status: passed
  - id: ac-005
    summary: reg_gpw-d4-atm-forces green
    type: runtime
    pass_when: |
      regressions/gpw/gpw-d4-atm-forces builds behind ATV_BUILD_REGRESSIONS,
      runs <= 60 s, checks e_disp vs the ATM pin (1e-8 Ha), standalone
      forces vs -grad pin (5e-8 direct + 1e-9 FD-parity Ha/Bohr), the ATM
      increment, and product wiring (1e-7 Ha/Bohr).
    status: passed
  - id: ac-006
    summary: Queue row closed
    type: docs
    pass_when: |
      docs/developer/cp2k-method-queue.md marks gpw-d4-atm-forces done with
      the measured numbers.
    status: passed
out_of_scope:
  - D3 ATM forces (sibling slug)
  - ATM stress route changes
  - CUDA
