---
spec: gaff-md-path-04-ff-compose
created: 2026-07-21
criteria:
  - id: ac-001
    summary: kcal/mol to eV factor is single-sourced and tested
    type: scientific
    pass_when: |
      test_lammps_real pins convert(1 kcal/mol) = 1/23.060547830619026 eV to
      relative 1e-15; the factor is not redefined as a different literal in
      the new public header.
    status: pending
  - id: ac-002
    summary: Single-term fourier maps to DihedralCharmm parameters
    type: runtime
    pass_when: |
      A fourier term (K=0.25 kcal/mol, n=1, d=0 deg) converts to K_ev, n=1,
      delta=0 rad usable as DihedralCharmm constructor args; energy form
      identity is documented in the header.
    status: pending
  - id: ac-003
    summary: Multi-term fourier is rejected in this phase
    type: runtime
    pass_when: |
      Conversion API throws or returns error when m>1; test asserts the failure
      mode (no silent drop of higher terms).
    status: pending
  - id: ac-004
    summary: Angle degree to radian conversion is correct
    type: runtime
    pass_when: |
      theta0 = 109.5 deg converts to radians within 1e-15 relative of
      109.5 * pi/180.
    status: pending
out_of_scope:
  - Full peo25.ff text parser
  - Multi-term Fourier component
  - MD driver wiring
---

# Acceptance — gaff-md-path-04-ff-compose

Done means real-unit GAFF tables can be converted into Atomiverse potential
constructor arguments with tested unit factors and an explicit single-term
Fourier→Charmm map.
