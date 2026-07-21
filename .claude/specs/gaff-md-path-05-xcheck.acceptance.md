---
spec: gaff-md-path-05-xcheck
created: 2026-07-21
criteria:
  - id: ac-001
    summary: GAFF fragment single-point total energy matches pin
    type: scientific
    pass_when: |
      regressions target gaff_peo_fragment exits 0; |E - E_pin| / max(|E_pin|,1e-8)
      ≤ 1e-4 (or the tighter absolute bound stated in REFERENCE.md in eV).
    status: pending
  - id: ac-002
    summary: Bonded and nonbonded partial energies reported and within term tolerance
    type: scientific
    pass_when: |
      REFERENCE lists bonded vs nonbonded (or finer) pins; regression checks each
      published term at relative ≤ 5e-4 unless REFERENCE documents a looser PME
      term with justification.
    status: pending
  - id: ac-003
    summary: Fixture is in-repo and path-stable
    type: code
    pass_when: |
      fragment.data (or equivalent) lives under regressions/forcefield/gaff_peo_fragment/
      with no absolute /home or /nobackup paths required to run the binary.
    status: pending
  - id: ac-004
    summary: REFERENCE documents how pins were produced
    type: docs
    pass_when: |
      REFERENCE.md states cutoff, alpha or kspace settings, special_bonds amber,
      unit conversion, and either LAMMPS version/command or analytic derivation.
    status: pending
out_of_scope:
  - Full Gudla-Zhang multi-ns protocol
  - Live LAMMPS in CI
  - N=100 production parity
---

# Acceptance — gaff-md-path-05-xcheck

Done means a small GAFF PEO-like system has a checked-in single-point energy
oracle that the Atomiverse classical stack matches within documented tolerances.
