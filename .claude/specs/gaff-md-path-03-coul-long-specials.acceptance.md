---
spec: gaff-md-path-03-coul-long-specials
created: 2026-07-21
criteria:
  - id: ac-001
    summary: CoulLong set_special_bonds API exists and scales real-space 1-4
    type: scientific
    pass_when: |
      test_coul_long_special is green: for a 1-4 charge pair,
      |E_real(s14=0.5) - 0.5*E_real(s14=1)| / max(|E_real(1)|,1e-12) ≤ 1e-10.
    status: pending
  - id: ac-002
    summary: 1-2 pairs with s12=0 contribute zero real-space Coulomb
    type: runtime
    pass_when: |
      Bonded 1-2 pair under set_special_bonds(0,0,0.5) yields real-space pair
      energy contribution 0 for that pair (isolated two-atom bonded fixture).
    status: pending
  - id: ac-003
    summary: NaCl PME regression unchanged without specials
    type: runtime
    pass_when: |
      regressions/md/npt_pme_nacl (or equivalent ctest target) still passes with
      the same tolerances documented in its REFERENCE.md.
    status: pending
  - id: ac-004
    summary: Header documents Amber usage with PME
    type: docs
    pass_when: |
      coul_long.h documents set_special_bonds and states that scales apply to
      real-space erfc only; Amber (0,0,0.5) example present.
    status: pending
out_of_scope:
  - Reciprocal-space special weighting
  - Full GAFF energy parity vs LAMMPS
  - Mandatory CUDA CoulLong specials
---

# Acceptance — gaff-md-path-03-coul-long-specials

Done means the real-space half of long-range Coulomb respects Amber special
weights, unblocking charged GAFF polymers under PME.
