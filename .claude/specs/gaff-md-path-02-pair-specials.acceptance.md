---
spec: gaff-md-path-02-pair-specials
created: 2026-07-21
criteria:
  - id: ac-001
    summary: PairLJCut scales 1-4 energy by s14 including 0.5
    type: scientific
    pass_when: |
      test_pair_special_scale reports |E(s14=0.5) - 0.5*E(s14=1)| / |E(s14=1)| ≤ 1e-12
      for a bonded topology that marks the interacting pair as 1-4 only; forces
      match the same scale to 1e-10 relative on the non-zero components.
    status: pending
  - id: ac-002
    summary: PairCoulCut scales 1-4 Coulomb the same way
    type: scientific
    pass_when: |
      Same relative 1e-12 energy / 1e-10 force test as ac-001 for PairCoulCut
      with non-zero charges on the 1-4 pair.
    status: pending
  - id: ac-003
    summary: s14=0 still fully excludes the pair
    type: runtime
    pass_when: |
      With set_special_bonds(0,0,0) or s14=0 on a 1-4-only pair, pair energy is
      exactly 0 and forces are 0 (abs ≤ 1e-15 eV / eV/Å).
    status: pending
  - id: ac-004
    summary: Existing (0,0,0) TIP3P-style regressions stay green
    type: runtime
    pass_when: |
      regressions/forcefield/tip3p_two_water and regressions/md/rigid_water (or
      their unit stand-ins that set_special_bonds(0,0,0)) still pass under ctest.
    status: pending
  - id: ac-005
    summary: Headers document Amber weights
    type: docs
    pass_when: |
      pair_lj_cut.h and pair_coul_cut.h mention that scales multiply energy and
      force, and give set_special_bonds(0.0, 0.0, 0.5) as the Amber/GAFF example.
    status: pending
out_of_scope:
  - CoulLong specials
  - CUDA mixed special scales
  - Separate lj vs coul weight tables inside one component
---

# Acceptance — gaff-md-path-02-pair-specials

Done means CPU cut pair styles implement true LAMMPS special_bonds scaling so
Amber 1–4 half-weight is expressible and tested.
