---
slug: multi-gpu-atom-cell-sort-01-permute
criteria:
  - id: ac-001
    summary: Cell-key permutation reorders all owned per-atom arrays in place
    type: code
    pass_when: |
      A unit test in tests/cuda/md/ builds a frame with known positions +
      full bonded SoA, applies the cell-sort permute on one device, and
      asserts every LJ array (x,y,z,vx,vy,vz,fx,fy,fz,mass,charge,type,tag,
      image) AND every bonded member (num_bond, bond_atom, bond_type,
      num_angle, angle_atom, angle_type, num_dihedral, dihedral_atom,
      dihedral_type, nspecial, special) is reordered by the same permutation.
    status: verified
    last_checked: 2026-06-08
  - id: ac-002
    summary: Permute round-trips to identity under inverse permutation
    type: code
    pass_when: |
      Applying the permutation then its inverse restores the original frame
      byte-for-byte (asserted in the same test).
    status: verified
    last_checked: 2026-06-08
  - id: ac-003
    summary: Sorted atoms are non-decreasing in neighbor-grid cell index
    type: code
    pass_when: |
      After the permute, the test reads back per-atom cell indices and
      asserts they are sorted non-decreasing.
    status: verified
    last_checked: 2026-06-08
---
