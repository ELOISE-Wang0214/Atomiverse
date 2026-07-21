---
spec: gaff-md-path-01-io-full
created: 2026-07-21
criteria:
  - id: ac-001
    summary: full atom_style rows populate charge, type, coords
    type: runtime
    pass_when: |
      tests/cpu/io/test_lammps_data_full.cpp is registered and green; mini_full.data
      loads with nlocal matching file, each kAtom_charge equal to the file q to
      1e-12 abs, types 0-indexed, and xyz matching file Angstroms.
    status: pending
  - id: ac-002
    summary: molecular path remains byte-compatible for 6-column rows
    type: runtime
    pass_when: |
      mini_molecular.data loads without error; kAtom_type and coordinates match
      the molecular layout; no throw; existing consumers of load_lammps_data for
      molecular files still compile against the same API.
    status: pending
  - id: ac-003
    summary: bonds trigger build_special_list with non-trivial 1-2 entries
    type: runtime
    pass_when: |
      After loading a full fixture that has at least one bond, frame.span(kAtom_nspecial)
      has at least one atom with nspecial(*,0) > 0 (1-2 count), and kAtom_special
      partners are global tags of bonded neighbors.
    status: pending
  - id: ac-004
    summary: malformed Atoms rows throw with path context
    type: runtime
    pass_when: |
      A fixture with a truncated Atoms line causes load_lammps_data to throw
      std::runtime_error whose what() contains the path or "malformed".
    status: pending
  - id: ac-005
    summary: public header documents full and molecular support
    type: docs
    pass_when: |
      include/atomiverse/cpu/io/lammps_data.h @brief or @note names both
      atom_style full (with charge) and molecular; Doxygen-style @param path
      remains present.
    status: pending
out_of_scope:
  - Parsing pair_coeff / bond_coeff force-field include files
  - special_bonds energy scaling
  - molecule-ID frame column
  - CUDA I/O
---

# Acceptance — gaff-md-path-01-io-full

Done means a GAFF-style LAMMPS `atom_style full` data file can be loaded into
`cpu::Frame` with correct charges and topology specials scaffolding, while the
legacy molecular PE path still works.
