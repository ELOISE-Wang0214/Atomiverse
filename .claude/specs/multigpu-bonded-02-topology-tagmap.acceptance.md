---
slug: multigpu-bonded-02-topology-tagmap
criteria:
  - id: ac-001
    summary: Per-atom bonded SoA on cuda::Frame::Atoms, allocated/freed by capacity
    type: code
    pass_when: |
      cuda::Frame::Atoms declares num_bond/bond_atom(tags)/bond_type, num_angle/
      angle_atom/angle_type, num_dihedral/dihedral_atom/dihedral_type, nspecial/
      special; they are allocated to capacity_ on the multi-GPU bonded path and
      freed in destroy. The single-GPU flat path leaves them null and is unchanged.
    status: pending
  - id: ac-002
    summary: migrate carries per-atom bonded topology to the new owner
    type: scientific
    evaluator_hint: operator-action GH200, >=2 GPUs
    pass_when: |
      After an atom changes owner GPU via migrate, its num_bond/bond_atom(tags)/
      angle/dihedral/special arrays are present and identical (tag-wise) on the
      new owner — verified by a migrate-then-readback test.
    status: pending
  - id: ac-003
    summary: device tag->local map resolves owned+ghost, -1 for absent
    type: code
    pass_when: |
      The per-subdomain tag->local map (built over [0, n_local+n_ghost) after the
      halo) resolves a known owned tag and a known ghost tag to their correct
      local rows and returns -1 for a tag not present on the subdomain (unit test).
    status: pending
  - id: ac-004
    summary: re-materialized flat tuples match single-GPU set on 2 GPUs (no dup/loss)
    type: scientific
    evaluator_hint: operator-action GH200, >=2 GPUs
    pass_when: |
      For a small gen_pe system, the union over GPUs of the re-materialized
      bond/angle/dihedral tuples (partners mapped back to global tags) equals the
      single-GPU materialized tuple set; each tuple is owned by exactly one GPU
      (bond=min-tag, angle=center, dihedral=2nd atom); drop tally is 0 with an
      adequate halo cutoff.
    status: pending
  - id: ac-005
    summary: existing LJ multi-GPU tests stay green
    type: code
    pass_when: |
      test_cuda_multi_gpu (LJ, 13/13) remains green; the per-atom bonded SoA does
      not perturb the single-GPU or LJ multi-GPU paths.
    status: pending
---

# Acceptance criteria

- **ac-001** — per-atom bonded SoA on the CUDA frame, capacity-allocated, single-GPU path null/unchanged (code).
- **ac-002** — migrate carries the per-atom bonded topology to the new owner (scientific; operator GH200).
- **ac-003** — device tag→local map resolves owned+ghost, −1 for absent (code).
- **ac-004** — 2-GPU re-materialized tuple set == single-GPU set, exactly-one-owner, drop=0 (scientific; operator GH200).
- **ac-005** — existing LJ multi-GPU tests stay green (code).
