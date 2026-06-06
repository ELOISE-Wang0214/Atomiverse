---
slug: multigpu-bonded-04-angle-dihedral
criteria:
  - id: ac-001
    summary: Angle + Dihedral evaluated cross-GPU with reverse_comm
    type: code
    pass_when: |
      MultiGpu runs harmonic_angle_k and multi_harmonic_dihedral_k over the
      re-materialized per-GPU tuples (deposit on all members incl. ghosts) with
      reverse_comm bringing ghost forces home.
    status: pending
  - id: ac-002
    summary: special-bonds cross-domain exclusion enabled on multi-GPU
    type: code
    pass_when: |
      multi-GPU setup_force_engine calls set_special_bonds and the exclusion list
      resolves excluded partner tags (including ghosts) to local rows, rebuilt
      each reneighbor.
    status: pending
  - id: ac-003
    summary: full PE energy + forces match single-GPU and CPU within tolerance
    type: scientific
    evaluator_hint: operator-action GH200, >=2 GPUs
    pass_when: |
      For a gen_pe system, 2/4-GPU full-PE energy decomposition
      (E_pair/E_bond/E_angle/E_dihedral) and per-atom forces match single-GPU and
      --backend cpu within tolerance; special-bonds-excluded pairs contribute
      identically; ΣF≈0 each step.
    status: pending
  - id: ac-004
    summary: LJ multi-GPU tests stay green
    type: code
    pass_when: test_cuda_multi_gpu (LJ, 13/13) remains green.
    status: pending
---

# Acceptance criteria

- **ac-001** — Angle + Dihedral cross-GPU + reverse_comm (code).
- **ac-002** — special-bonds cross-domain exclusion on multi-GPU (code).
- **ac-003** — full-PE energy+forces 2/4-GPU == single == CPU (scientific; operator GH200).
- **ac-004** — LJ multi-GPU tests green (code).
