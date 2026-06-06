---
slug: multigpu-bonded-03-bond
criteria:
  - id: ac-001
    summary: BondHarmonic evaluated over per-GPU tuples + reverse_comm wired
    type: code
    pass_when: |
      MultiGpu runs harmonic_bond_k over the re-materialized per-GPU bond tuples
      (deposit on owned + ghost partner) and calls reverse_comm() exactly once
      after the bonded force, before integration.
    status: verified
    last_checked: 2026-06-06
  - id: ac-002
    summary: 2/4-GPU bond energy + forces match single-GPU and CPU within tolerance
    type: scientific
    evaluator_hint: operator-action GH200, >=2 GPUs
    pass_when: |
      For a gen_pe system, 2- and 4-GPU E_bond and per-atom bond forces match the
      single-GPU run and the --backend cpu reference to <=1e-6 relative (or the
      documented mixed-precision tolerance).
    status: verified
    last_checked: 2026-06-06
  - id: ac-003
    summary: momentum conserved (reverse-once)
    type: scientific
    pass_when: |
      Sum of forces over all atoms stays ~0 each step on 2/4-GPU bond evaluation
      (no double application; ghost bond forces reversed exactly once).
    status: verified
    last_checked: 2026-06-06
  - id: ac-004
    summary: LJ multi-GPU tests stay green
    type: code
    pass_when: test_cuda_multi_gpu (LJ, 13/13) remains green.
    status: verified
    last_checked: 2026-06-06
---

# Acceptance criteria

- **ac-001** — BondHarmonic over per-GPU tuples + reverse_comm wired (code).
- **ac-002** — 2/4-GPU bond energy+forces == single == CPU (scientific; operator GH200).
- **ac-003** — ΣF≈0 reverse-once (scientific).
- **ac-004** — LJ multi-GPU tests green (code).
