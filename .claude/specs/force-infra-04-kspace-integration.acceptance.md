---
slug: force-infra-04-kspace-integration
criteria:
  - id: ac-001
    summary: Real/reciprocal split formalized with single g_ewald/alpha source
    type: code
    pass_when: |
      Real-space coulomb is a -03 pair-kernel; the KSpace family owns reciprocal only; the
      engine sets g_ewald/alpha once (kspace before pair). No duplicated authoritative alpha
      passed separately to CoulLong and Ewald.
    status: pending
  - id: ac-002
    summary: KSpace comm contract distinct from atom comm
    type: code
    pass_when: |
      The KSpace base defines a reciprocal comm contract: a global structure-factor reduction
      hook (Ewald) and grid forward/reverse hooks (PME), separate from the -02 atom
      forward/reverse comm. Ewald implements MPI_Allreduce (CPU) and NCCL allreduce (CUDA
      multi-GPU; single-GPU = no-op).
    status: pending
  - id: ac-003
    summary: Multi-GPU Ewald wired on MultiGpu with net-charge correction
    type: code
    pass_when: |
      MultiGpu registers KSpaceEwald via the -01 factory; charge is supplied by the -02
      contract; the force phase runs pair→bonded→kspace; the CUDA Ewald net-charge correction
      (previously omitted) is implemented.
    status: pending
  - id: ac-004
    summary: Single-GPU + CPU Ewald+Coul unchanged
    type: scientific
    evaluator_hint: operator-action — CPU + GH200 single-GPU
    pass_when: |
      test_ewald_coul total electrostatic energy, per-atom forces, and stress match baseline
      within documented tolerance on CPU and single-GPU.
    status: pending
  - id: ac-005
    summary: Multi-GPU Ewald == single-GPU == CPU (the kspace cross-domain proof)
    type: scientific
    evaluator_hint: operator-action — GH200 >=2 GPUs + CPU
    pass_when: |
      For a charged system, 2/4-GPU total electrostatic energy and per-atom forces match
      single-GPU and CPU within tolerance, proving charge-forward (-02) + cross-GPU S(k)
      reduction are correct.
    status: pending
  - id: ac-006
    summary: PME grid-comm interface defined; deferral logged; suites green
    type: code
    pass_when: |
      The PME grid forward/reverse interface exists and current single-GPU/replicated-grid
      behavior is preserved; multi-GPU distributed-grid PME is explicitly documented as a
      follow-up (no silent cap). Full CPU + CUDA ctest green.
    status: pending
---

# Acceptance criteria

- **ac-001** — real/recip split + single g_ewald/alpha source (code).
- **ac-002** — KSpace comm contract (S(k) reduction / grid) distinct from atom comm (code).
- **ac-003** — multi-GPU Ewald wired + net-charge correction (code).
- **ac-004** — single-GPU + CPU Ewald+Coul unchanged (scientific; operator).
- **ac-005** — multi-GPU Ewald == single == CPU (scientific; operator).
- **ac-006** — PME grid-comm interface + logged deferral; suites green (code).
