---
slug: force-infra-03-pair-kernel-seam
criteria:
  - id: ac-001
    summary: Backend-agnostic, potential-neutral PairKernel concept (multi-contribution, templated T)
    type: code
    pass_when: |
      A PairKernel concept + contribution type exists in a CUDA-free header, templated on scalar
      T (double, xsimd::batch<double>, CUDA real_t). Its name/shape is generic (NOT LJExpr / not
      LJ-specific) and supports one-or-more force/energy contributions per pair.
    status: pending
  - id: ac-002
    summary: Existing LJ AND Coulomb-real-space expressed as kernel instances
    type: code
    pass_when: |
      Both the LJ and the real-space coulomb (erfc) force math are kernel instances of the same
      concept — proving the seam is not bound to LJ. No new physics is introduced.
    status: pending
  - id: ac-003
    summary: Shared per-backend launchers; PairLJCut + CoulLong ported on both
    type: code
    pass_when: |
      CUDA pair_compute<Kernel> (Verlet traversal+reduction) and CPU pair_compute<Kernel>
      (xsimd ket-loop) exist; PairLJCut and CoulLong on both backends are reduced to a kernel +
      launcher call (no per-style neighbor loop / reduction). The CPU(ComponentBase)/CUDA(Pair)
      CoulLong base asymmetry is unified.
    status: pending
  - id: ac-004
    summary: Superposition — combined LJ+coul kernel equals separate runs
    type: scientific
    evaluator_hint: operator-action — CPU + GH200
    pass_when: |
      A combined LJ+coul kernel evaluated in one neighbor walk yields per-atom forces and
      energy equal (within tol) to running LJ and CoulLong separately and summing.
    status: pending
  - id: ac-005
    summary: Migrated LJ + Coul bit-identical to baseline
    type: scientific
    evaluator_hint: operator-action — CPU n1/n8 + GH200 1/2/4-GPU
    pass_when: |
      Per-atom forces, energy, virial from kernel-based PairLJCut and CoulLong match the
      pre-change hand-written kernels bit-for-bit (FP64) or within documented mixed precision.
    status: pending
  - id: ac-006
    summary: No throughput regression; suites green
    type: performance
    evaluator_hint: operator-action — GH200 single-GPU + CPU node
    pass_when: |
      Single-GPU 1M PE step/s and CPU per-core LJ throughput within tolerance (>= baseline*(1-tol))
      of the recorded AV pre-seam baseline; full CPU + CUDA ctest green.
    status: pending
---

# Acceptance criteria

- **ac-001** — generic, potential-neutral, multi-contribution PairKernel concept (code).
- **ac-002** — existing LJ AND Coulomb-real expressed as kernels (code).
- **ac-003** — shared launchers; PairLJCut + CoulLong ported; base asymmetry unified (code).
- **ac-004** — combined LJ+coul == separate runs (superposition) (scientific; operator).
- **ac-005** — migrated LJ + Coul bit-identical to baseline (scientific; operator).
- **ac-006** — no throughput regression; suites green (performance; operator).
