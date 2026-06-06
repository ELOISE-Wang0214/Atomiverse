---
slug: multigpu-bonded-05-optimize
criteria:
  - id: ac-001
    summary: P2P batching + launch bounds + no per-step alloc
    type: code
    pass_when: |
      reverse/forward halo packs xyz/forces into a single P2P buffer; bonded
      kernels carry __launch_bounds__; the re-materialize/tag-map path does no
      per-step cudaMalloc/Free (buffers reused).
    status: pending
  - id: ac-002
    summary: correctness unchanged after optimization
    type: scientific
    evaluator_hint: operator-action GH200, >=2 GPUs
    pass_when: |
      The -02..-04 correctness tests (energy/force parity, ΣF≈0) remain green
      after the optimizations.
    status: pending
  - id: ac-003
    summary: measured multi-GPU PE step-time speedup
    type: performance
    evaluator_hint: operator-action GH200, before/after
    pass_when: |
      Multi-GPU PE step-time improves measurably over the first correct version
      (-04), reported with before/after numbers and the GPU count.
    status: pending
---

# Acceptance criteria

- **ac-001** — P2P batching + launch bounds + no per-step alloc (code).
- **ac-002** — correctness unchanged (scientific; operator GH200).
- **ac-003** — measured multi-GPU PE step-time speedup before/after (performance; operator GH200).
