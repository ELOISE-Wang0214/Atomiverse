---
slug: multigpu-bonded-01-reverse-halo
criteria:
  - id: ac-001
    summary: reverse_comm + gather/scatter kernels + HaloPlan reverse buffers exist
    type: code
    pass_when: |
      MultiGpu::reverse_comm() is declared/defined, gather_ghost_forces and
      scatter_add_forces kernels exist, HaloPlan has dst-side reverse buffers
      allocated in build_halo and freed in destroy_halo_plan.
    status: verified
    last_checked: 2026-06-06
  - id: ac-002
    summary: GH200 build green + LJ multi-GPU tests unchanged
    type: performance
    evaluator_hint: operator-action GH200 (verified — job 9926)
    pass_when: |
      build-cuda builds the multi_gpu.cu change on GH200 and test_cuda_multi_gpu
      stays 13/13 (reverse_comm not on the LJ path → no regression).
    status: verified
    last_checked: 2026-06-06
---

# Acceptance criteria

- **ac-001** — reverse_comm + kernels + HaloPlan reverse buffers (code, verified).
- **ac-002** — GH200 build green + test_cuda_multi_gpu 13/13 (performance, verified — job 9926, commit bdc1645).
