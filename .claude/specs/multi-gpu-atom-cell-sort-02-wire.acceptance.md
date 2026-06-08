---
slug: multi-gpu-atom-cell-sort-02-wire
criteria:
  - id: ac-001
    summary: Multi-GPU forces stay bit-identical to single-GPU with sort on
    type: code
    pass_when: |
      All 19 tests in tests/cuda/md/test_multi_gpu.cu pass on a GH200,
      including FullPotentialParityVsSingleGpu, BondForceParityVsSingleGpu,
      NveTrajectoryParityNlist, NveMigrationParity.
    status: verified
    last_checked: 2026-06-08
  - id: ac-002
    summary: MD and skin-gating suites unregressed
    type: code
    pass_when: |
      test_cuda_md (9) and test_cuda_neighborlist_skin_gating (2) pass on
      GH200 after the sort is wired into reneighbor().
    status: verified
    last_checked: 2026-06-08
  - id: ac-003
    summary: Sort precedes build_halo / tag-map / bonded rebuild in reneighbor
    type: code
    pass_when: |
      The permute call appears in MultiGpu::reneighbor() after migrate() and
      before build_halo(), build_tag_maps(), and rebuild_bonded_tuples();
      confirmed in src/cuda/md/multi_gpu.cu.
    status: verified
    last_checked: 2026-06-08
---
