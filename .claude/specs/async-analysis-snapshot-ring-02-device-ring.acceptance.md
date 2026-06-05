---
slug: async-analysis-snapshot-ring-02-device-ring
criteria:
  - id: ac-001
    summary: DevicePublish exists, binds end_of_step@+10, packs on compute stream into device_ring
    type: code
    pass_when: |
      include/atomiverse/cuda/device_publish.cuh and src/cuda/device_publish.cu
      define cuda::DevicePublish whose launch() records pack_done_event on the
      compute/main stream after a D2D pack into device_ring_[slot] and BEFORE
      enqueuing the defer-stream D2H; default binding is md::Phase::end_of_step
      at order +10.
    status: pending
  - id: ac-002
    summary: Tearing regression proves every delivered device-ring snapshot is single-step + correctly stamped
    type: code
    evaluator_hint: operator-run on GH200 (ATV_BUILD_CUDA=ON)
    pass_when: |
      tests/cuda/async/test_async_record_correctness.cu drives an analytic
      trajectory x_i(step)=f(tag_i,step) with every-step capture and a slowed
      host worker, and asserts for every delivered snapshot that all atoms
      satisfy x_i==f(tag_i,snap.step); test FAILS against the pre-fix live-pointer
      path and PASSES against DevicePublish.
    status: pending
  - id: ac-003
    summary: Thermo scalars->energy race is covered explicitly
    type: code
    evaluator_hint: operator-run on GH200
    pass_when: |
      a case in test_async_record_correctness.cu asserts the thermo energy
      delivered for step n equals the frozen step-n value and is not corrupted
      by the next step's ZeroEnergy, exercising the cuda/recorder/thermo.cuh path.
    status: pending
  - id: ac-004
    summary: AsyncRecordBase device_ring[K] never overwrites an in-flight snapshot
    type: code
    evaluator_hint: operator-run on GH200
    pass_when: |
      a back-to-back capture test asserts that with K ring slots and a slowed
      worker, no snapshot whose D2H has not yet drained is overwritten; the
      single-buffer pre-fix code FAILS this and the device_ring[K] code PASSES.
    status: pending
  - id: ac-005
    summary: DevicePublish/device-ring path adds zero hot-path stream syncs
    type: code
    evaluator_hint: operator-run on GH200
    pass_when: |
      test_hot_path_no_sync.cu reset_sync_count() before the step loop and
      asserts cuda::detail::sync_count()==0 after stepping with DevicePublish
      and at least one ring subscriber active.
    status: pending
  - id: ac-006
    summary: issue_d2h consumers no longer read live device pointers
    type: code
    pass_when: |
      grep over include/atomiverse/cuda/recorder/recorder.cuh, recorder/thermo.cuh,
      observer/rdf.cuh, observer/mulliken.cuh shows no issue_d2h (or its
      replacement) reading ctx.frame.atoms.* , &ctx.frame.scalars->energy,
      pairs-> live fields, or workspace().P directly; reads target the frozen
      device-ring slot instead.
    status: pending
  - id: ac-007
    summary: DevicePublish and device_ring carry Doxygen with units
    type: docs
    pass_when: |
      cuda::DevicePublish, its device_ring_ member, and the pack kernel each
      carry @brief and document units (positions Angstrom, energy eV, step
      dimensionless) per the project Doxygen convention.
    status: pending
  - id: ac-008
    summary: CUDA ctest suite passes on GH200
    type: runtime
    evaluator_hint: operator-action on GH200 node (gpu_required)
    pass_when: |
      on the aarch64 GH200 node with ATV_BUILD_CUDA=ON per CLAUDE.md, cmake
      build succeeds and `ctest --test-dir build-cuda --output-on-failure`
      passes including the cuda;async;gpu_required labelled tests.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-004 / ac-006**: 结构与设计 contract——计算流冻结定序、设备环深度 K、消费源切换为冻结 slot。
- **ac-002 / ac-003**: 撕裂修复的可证伪核心，分别覆盖坐标路径与 thermo 能量竞争；要求 RED-先行（对旧活跃指针路径必须失败）。
- **ac-005**: hot-path 零同步不变量，沿用 `sync_guard.cuh` 计数器。
- **ac-007**: Doxygen + 单位。
- **ac-008**: GH200 上 operator-action 的整体运行门。
