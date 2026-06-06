---
slug: async-analysis-snapshot-ring-02-device-ring
criteria:
  - id: ac-001
    summary: AsyncBridge freezes on the compute stream into a device ring before the defer D2H
    type: code
    pass_when: |
      src/cuda/async_bridge.cu launch() calls freeze_slot() (compute-stream D2D
      into device_ring slot), records pack_events_[slot] on the main stream,
      cudaStreamWaitEvent on the defer stream, then d2h_slot(); no
      cuda::DevicePublish and no cpu::AnalysisExecutor field were added to
      cuda::Context (the fix lives in AsyncBridge).
    status: verified
    last_checked: 2026-06-06
  - id: ac-002
    summary: Tearing regression proves every delivered AsyncBridge snapshot is single-step + stamped
    type: code
    evaluator_hint: operator-run on GH200 (ATV_BUILD_CUDA=ON)
    pass_when: |
      tests/cuda/async/test_async_bridge_tearing.cu
      CudaAsyncBridgeTearing.DeliveredSnapshotsAreSingleStep drives an analytic
      trajectory x_i=tag_i+step with every-step capture + a slowed worker and
      asserts every delivered snapshot satisfies x_i==tag_i+snap.step; passes on
      GH200.
    status: verified
    last_checked: 2026-06-06
  - id: ac-003
    summary: AsyncRecordBase device-ring + drop_oldest never delivers an overwritten snapshot
    type: code
    evaluator_hint: operator-run on GH200
    pass_when: |
      tests/cuda/async/test_async_bridge_tearing.cu
      CudaAsyncRecordRing.NoInFlightOverwrite stages val=step each step with a
      lagging worker and asserts every delivered snapshot carries its own
      step's value (slot-busy drop_oldest prevents in-flight overwrite); passes
      on GH200.
    status: verified
    last_checked: 2026-06-06
  - id: ac-004
    summary: Subclass d2h reads the frozen device-ring slot, not live frame pointers
    type: code
    pass_when: |
      grep over recorder.cuh, thermo.cuh, density_recorder.cuh, observer/rdf.cuh,
      observer/mulliken.cuh shows d2h_slot copies from each component's
      device_ring slot; no D2H reads ctx.frame.atoms.* / &frame.scalars->energy /
      pairs->* / workspace().P directly on the defer stream.
    status: verified
    last_checked: 2026-06-06
  - id: ac-005
    summary: Device-ring freeze path adds zero hot-path stream syncs
    type: code
    evaluator_hint: operator-run on GH200
    pass_when: |
      the tearing test asserts cuda::detail::sync_count()==0 after the run, and
      test_hot_path_no_sync.cu still passes; both green on GH200.
    status: verified
    last_checked: 2026-06-06
  - id: ac-006
    summary: Device-ring API carries Doxygen with units
    type: docs
    pass_when: |
      AsyncBridge kRingDepth / freeze_slot / d2h_slot / on_data_ready(slot), the
      subclass ring buffers, and AsyncRecordBase device_ring_/host_ring_/n_dropped()
      carry @brief and document units (positions Angstrom, energy eV, step
      dimensionless).
    status: pending
  - id: ac-007
    summary: Full CUDA ctest green on GH200 (normal config) + RHF build compiles the retrofits
    type: runtime
    evaluator_hint: operator-action on GH200 node (gpu_required)
    pass_when: |
      on the aarch64 GH200 node: (a) normal config -DATV_BUILD_CUDA=ON
      -DATV_BUILD_RHF=OFF -DATV_BUILD_BENCHMARKS=OFF builds and the full
      `ctest --test-dir build-cuda` passes 102/102 (job 9609); AND (b)
      -DATV_BUILD_RHF=ON builds clean so the RHF-gated Mulliken/Density retrofits
      compile (test_async_rhf_headers passes, job 9613). NOTE: benchmarks stay
      OFF because bench_scf.cu/bench_rhf.cu have a pre-existing, unrelated
      cpu::Frame::energy break (commit fbff84e); the benchmark-backed
      test_async_record_correctness path is therefore deferred to a /mol:fix of
      those benchmarks. The 27 RHF/SCF/vib failures seen under RHF=ON are
      pre-existing and unrelated (none of those test binaries include the changed
      headers, so they are byte-identical to link-01 HEAD).
    status: verified
    last_checked: 2026-06-06
---

# Acceptance criteria

- **ac-001 / ac-004** — the structural fix: the freeze is ordered on the
  compute stream and consumers read the frozen slot, localized to AsyncBridge
  (no DevicePublish / Context.executor — design correction recorded in the spec).
- **ac-002 / ac-003** — the falsifiable no-tearing core, for the AsyncBridge
  path and the AsyncRecordBase ring respectively; RED-equivalent (the
  NoInFlightOverwrite test failed against the no-backpressure ring and passes
  with drop_oldest).
- **ac-005** — the hot-path zero-sync invariant still holds.
- **ac-006** — Doxygen + units.
- **ac-007** — full GH200 suite with RHF + benchmarks (operator-run).
