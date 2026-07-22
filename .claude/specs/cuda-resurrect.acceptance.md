---
spec: cuda-resurrect
created: 2026-07-22
criteria:
  - id: ac-001
    summary: cuda no longer references deleted Frame members
    type: code
    pass_when: |
      grep over include/atomiverse/cuda/ + src/cuda/ finds no Frame member access
      via .atoms / .box / .scalars (host cpu::Frame is blocks+meta+simbox; cuda
      reads it through span/meta/simbox).
    status: pending
  - id: ac-002
    summary: cuda no longer references the deleted Real type or core/real.h
    type: code
    pass_when: |
      No `core/real.h` include and no catch-all `Real` type in the cuda trees;
      host reals are StorageReal/AccReal, device reals are real_t/work_t/prec_t.
    status: pending
  - id: ac-003
    summary: cuda uses current units constant names
    type: code
    pass_when: |
      No kFtm2v / kMvv2e / kBoltzmann / kPI / kCoulomb in the cuda trees; unit
      conversions use core/units.h current names (units::force_to_accel, …).
    status: pending
  - id: ac-004
    summary: CUDA configures and builds on the GPU node
    type: runtime
    pass_when: |
      On the GH200 node, cmake -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON (+HF/KSPACE)
      configures and `cmake --build` completes with no errors (libatomiverse +
      cuda test targets link).
    status: pending
  - id: ac-005
    summary: CUDA ctest green on the GPU node; CPU suite still green
    type: runtime
    pass_when: |
      ctest on the GPU-node CUDA build passes; and the login-node CPU suite
      (cmake -DATV_BUILD_CUDA=OFF) still reports 236/236 — proof no CPU code was
      touched. The cuda_no_stale_api_scrub is registered and passes.
    status: pending
out_of_scope:
  - CPU code changes; namespace-to-directory alignment; new CUDA features/kernels;
    torch backend; precision default changes.
---

# Acceptance — cuda-resurrect

Done means: the CUDA backend, dead on `dev`, compiles and passes its tests on the
GH200 GPU node again — achieved purely by mirroring the current CPU Frame /
precision / units API on the CUDA side, with zero CPU changes (the login-node CPU
suite stays 236/236) and a grep gate that keeps the stale API from creeping back.
