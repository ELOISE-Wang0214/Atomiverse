---
slug: driver-device-unification-02-backend-driver-merge
criteria:
  - id: ac-001
    summary: backend_traits exports pair_evaluator_base / thermo_type / virial_stress_type
    type: code
    pass_when: |
      Both `backend_traits<cpu::device_tag>` and `<cuda::device_tag>` export
      `pair_evaluator_base`, `thermo_type`, `virial_stress_type` typedefs to the
      backend's concrete types. Project builds clean.
    status: pending
  - id: ac-002
    summary: BackendDriver<P,DT> exists and replaces CPUDriver + CUDADriver
    type: code
    pass_when: |
      `include/atomiverse/facade/backend_driver.h` defines
      `ATV::BackendDriver<PhaseEnum, DeviceTag>` carrying build / setup_pass /
      step_once / teardown_pass (via stage_in/stage_out/drain) / outer_status /
      for_each_pending / throw_missing / install_defaults / scan_pending_infra_
      + the 8 shared categorical setters. `include/atomiverse/cpu/cpu_driver.h`
      and `include/atomiverse/cuda/cuda_driver.cuh` no longer exist.
    status: pending
  - id: ac-003
    summary: All 7 category classes reparented onto BackendDriver
    type: code
    pass_when: |
      `cpu::{Evaluate,Dynamics,Optimize,Response}Driver` and
      `cuda::{Evaluate,Dynamics,Optimize}Driver` each inherit
      `ATV::BackendDriver<TheirPhase, DeviceTag>` (grep of the class headers).
      No reference to CPUDriver/CUDADriver remains.
    status: pending
  - id: ac-004
    summary: cuda::md::MD no longer redefines shared setters or install_infra_
    type: code
    pass_when: |
      `cuda::md::MD` defines only Dynamics-specific set_integrator/
      set_thermostat/set_barostat locally; set_pair/bond/angle/dihedral/
      neighbor/observer/record/thermo and the auto-injection pass are inherited
      from BackendDriver. `grep -n 'install_infra_' include/atomiverse/cuda/md/md.cuh`
      returns 0.
    status: pending
  - id: ac-005
    summary: Auto-injection is typed-query only — no name strings, ≤1 pair dynamic_cast
    type: code
    pass_when: |
      `grep -rn 'name() == "virial_stress"\|name() == "neighbor_list"'
      include/atomiverse` returns 0 lines. The shared scan consults reads() /
      capabilities() / needs_neighbor_list() / provides_neighbor_list(); any
      dynamic_cast to a pair-evaluator base appears at most once (cutoff query).
    status: pending
  - id: ac-006
    summary: New BackendDriver test + full regression pass
    type: runtime
    pass_when: |
      `ctest --test-dir build -R test_backend_driver` exits 0 (shared-setter
      Plan build + typed VirialStress injection), and the full
      `ctest --test-dir build --output-on-failure` exits 0 with no behavior
      change vs link 01 (same Plans, same validation throws).
    status: pending
---
