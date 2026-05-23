---
slug: msu-pe-cuda-components-03-cuda-integrators-stage1
criteria:
  - id: ac-001
    summary: cuda::md::NVELimit caps per-step displacement at xmax
    type: code
    pass_when: |
      tests/cuda/md/test_nve_limit.cu's cap-engages test sets one atom
      with f chosen so the half-step kick gives |v|·dt = 10·xmax; after
      a single launch with stage = initial_integrate, the atom's
      |Δposition| over that step equals xmax to 1e-5 absolute (single
      precision); GTest reports the test passing.
    status: pending
  - id: ac-002
    summary: cuda::md::NVELimit survives huge-force initial conditions without NaN
    type: code
    pass_when: |
      tests/cuda/md/test_nve_limit.cu runs 100 steps on a 16-atom
      configuration with fx[0] = 1e8 eV/Å (msu_pe.data-style overlap);
      every position and velocity component on every atom is finite at
      every recorded step; GTest reports the test passing.
    status: pending
  - id: ac-003
    summary: cuda::md::NVELimit matches cpu::md::NVELimit on a no-force drift
    type: scientific
    pass_when: |
      A 64-atom periodic configuration with random velocities and zero
      forces is propagated 50 steps by both cuda::md::NVELimit{dt=1.0,
      xmax=0.05} and cpu::md::NVELimit{dt=1.0, xmax=0.05}; max
      per-component position / velocity difference < 1e-5 in single
      precision; GTest reports the cross-backend comparison test
      passing.
    status: pending
  - id: ac-004
    summary: cuda::md::FixLangevin reaches the target temperature within ±10%
    type: scientific
    pass_when: |
      tests/cuda/md/test_fix_langevin.cu runs 5000 steps on a 256-atom
      synthetic system at fixed seed with T_start = T_stop = 300 K, no
      LJ, dt = 1.0 fs, damp = 10.0 fs; the kinetic temperature averaged
      over the last 2500 steps is within ±30 K of 300 K; GTest reports
      the test passing.
    status: pending
  - id: ac-005
    summary: cuda::md::FixLangevin honours the linear T_start→T_stop ramp
    type: code
    pass_when: |
      With T_start = 200, T_stop = 500, set_run_span(0, 1000), a single
      launch at ctx.step = 500 leaves comp.T_target() == 350.0 ± 1e-9;
      another launch at ctx.step = 2000 leaves comp.T_target() ==
      500.0 ± 1e-9 (clamp at nlast); GTest reports the ramp test passing.
    status: pending
  - id: ac-006
    summary: Both components advertise the correct contract triple
    type: code
    pass_when: |
      cuda::md::NVELimit::role() == ComponentRole::Integrator and its
      default_bindings() returns exactly {(initial_integrate, 0),
      (final_integrate, 0)}; cuda::md::FixLangevin::role() ==
      ComponentRole::Thermostat, capabilities() == Capability::Force,
      default_bindings() == {(post_force, +500)}; GTest contract test
      passes.
    status: pending
  - id: ac-007
    summary: A Langevin-only MD driver passes validate_components()
    type: code
    pass_when: |
      A cuda::md::MD with only NVELimit + FixLangevin added (no LJ, no
      bonded evaluator) builds successfully — validate_components()
      sees ≥ 1 Integrator role and ≥ 1 Force-capable component (from
      FixLangevin); MD::build() does not throw; GTest test passes.
    status: pending
  - id: ac-008
    summary: New CUDA sources and tests are wired and the suite is green
    type: runtime
    pass_when: |
      cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON succeeds;
      cmake --build build -j compiles the new .cu sources without
      warnings escalated to errors; ctest --test-dir build
      --output-on-failure runs both test_cuda_md_nve_limit and
      test_cuda_md_fix_langevin and the entire suite is green
      (0 FAILED).
    status: pending
  - id: ac-009
    summary: Doxygen on every new public symbol with LAMMPS audit table
    type: docs
    pass_when: |
      Every public symbol in nve_limit.cuh and fix_langevin.cuh has
      @brief, @param with units, @return; both headers carry a LAMMPS-
      correspondence audit table matching the CPU header structure
      (Aspect | LAMMPS | This class | Verdict) covering at minimum
      xmax meaning / cap norm / rescale formula / both halves for
      NVELimit and argument order / damp meaning / ramp formula / RNG
      divergence note for FixLangevin; visual inspection on the PR
      confirms presence.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** pin `NVELimit`: the cap is exact, huge forces don't NaN, and the CUDA port matches the already-LAMMPS-aligned CPU implementation on a closed-form no-force drift.
- **ac-004 / ac-005** pin `FixLangevin` physics: temperature setpoint convergence (FDT-consistent kick + drag) and the linear temperature ramp.
- **ac-006 / ac-007** pin the structural contract — `role()`, `capabilities()`, `default_bindings()` — and verify that the `Capability::Force` advertisement on `FixLangevin` actually lets a Langevin-only `cuda::md::MD` build.
- **ac-008** is the standard build + ctest green gate.
- **ac-009** enforces the project's doxygen + LAMMPS-correspondence-table convention already used by `cpu/md/nve_limit.h` and `cpu/md/fix_langevin.h`.
