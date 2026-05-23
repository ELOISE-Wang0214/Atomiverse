---
slug: msu-pe-cuda-components-04-cuda-nose-hoover-npt
criteria:
  - id: ac-001
    summary: cuda::md::NoseHooverNPT header declares Config + role()==Integrator + default_bindings + set_run_span
    type: code
    pass_when: |
      include/atomiverse/cuda/md/nose_hoover_npt.cuh exists, declares
      class ATV::cuda::md::NoseHooverNPT : public cuda::ComponentBase
      with public Config aggregate (T_target, tau_T, P_target, tau_P, dt,
      drag=0.0, mtchain=3, mpchain=3, couple_axis[3]={true,true,true},
      T_stop=-1.0), explicit Config ctor, set_run_span(int64, int64)
      returning NoseHooverNPT&, role() returning ComponentRole::Integrator,
      and default_bindings() returning entries for
      initial_integrate@0 + final_integrate@0 + end_of_step@-40.
    status: pending
  - id: ac-002
    summary: build/destroy lifecycle manages all device allocations cleanly
    type: code
    pass_when: |
      src/cuda/md/nose_hoover_npt.cu's build() cudaMallocs d_eta_t,
      d_eta_t_dot, d_eta_b, d_eta_b_dot (sized mtchain / mpchain),
      d_eps_dot[3], d_ke_global, d_ke_axis[3], d_diag_axis[3], d_scale[3],
      d_block_count; sets eps_dot[a]=0 for !couple_axis[a]; initializes
      Q_t / Q_b per LAMMPS MTK masses; sets tdrag = 1 − drag·dt/tau_T,
      pdrag = 1 − drag·dt/tau_P; destroy() cudaFrees each allocation and
      leaves State zeroed. Verified by AddressSanitizer-clean build +
      cuda-memcheck on the test suite reporting zero leaks/errors.
    status: pending
  - id: ac-003
    summary: launch() reproduces CPU Trotter ordering across the three stages
    type: code
    pass_when: |
      src/cuda/md/nose_hoover_npt.cu's launch() invokes, in order:
      initial_integrate → nhc_baro_half, nhc_atoms_half, update_eps_dot,
      scale_v_by_baro, kick_v_by_force, propagate_x_and_box (P1–P6);
      final_integrate → kick_v_by_force only (P7);
      end_of_step → recompute KE_global, scale_v_by_baro, update_eps_dot,
      nhc_atoms_half, nhc_baro_half, write frame.scalars->kinetic and
      frame.scalars->conserved (P8–P11). Order verified by code review
      against cpu/md/nose_hoover_npt.cpp::launch.
    status: pending
  - id: ac-004
    summary: Isotropic NPT relaxes a perturbed LJ box to P_target within tolerance
    type: runtime
    pass_when: |
      tests/cuda/md/test_cuda_nose_hoover_npt.cu test
      `NoseHooverNPT_IsoRelaxesPressure` runs 2000 steps on an 8-atom
      cubic LJ block with couple_axis={true,true,true}, P_target=0,
      tau_P=1000 fs, dt=0.001 fs and asserts block-average ⟨P⟩ over the
      last 500 steps satisfies |⟨P⟩| < 5e-5 eV/Å³ (~80 bar);
      frame.box(0,0) finite and positive; |frame.box(0,0) −
      frame.box(1,1)| < 1e-3 Å.
    status: pending
  - id: ac-005
    summary: Anisotropic y/z-only coupling leaves x rigid while y/z relax
    type: runtime
    pass_when: |
      Test `NoseHooverNPT_AnisoYZOnlyKeepsXRigid` runs the same 8-atom
      LJ block with couple_axis={false,true,true}, x-direction +1%
      pre-strain, 2000 steps. Asserts frame.box(0,0) ==
      box_x_initial (bit-equal, since ε̇_x≡0); ⟨Pyy⟩ and ⟨Pzz⟩ block
      averages over last 500 steps each |·| < 5e-5 eV/Å³;
      ⟨Pxx⟩ separated from ⟨Pyy⟩/⟨Pzz⟩ by ≥ 1e-4 eV/Å³ (confirms x is
      not passively pulled by the y/z barostat).
    status: pending
  - id: ac-006
    summary: Temperature ramp tracks T_start → T_stop linear schedule
    type: runtime
    pass_when: |
      Test `NoseHooverNPT_TemperatureRampTracks500to100` runs 1000 steps
      with T_target=500, T_stop=100, set_run_span(0, 999), tau_T=100 fs,
      sampling ⟨KE⟩ every 100 steps. For each sample step `s`, asserts
      |T_measured(s) − (500 − 400·s/999)| / T_expected(s) < 0.08
      (8% small-system fluctuation envelope).
    status: pending
  - id: ac-007
    summary: drag parameter damps barostat oscillations
    type: runtime
    pass_when: |
      Test `NoseHooverNPT_DragDampsBoxOscillation` runs two identical
      configurations (8-atom LJ, eps_dot seeded to 1e-3 /fs, 500 steps)
      with drag=0 and drag=2 respectively. Compares peak-to-peak
      amplitude of frame.box(0,0) over steps [100, 500]: asserts
      amplitude(drag=2) ≤ 0.5 · amplitude(drag=0).
    status: pending
  - id: ac-008
    summary: CPU↔CUDA block-averaged parity on 8-atom LJ within single-precision noise
    type: scientific
    pass_when: |
      Test `NoseHooverNPT_CPUvsCUDAParity` runs identical 8-atom LJ block
      with cpu::md::NoseHooverNPT and cuda::md::NoseHooverNPT (same
      FF, same dt, same InitVelocity seed, 2000 steps,
      couple_axis={true,true,true}). Block-averaged ⟨T⟩, ⟨P⟩, ⟨box_x⟩
      over the last 1000 steps each satisfy relative error
      |X_cuda − X_cpu| / |X_cpu| < 0.03 (3% single-precision band).
    status: pending
  - id: ac-009
    summary: Extended MTK Hamiltonian drift bounded under weak coupling
    type: scientific
    pass_when: |
      Test `NoseHooverNPT_ConservedHamiltonianFlat` runs 500 steps with
      tau_T=tau_P=1e6 fs (NVE limit). Records frame.scalars->conserved
      at step 0 and step 499; asserts |H(499) − H(0)| / |H(0)| < 5e-3.
    status: pending
  - id: ac-010
    summary: Doxygen and LAMMPS cross-references present on all public symbols
    type: docs
    pass_when: |
      doxygen build over include/atomiverse/cuda/md/nose_hoover_npt.cuh
      emits zero warnings; every public member (Config fields, ctor,
      set_run_span, role, default_bindings, build, launch, destroy)
      carries @brief, units-annotated @param, and @return; the class
      docstring cites Tuckerman 2006 / Shinoda 2004 / LAMMPS
      fix_nh.cpp lines 2541–2569.
    status: pending
  - id: ac-011
    summary: Full check + test suite green on CUDA build
    type: runtime
    pass_when: |
      With -DATV_BUILD_TESTS=ON and CUDA toolkit available,
      `cmake --build build -j$(nproc)` completes with zero errors and
      `ctest --test-dir build --output-on-failure -R cuda_nose_hoover_npt`
      reports all `NoseHooverNPT_*` tests passing; full
      `ctest --test-dir build --output-on-failure` regression-free
      (no prior CUDA test starts failing).
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-003** are static (code-shape) gates ensuring the API surface and Trotter ordering match the CPU reference before any numerical check fires.
- **ac-004 / ac-005** lock the two MSU stage paths the parent validation needs: stages 2–4 (iso) and stage 5 (aniso y/z).
- **ac-006** verifies the stage-3 cooling ramp (`500 → 100 K`) — the spec G.7 line item that the CUDA backend explicitly lacks today.
- **ac-007** covers the `drag` knob (LAMMPS `fix npt … drag 2`); without this the over-compressed amorphous-PE start state rings.
- **ac-008** is the load-bearing parity check link 06 will lean on when comparing CPU vs CUDA per-engine benchmark traces; 3% relative tolerance is the single-precision noise band documented in the parent `msu-pe-validation.md` G.* phase.
- **ac-009** is the MTK extended-Hamiltonian conservation regression — the textbook diagnostic that a symplectic NPT integrator is actually symplectic; the 5×10⁻³ bound is the single-precision-MTK upper limit (CPU double-precision typically holds ≤ 1×10⁻⁴ on the same test).
- **ac-010 / ac-011** are the standard documentation + suite-green gates, matching the per-link chain contract.
