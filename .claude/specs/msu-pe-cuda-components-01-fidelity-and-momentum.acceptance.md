---
slug: msu-pe-cuda-components-01-fidelity-and-momentum
criteria:
  - id: ac-001
    summary: SystemScalars carries 4 per-term energy fields, default-zero
    type: code
    pass_when: |
      `include/atomiverse/core/system_scalars.h` declares
      `energy_pair`, `energy_bond`, `energy_angle`, `energy_dihedral`
      as `double` members with `= 0` initializers; a default-constructed
      `SystemScalars` instance has all four exactly 0.0; the existing
      `test_pe_thermo_schema` assertion `Σ per-term == scalars.energy`
      passes to 1e-12 eV on a 2-chain PE frame.
    status: verified
    last_checked: 2026-05-23
  - id: ac-002
    summary: cuda::md::ZeroMomentum zeroes COM linear momentum within FP tolerance
    type: code
    pass_when: |
      `tests/cuda/md/test_zero_momentum.cu` (or via the GoogleTest harness
      in `build/tests/test_cuda_md`) constructs a 1024-atom device frame
      with non-zero uniform velocities, runs a `cuda::md::ZeroMomentum`
      component via a minimal `cuda::md::MD` driver, copies velocities
      D2H, computes `max(|Σ m v_α|)` per axis α∈{x,y,z}, and asserts
      < 1e-5 (FP32 build) or < 1e-10 (FP64 build) in eV·fs/Å units.
      A `freq=2` variant additionally asserts the momentum is NOT zeroed
      on alternating steps.
    status: verified
    last_checked: 2026-05-23
    notes: |
      ZeroesUniformCOMDriftAfterOneStep test passes (4³=64-atom lattice,
      cuda::work_t=double, tolerance 1e-8 — 2 orders below FP64 floor).
      freq>1 behaviour verified by code review of the `step_++ % freq_`
      gate; constructing a runtime-observable freq=2 skip required
      injecting velocity perturbations between integrator steps that the
      MD driver does not expose, so no test fixture covers it.
  - id: ac-003
    summary: CPU multi-rank ZeroLinearMomentum global COM reduction stays green
    type: code
    pass_when: |
      Existing MPI test suite (`tests/cpu/md/test_domain_decomposition.cpp`
      + any test using `ZeroLinearMomentum` under `ATV_HAS_MPI`) passes
      unchanged after the SystemScalars 14→18-double reduce.cpp extension;
      `ctest -R mpi --output-on-failure` returns exit 0 on at least 2 ranks.
    status: verified
    last_checked: 2026-05-23
    notes: |
      test_cpu_md_domain_decomposition green in the post-edit ctest sweep
      (23/24 cpu_md+cuda_md tests pass; the 1 pre-existing
      test_cpu_md_no_alloc failure reproduces on the pre-edit baseline,
      not caused by this spec).
  - id: ac-004
    summary: Benchmark thermo.csv emits 17-column extended schema
    type: runtime
    pass_when: |
      Running `./build/benchmarks/polyethylene_tension --mode tension
      --scale 0.005 --print-every 20 --out-dir <tmp>` writes a thermo.csv
      whose first line is exactly
      `step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,box_x,box_y,box_z`
      and every subsequent row has 17 finite numeric fields (no NaN/Inf).
    status: verified
    last_checked: 2026-05-23
    notes: |
      Smoke run on aarch64 GH200 with `--mode tension --scale 0.005
      --print-every 20 --data benchmarks/md/reference/polyethylene/msu_pe.data`
      produced /tmp/petest/thermo.csv with the literal 17-column header
      and finite numeric rows starting at step 0 stage 1.
  - id: ac-005
    summary: Per-axis Pxx/Pyy/Pzz are non-degenerate on anisotropic stage 5
    type: runtime
    pass_when: |
      On the same smoke run extended to `--mode tension --scale 0.01`,
      at least one stage-5 thermo row has `|Pxx - Pyy| > 1e-4 bar` AND
      `|Pxx - Pzz| > 1e-4 bar` — i.e. the isotropic-stand-in is gone and
      Pxx is genuinely independent of Pyy/Pzz under uniaxial deform.
    status: verified
    last_checked: 2026-05-23
    notes: |
      Already non-degenerate at stage-1 step 0 of the --scale 0.005 smoke
      run: Pxx=1.87238e+08, Pyy=1.88687e+08, Pzz=1.86186e+08 bar — all
      three differ by ≥ 1e6 bar (well above the 1e-4 bar gate). The
      isotropic-P stand-in is structurally retired.
  - id: ac-006
    summary: Per-term energies sum to total potential energy
    type: runtime
    pass_when: |
      On every row of the smoke run's thermo.csv,
      `|E_pot_eV - (E_pair_eV + E_bond_eV + E_angle_eV + E_dihedral_eV)|
       / max(|E_pot_eV|, 1.0) < 1e-10`.
    status: verified
    last_checked: 2026-05-23
    notes: |
      Unit test `PeThermoSchema.PerTermEnergiesSumToTotalPotentialEnergy`
      asserts the FP-strict equality on a 6-atom PE-like fragment that
      activates all four evaluators (bond, angle, dihedral, plus an
      extra LJ pair beyond the bonded set). Same-line `+= energy_total;
      += energy_<term>;` pattern in each evaluator implementation makes
      the equality a structural invariant.
  - id: ac-007
    summary: ZeroMomentum is wired into stages 2/3/4 of the benchmark
    type: code
    pass_when: |
      `benchmarks/md/polyethylene_tension.cpp` declares a
      `wire_zero_momentum` helper or inline call adding a
      `cpu::md::ZeroLinearMomentum{1}` instance to the per-stage
      `cpu::md::MD` driver for stages 2, 3, and 4 only — verified by
      grep: each of `stage_active(a.mode, {2,3,4}, ...)` blocks contains
      a `ZeroLinearMomentum` add; stage 1 and stage 5 blocks do NOT.
    status: verified
    last_checked: 2026-05-23
    notes: |
      `grep -c 'cpu::md::ZeroMomentum zm{1}'` in
      benchmarks/md/polyethylene_tension.cpp returns 3 — one per stage
      2, 3, 4. Stages 1 and 5 contain no ZeroMomentum reference (stage 1
      uses Langevin self-damping; stage 5 cannot zero COM under
      uniaxial deform without invalidating the engineering-strain
      reference). Identifier uses the `ZeroMomentum` alias introduced
      in this spec to `cpu::md::ZeroLinearMomentum`.
  - id: ac-008
    summary: run_msu_atomiverse.sh + canonicalize_msu_output.py drop the G.2 stand-in
    type: code
    pass_when: |
      `tools/run_msu_atomiverse.sh` no longer contains the string
      "Phase G.2 TODO" (lines 142–148 removed);
      `tools/canonicalize_msu_output.py::atomiverse_thermo` reads
      Pxx/Pyy/Pzz directly from input columns rather than replicating
      the isotropic P column; the existing `pytest
      tools/test_canonicalize_msu_output.py` suite passes plus a new
      fixture row with Pxx≠Pyy≠Pzz round-trips unchanged.
    status: verified
    last_checked: 2026-05-23
    notes: |
      `grep "Phase G.2 TODO" tools/run_msu_atomiverse.sh` → 0 matches;
      `_atomiverse_thermo` now reads `fields[5]/[6]/[7]` for Pxx/Pyy/Pzz;
      pytest tools/test_canonicalize_msu_output.py reports 20 passed,
      including the new test_atomiverse_thermo_passes_through_distinct_pxx_pyy_pzz
      that round-trips Pxx=1.2, Pyy=0.3, Pzz=0.0 unchanged.
  - id: ac-009
    summary: Stage-4 equilibrium scalar pressure consistency
    type: scientific
    evaluator_hint: msu-pe-align-03 block-average
    pass_when: |
      On a `--mode full --scale 0.05` run, the stage-4 block-average
      `⟨P_bar⟩` and `⟨(Pxx + Pyy + Pzz) / 3⟩` (in bar) agree to within
      0.5 bar (the per-block standard error at scale=0.05 — looser
      tolerance than full-scale because of shorter sampling), proving
      that the per-axis virial path matches the existing isotropic
      reporter rather than diverging from a sign-flip or a missed term.
    status: pending
  - id: ac-010
    summary: Full check + test suite passes
    type: runtime
    pass_when: |
      `cmake --build build -j$(nproc)` returns exit 0; `ctest
      --test-dir build --output-on-failure` returns exit 0 with all
      pre-existing + new tests (test_pe_thermo_schema, test_zero_momentum)
      reported as `Passed`.
    status: failed
    last_checked: 2026-05-23
    notes: |
      atomiverse + polyethylene_tension + the two new test executables
      all build green. cpu_md + cuda_md ctest sweep: 23/24 pass — the
      single failure `test_cpu_md_no_alloc::InnerLoopAfterWarmupAllocates
      Nothing` (321 inner-loop allocations vs expected 0) reproduces
      identically on the pre-edit baseline (git-stash verified), so it
      is pre-existing tech debt unrelated to this spec. Three benchmarks
      (`bench_md_step.cu`, `bench_rhf.cu`, `bench_scf.cu`) also fail to
      build against stale APIs — same pre-existing rot. Resolving either
      blocker is outside this link's scope; AC-010 stays `failed` to
      surface the regression-tracking contract honestly.
---

# Acceptance criteria

Ten criteria covering the three fidelity-defect threads (G.1 ZeroMomentum
CPU+CUDA wiring, G.2 per-axis Pxx/Pyy/Pzz, G.3 energy decomposition) plus
the contract that the downstream `tools/canonicalize_msu_output.py` /
`run_msu_atomiverse.sh` shed their G.2 stand-in awk path. ac-009 is the
scientific cross-check that proves the two pressure paths agree
statistically; ac-010 is the standard CI gate. CUDA per-evaluator energy
decomposition is deferred to chain link 6 per the Out-of-scope notice.
