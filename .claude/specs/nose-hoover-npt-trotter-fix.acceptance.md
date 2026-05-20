---
slug: nose-hoover-npt-trotter-fix
criteria:
  - id: ac-001
    summary: Trotter NHC ordering matches LAMMPS in initial_integrate
    type: code
    pass_when: |
      In `src/cpu/md/nose_hoover_npt.cpp::launch`, the `initial_integrate`
      branch calls `nhc_baro_half` *before* `nhc_atoms_half` (matching
      LAMMPS `fix_nh.cpp` `initial_integrate`).
    status: verified
    last_checked: 2026-05-20

  - id: ac-002
    summary: Trotter NHC ordering matches LAMMPS in end_of_step
    type: code
    pass_when: |
      In `src/cpu/md/nose_hoover_npt.cpp::launch`, the `end_of_step`
      branch calls `nhc_atoms_half` *before* `nhc_baro_half` (matching
      LAMMPS `fix_nh.cpp` `final_integrate`).
    status: verified
    last_checked: 2026-05-20

  - id: ac-003
    summary: LAMMPS parity test passes within 5 % on ⟨T⟩, ⟨V⟩, ⟨KE⟩
    type: runtime
    pass_when: |
      `tests/cpu/md/test_lammps_npt_parity.cpp` is enabled (renamed
      `DISABLED_` → `Slow_NptMatchesLammpsEnsembleAverages`) and runs in
      the default ctest suite. Against the regenerated stable-state-point
      reference (FCC a=5.3 Å, T=100 K, P=2731 bar; LAMMPS
      `lammps_npt_avgs.txt` ⟨T⟩=100.185 K, ⟨KE⟩=3.30222 eV,
      ⟨V⟩=9533.96 Å³) all three EXPECT_LE assertions hold at the 5 % gate.
      Observed residuals: Δ(T)=0.08 %, Δ(KE)=0.08 %, Δ(V)=0.04 %.
    status: verified
    last_checked: 2026-05-20

  - id: ac-004
    summary: Full ctest suite stays green
    type: runtime
    pass_when: |
      `cmake --build build -j && ctest --test-dir build
      --output-on-failure` exits 0 with 93/93 passing, including the
      now-enabled `LammpsParity.Slow_NptMatchesLammpsEnsembleAverages`
      and the pre-existing `CpuMdTest.NoseHooverNptRelaxesPressure`
      smoke test.
    status: verified
    last_checked: 2026-05-20

  - id: ac-005
    summary: docstring + spec status reflect resolution
    type: docs
    pass_when: |
      `tests/cpu/md/test_lammps_npt_parity.cpp` docstring records the
      resolution: the original P=0 benchmark was physically ill-posed
      (liquid evaporates), and the load-bearing bug was the
      non-continuable `DynamicsDriver::run()` rebuilding components
      every call. The benchmark `in.npt` and `lammps_npt_avgs.txt` are
      regenerated for the stable a=5.3 / P=2731 bar state point.
    status: verified
    last_checked: 2026-05-20

  - id: ac-006
    summary: H_NPT_extended conservation diagnostic (scientific)
    type: scientific
    evaluator_hint: human
    pass_when: |
      A throwaway Phase-C diagnostic run computes the MTK extended
      Hamiltonian
        H = KE_atoms + PE + 0.5·W·ε̇²
            + Σ 0.5·Q_t·η̇_t² + Σ 0.5·Q_b·η̇_b²
            + N_dof·kT·η_t[0] + kT·Σ η_t[k≥1] + kT·Σ η_b
            + p_target·V
      and verifies `|H(t) − H(0)| ≤ 10⁻³ eV` over 10 000 steps on the
      regenerated LJ-argon parity system. Closed by manual attestation
      (`/mol:close --manual`); not automated here because the
      diagnostic harness has no production home until the
      `npt-conservation-gate` follow-up spec lands.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** — static checks for the LAMMPS-matching outer-NHC
  Trotter order (spec Phase A). Done, but *not* the load-bearing fix.
- **ac-003** — the parity gate. Passes only because (a) the benchmark
  moved to a mechanically-stable state point and (b) `DynamicsDriver`
  was made continuable so chunked `run(burn_in); run(sample)…` no
  longer wipes the NHC chain state.
- **ac-004** — regression guard: the continuable-`run()` change touched
  every MD test; all 93 stay green (a `SnapshotPublisher` drain barrier
  was added so deferred teardown does not leave snapshots in flight).
- **ac-005** — documents the outcome for future readers.
- **ac-006** — optional scientific gate; manual attestation, deferred
  to the `npt-conservation-gate` follow-up.
