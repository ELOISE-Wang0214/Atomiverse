---
slug: phase1-productization-03-denser-grid
criteria:
  - id: ac-001
    summary: New Lebedev orders 38, 86, 770 integrate spherical harmonics exactly
    type: runtime
    pass_when: |
      Extended test_grid.cpp Lebedev-exactness test integrates every
      spherical harmonic up to each grid's design algebraic degree for
      orders 6, 38, 50, 86, 110, 194, 302, 590, 770 and the result
      matches the analytic value to <= 1e-12; the test binary passes.
    status: pending
  - id: ac-002
    summary: SG-1 grid recovers water electron count to <= 1e-4
    type: runtime
    pass_when: |
      In test_grid.cpp, the SG-1 grid-integrated electron count
      sum(w_i * rho(r_i)) for H2O at def2-SVP equals N_electrons (10)
      with absolute error <= 1e-4; the test passes.
    status: pending
  - id: ac-003
    summary: match-PySCF grid recovers water electron count to <= 1e-4
    type: runtime
    pass_when: |
      In test_grid.cpp, the match-PySCF grid-integrated electron count
      for H2O at def2-SVP equals N_electrons (10) with absolute error
      <= 1e-4; the test passes.
    status: pending
  - id: ac-004
    summary: SG-1 brings def2-SVP energy residue on H2O/NH3/CH4 to <= 5e-4 Ha
    type: runtime
    pass_when: |
      In test_grid.cpp, the ATV single-point KS/PBE energy with the
      SG-1 grid differs from the stored PySCF def2-SVP reference by
      <= 5e-4 Ha for each of H2O, NH3 and CH4; the test passes.
    status: pending
  - id: ac-005
    summary: match-PySCF grid brings the same residue to <= 1e-6 Ha
    type: runtime
    pass_when: |
      In test_grid.cpp, the ATV single-point KS/PBE energy with the
      match-PySCF grid differs from the stored PySCF def2-SVP reference
      by <= 1e-6 Ha for each of H2O, NH3 and CH4; the test passes.
    status: pending
  - id: ac-006
    summary: Becke partition weights sum to unity at every grid point
    type: runtime
    pass_when: |
      In test_grid.cpp, for a multi-atom molecule (H2O) the sum of
      per-atom Becke cell weights at each grid point equals 1 within
      1e-12; the test passes.
    status: pending
  - id: ac-007
    summary: default_config and scf_solve.h doc comment updated to the new grid
    type: code
    pass_when: |
      include/atomiverse/cpu/scf/scf_solve.h: SCFConfig has a
      `grid_level` field, `default_config` sets it to
      GridLevel::MatchPySCF, and the former "SG-0-ish" doc comment is
      replaced with text describing the match-PySCF level-3 default.
    status: pending
  - id: ac-008
    summary: bench_qc_suite energy threshold tightened from 1e-2 Ha
    type: code
    pass_when: |
      benchmarks/qc/bench_qc_suite.cpp: `threshold_ha` is no longer
      1e-2 and is set to a value <= 5e-4 Ha, with an accompanying
      comment naming the grid mode that justifies it.
    status: pending
  - id: ac-009
    summary: GridLevel::SG0 legacy path reproduces prior (50,110) behavior
    type: runtime
    pass_when: |
      In test_grid.cpp, building a grid with GridLevel::SG0 produces a
      MolecularGrid whose point count and weights match the pre-change
      (n_radial=50, n_angular=110) uniform-order grid; the test passes.
    status: pending
  - id: ac-010
    summary: Full check and ctest suite stay green
    type: runtime
    pass_when: |
      `cmake --build build -j` succeeds and `ctest --test-dir build
      --output-on-failure` reports all tests passed, including the
      extended test_grid.cpp.
    status: pending
---

# Acceptance criteria

- **ac-001** — Lebedev exactness for the three newly tabulated orders is the
  GREEN check for the table-addition task; degrees verified are 38→degree 11+
  (per the SG-1 order set), 86→degree 17+, 770→degree 47, alongside the
  existing orders.
- **ac-002 / ac-003** — Electron-count recovery is the cheapest physically
  meaningful grid-quality probe; checked independently for both new modes so a
  regression in either is isolated.
- **ac-004 / ac-005** — The honesty-constrained energy targets: SG-1 is held to
  its ~320 microHartree design accuracy (rounded to a 5e-4 Ha test bar),
  match-PySCF to the sub-microHartree regime (1e-6 Ha). No "1 microHartree with
  SG-1" criterion is written, per the verified domain reality check.
- **ac-006** — Guards the Becke-partition audit task; a unity-sum failure is the
  signature of an incorrect switching-function or size-adjustment term.
- **ac-007 / ac-008** — Static-source checks for the config and benchmark
  threshold changes.
- **ac-009** — Ensures the `SG0` branch is a true backward-compatible path so
  legacy regression baselines that pin (50,110) do not silently shift.
- **ac-010** — Whole-suite regression gate.
