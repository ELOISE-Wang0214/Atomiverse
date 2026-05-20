---
slug: msu-pe-phase-d-03-equilibrium
criteria:
  - id: ac-001
    summary: validate_polyethylene.py exposes a --check equilibrium entry point
    type: code
    pass_when: |
      tools/validate_polyethylene.py defines a --check argument accepting
      'equilibrium' and a _check_equilibrium function; running
      `validate_polyethylene.py --check equilibrium --help` exits 0.
    status: verified
    last_checked: 2026-05-20
  - id: ac-002
    summary: equilibrium helpers compute RDF first peak and end-to-end R-squared
    type: code
    pass_when: |
      tools/validate_polyethylene.py defines _rdf_first_peak, _end_to_end_r2,
      _read_xyz, EquilibriumStats, and compute_equilibrium_stats; the pytest
      suite in tools/test_validate_polyethylene.py covering them passes.
    status: verified
    last_checked: 2026-05-20
  - id: ac-003
    summary: equilibrium check gates all five D.3 observables and exits non-zero on breach
    type: runtime
    evaluator_hint: pytest tools/test_validate_polyethylene.py
    pass_when: |
      With link-01 reference artifacts present (or seeded fixtures),
      `validate_polyethylene.py --check equilibrium` exits 0 when all five
      observables (density 1%, E_pot 0.5%, RDF peak pos <0.05 A, RDF peak
      height 3%, R^2 10%) are within tolerance and exits 1 naming the
      breached observable when any one is exceeded.
    status: verified
    last_checked: 2026-05-20
  - id: ac-004
    summary: equilibrium tolerances match msu-pe-validation Phase D.3 table
    type: code
    pass_when: |
      Default Tolerances values are eq_density_rel=0.01, eq_epot_rel=0.005,
      eq_rdf_peak_pos_abs=0.05, eq_rdf_peak_height_rel=0.03, eq_r2_rel=0.10,
      matching the Phase D.3 table in msu-pe-validation.md.
    status: verified
    last_checked: 2026-05-20
  - id: ac-005
    summary: equilibrium check runs at 1-rank and 8-rank reference artifacts
    type: runtime
    pass_when: |
      `validate_polyethylene.py --check equilibrium` succeeds against both the
      1-rank and 8-rank stage-4 thermo / RDF / trajectory artifacts produced
      by msu-pe-phase-d-01, reporting per-rank pass/fail.
    status: pending
---

# Acceptance criteria

- **ac-001** — The equilibrium gate must be an independently invocable entry
  point (`--check equilibrium`), not buried in a monolithic `main()`. Verified
  by inspecting the source for the `--check` argument and `_check_equilibrium`
  symbol, plus a `--help` smoke run.
- **ac-002** — The structural and conformational observables (RDF first peak,
  end-to-end R²) require dedicated pure helpers, unit-tested with synthetic
  fixtures so their math is verified independently of any real run.
- **ac-003** — The load-bearing gate: the five-observable agreement is the D.3
  exit criterion. The script must exit non-zero on any breach and name the
  failing observable. Machine-verified once link-01 artifacts (or seeded
  fixtures) exist.
- **ac-004** — The tolerances are the binding scientific contract from
  `msu-pe-validation.md` Phase D.3; they must not drift from that table.
  Grep-verifiable against the default `Tolerances` fields, hence `type: code`.
- **ac-005** — Phase D.3 requires the check at 1 and 8 ranks; both reference
  artifact sets (from link 01) must be processable with per-rank verdicts.
