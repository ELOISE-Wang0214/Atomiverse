---
slug: msu-pe-phase-d-02-thermo-parity
criteria:
  - id: ac-001
    summary: validate_polyethylene.py exposes a thermo-check entry point
    type: code
    pass_when: |
      tools/validate_polyethylene.py defines compare_thermo_per_step
      and run_thermo_check, and main accepts a --check argument
      whose choices include "thermo".
    status: verified
    last_checked: 2026-05-20
  - id: ac-002
    summary: D.2 tolerances are the thermo-check defaults
    type: code
    pass_when: |
      ThermoTolerances default fields set T, E_pot, E_kin relative
      thresholds to 1e-4 and Pxx, Pyy, Pzz relative thresholds to
      1e-3.
    status: verified
    last_checked: 2026-05-20
  - id: ac-003
    summary: per-step compare returns worst breach with step number
    type: runtime
    pass_when: |
      pytest tools/test_validate_polyethylene.py passes; the
      breach-detection test confirms compare_thermo_per_step on a
      fixture deviating at one known step returns a StepDeviation
      naming that step and the correct rel_err.
    status: verified
    last_checked: 2026-05-20
  - id: ac-004
    summary: per-column threshold bands applied independently
    type: runtime
    pass_when: |
      pytest tools/test_validate_polyethylene.py passes the
      threshold-separation test: a ~5e-4 deviation passes on Pxx
      but fails on T.
    status: verified
    last_checked: 2026-05-20
  - id: ac-005
    summary: step misalignment is a failure, not a skip
    type: runtime
    pass_when: |
      pytest tools/test_validate_polyethylene.py passes the
      step-misalignment test: a missing/extra step within the first
      1000 steps produces a reported failure.
    status: verified
    last_checked: 2026-05-20
  - id: ac-006
    summary: thermo check exit codes are correct on fixtures
    type: runtime
    pass_when: |
      main(["--check","thermo", ...]) returns 0 on matching
      fixtures, 1 on breaching fixtures, and 2 on absent reference
      CSVs, as asserted by the integration test in
      tools/test_validate_polyethylene.py.
    status: verified
    last_checked: 2026-05-20
  - id: ac-007
    summary: cross-engine per-step thermo parity over first 1000 steps
    type: runtime
    evaluator_hint: bench harness — blocked until msu-pe-phase-d-01 commits reference CSVs
    pass_when: |
      with the four reference CSVs from msu-pe-phase-d-01 committed
      under benchmarks/md/reference/polyethylene/, running
      validate_polyethylene.py --check thermo exits 0 for both the
      1-rank and 8-rank pairs.
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 (code)** — Static structure of the script: the thermo-check API exists and carries the D.2 tolerances as defaults (no `--strict` gate). Verified by reading the file.
- **ac-003 / ac-004 / ac-005 (runtime)** — Unit-level behavior of `compare_thermo_per_step`, exercised by `pytest tools/test_validate_polyethylene.py` against synthetic fixtures: breach localization, per-column threshold separation, and step-alignment enforcement. Not blocked by the chain dependency.
- **ac-006 (runtime)** — Integration behavior of `main --check thermo`: exit codes 0/1/2 on fixture inputs. Not blocked.
- **ac-007 (runtime)** — The load-bearing science gate: real cross-engine parity over the first 1000 steps for both rank counts. Machine-verifiable (`validate_polyethylene.py --check thermo` exit code) but **blocked until chain link `msu-pe-phase-d-01-reference-runs` commits the four reference CSVs**. This is a recorded chain dependency; `status` stays `pending` until those artifacts land and the check is run.
