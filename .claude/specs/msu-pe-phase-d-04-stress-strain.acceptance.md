---
slug: msu-pe-phase-d-04-stress-strain
criteria:
  - id: ac-001
    summary: stress-strain and combined check entry points exist in validate script
    type: code
    pass_when: |
      tools/validate_polyethylene.py accepts a --check argument whose
      choices include "stress-strain" and "all", and defines a
      run_checks dispatch function routing each value to its check;
      `python tools/validate_polyethylene.py --help` lists both choices.
    status: pending
  - id: ac-002
    summary: three-descriptor stress-strain agreement gates within plus/minus 5%
    type: runtime
    pass_when: |
      With link 01's stress_strain_reference.csv present, running
      `python tools/validate_polyethylene.py --check stress-strain
      --atomiverse-dir <run> --reference-dir <ref>` exits 0 when
      elastic modulus, yield stress, and hardening slope each agree
      within +/-5% relative, and exits non-zero (1) when any one
      descriptor breaches +/-5%; verified on a synthetic +6% breach
      fixture and a synthetic +4% pass fixture.
    status: pending
  - id: ac-003
    summary: missing stress-strain reference yields required-input exit code 2
    type: runtime
    pass_when: |
      `python tools/validate_polyethylene.py --check stress-strain`
      run against a reference directory with no stress_strain_reference.csv
      exits with code 2 (not 0, not 1).
    status: pending
  - id: ac-004
    summary: combined Phase D gate runs D.2, D.3, D.4 in one invocation
    type: runtime
    pass_when: |
      `python tools/validate_polyethylene.py --check all --ranks 1` and
      the same with `--ranks 8` each dispatch the thermo, equilibrium,
      and stress-strain checks and exit 0 only when all three pass for
      that rank count; exit 1 if any sub-check breaches tolerance.
    status: pending
  - id: ac-005
    summary: Phase D status row in msu-pe-validation.md marked complete
    type: code
    pass_when: |
      The Phase D row of the status-snapshot table in
      .claude/specs/msu-pe-validation.md states Phase D is complete
      (D.2, D.3, D.4 pass at 1 rank and 8 ranks), replacing the prior
      "ready to run" wording.
    status: pending
---

# Acceptance criteria

- **ac-001** (code) — The CLI surface exists: `--check` exposes `thermo`, `equilibrium`, `stress-strain`, and `all`; `run_checks` dispatches to the corresponding check functions. Verified by reading the script and `--help` output.
- **ac-002** (runtime) — The load-bearing gate. Once chain link 01's `stress_strain_reference.csv` is in the tree, `--check stress-strain` is fully machine-verifiable: it exits 0 on ±5%-agreeing descriptors and non-zero on any breach. Pinned with synthetic +6%-breach and +4%-pass fixtures so it is testable before the real reference CSVs land.
- **ac-003** (runtime) — A missing stress-strain reference is treated as an operator error (exit 2), not a silent skip, so the Phase D.4 gate cannot pass vacuously.
- **ac-004** (runtime) — `--check all` is the Phase D exit-gate entry point: it runs D.2 + D.3 + D.4 together and is invoked once at `--ranks 1` and once at `--ranks 8`. Exit 0 at both rank counts is the Phase D pass condition.
- **ac-005** (code) — On landing, the operator updates the Phase D status row of `msu-pe-validation.md` to mark Phase D complete. Grep-verifiable on the status-snapshot table, hence `type: code`.
