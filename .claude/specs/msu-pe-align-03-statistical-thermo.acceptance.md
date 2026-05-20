---
slug: msu-pe-align-03-statistical-thermo
criteria:
  - id: ac-001
    summary: BlockStat and StatTolerances dataclasses replace the per-step types
    type: code
    pass_when: |
      tools/validate_polyethylene.py defines frozen dataclasses
      BlockStat {column,mean,sem,n_blocks} and StatTolerances with
      ~2% defaults on T and E_pot; ThermoTolerances and StepDeviation
      are removed.
    status: verified
    last_checked: 2026-05-20
  - id: ac-002
    summary: block_average partitions a series and drops the trailing partial block
    type: code
    pass_when: |
      tools/validate_polyethylene.py block_average partitions values
      into contiguous block_size blocks, returns mean of block means,
      SEM across blocks, and n_blocks == floor(len/block_size).
    status: verified
    last_checked: 2026-05-20
  - id: ac-003
    summary: compare_thermo_statistical compares block-averaged means per column
    type: code
    pass_when: |
      tools/validate_polyethylene.py compare_thermo_statistical
      block-averages each column for both engines and returns a breach
      per column whose relative mean gap exceeds its StatTolerances
      threshold; compare_thermo_per_step is removed.
    status: verified
    last_checked: 2026-05-20
  - id: ac-004
    summary: block_average unit tests pass (known mean/SEM, partial-block drop)
    type: runtime
    evaluator_hint: pytest tools/test_validate_polyethylene.py
    pass_when: |
      The block_average tests in tools/test_validate_polyethylene.py
      pass: constant series gives that mean and zero SEM, and a
      non-multiple length drops the trailing partial block.
    status: verified
    last_checked: 2026-05-20
  - id: ac-005
    summary: Statistical comparison breach detection and exit codes verified
    type: runtime
    evaluator_hint: pytest tools/test_validate_polyethylene.py
    pass_when: |
      The compare_thermo_statistical and main integration tests pass:
      agreeing fixtures exit 0, a >2% T-mean-gap pair exits 1 naming
      T, and absent reference CSVs exit 2.
    status: verified
    last_checked: 2026-05-20
  - id: ac-006
    summary: Per-column tolerance is applied independently for T and E_pot
    type: runtime
    evaluator_hint: pytest tools/test_validate_polyethylene.py
    pass_when: |
      The per-column threshold test passes: a mean gap breaching T but
      under the E_pot tolerance is reported as a T-only breach.
    status: verified
    last_checked: 2026-05-20
  - id: ac-007
    summary: D.2 spec and acceptance regenerated as the statistical-comparison spec
    type: code
    pass_when: |
      .claude/specs/msu-pe-phase-d-02-thermo-parity.md and
      .acceptance.md no longer mention per-step / 1e-4 comparison and
      instead describe block-average statistical comparison; the
      acceptance criteria restart at ac-001 with status pending.
    status: verified
    last_checked: 2026-05-20
  - id: ac-008
    summary: Real cross-engine block-averaged T and E_pot agree within ~2%
    type: scientific
    pass_when: |
      With link-01 reference CSVs present,
      `validate_polyethylene.py --check thermo` exits 0 for both the
      1-rank and 8-rank pairs, confirming block-averaged <T> and
      <E_pot> agree within ~2%.
    status: pending
---

# Acceptance criteria

- **ac-001..ac-003** (code) — the statistical machinery replaces the per-step machinery in `validate_polyethylene.py`: new dataclasses, `block_average`, `compare_thermo_statistical`; the per-step types are removed.
- **ac-004..ac-006** (runtime) — the pytest suite proves `block_average` arithmetic, breach detection, exit-code semantics, and per-column tolerance independence on synthetic fixtures.
- **ac-007** (code) — the superseded D.2 spec and acceptance files are regenerated as the statistical-comparison spec, grep-verifiable: no per-step / 1e-4 language remains and criteria restart at `ac-001 / pending`.
- **ac-008** (scientific) — the genuine cross-engine parity check; can only run once `msu-pe-phase-d-01-reference-runs` commits the operator reference runs, so it cannot be auto-verified by the implementer.
