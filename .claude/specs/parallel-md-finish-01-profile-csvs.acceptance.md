---
slug: parallel-md-finish-01-profile-csvs
criteria:
  - id: ac-001
    summary: profile_1.csv committed with required header and ≥ 10 data rows
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/profile_1.csv`
      exists in the working tree, first line is exactly
      `component,total_ms,pct_of_step,n_calls,mean_us`,
      and `wc -l < benchmarks/md/reference/polyethylene/profile_1.csv`
      reports ≥ 11 (header + ≥ 10 data rows). At least one
      data row's `component` field starts with `LJ` and has
      `total_ms > 0`. (Closes parent ac-009 single-rank half.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-002
    summary: profile_8.csv committed and reflects multi-rank MPI components
    type: code
    pass_when: |
      `benchmarks/md/reference/polyethylene/profile_8.csv`
      exists with the same header and ≥ 10 data rows; at least
      one row whose `component` starts with `HaloExchangeForward`
      has `total_ms > 0` AND `n_calls > 0` (proving MPI
      collectives actually ran). (Closes parent ac-009 multi-rank
      half.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-003
    summary: plot_profile.py emits non-empty PNGs from both CSVs
    type: code
    pass_when: |
      Running `python tools/plot_profile.py
      benchmarks/md/reference/polyethylene/profile_1.csv`
      produces `benchmarks/md/reference/polyethylene/profile_1.png`
      with file size > 0; same for `profile_8.csv` → `profile_8.png`.
      Both PNGs are committed. (Closes parent ac-010.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-004
    summary: Component-time accounting is plausible
    type: scientific
    pass_when: |
      Sum of the `pct_of_step` column in `profile_1.csv` and
      in `profile_8.csv` each falls within `[90.0, 105.0]`
      (allows up to 10% un-timed gap, up to 5% overhead
      double-counting from nested timers). Verified by a
      one-line awk on each CSV.
    status: pending
    observed_2026-05-19: profile_1.csv sums to 99.77%, profile_8.csv sums to 99.80% — both within [90.0, 105.0]. (Left at `pending` per /mol:impl rule that only code|runtime criteria auto-flip; user may promote to verified.)
---

# Acceptance criteria

- **ac-001 / ac-002** are pure file-existence + content-shape checks: they prove the profile binary ran end-to-end at both rank counts and produced data with the right schema. ac-001 verifies LJ timed (any non-MPI rank works); ac-002 verifies halo exchange timed (only meaningful at ≥ 2 ranks).
- **ac-003** locks in the plot pipeline so downstream readers can browse the data visually; PNGs are committed because the bench host may not be reproducible.
- **ac-004** is a sanity gate on the timer instrumentation itself — if pct_of_step doesn't sum to ≈ 100%, a scoped timer is missing or double-counted and the downstream scaling-report ratios will be misleading.

This spec closes parent acceptance criteria ac-008 (binary exists with scoped timers — already verified statically at source-tree level), ac-009 (CSVs committed), and ac-010 (plot script works).
