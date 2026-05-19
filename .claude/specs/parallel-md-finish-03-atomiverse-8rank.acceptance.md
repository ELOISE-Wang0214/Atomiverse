---
slug: parallel-md-finish-03-atomiverse-8rank
criteria:
  - id: ac-001
    summary: 8-rank polyethylene_tension full-mode run exits 0
    type: runtime
    pass_when: |
      `mpirun -n 8 ./build/benchmarks/md/polyethylene_tension
      --mode full --out-dir benchmarks/md/reference/polyethylene/runs/atomiverse_n8`
      exits with status 0 (verified post-hoc by the
      `exit_code: 0` line in
      `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/RUN.md`).
      (Closes parent ac-004 exit-code portion.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-002
    summary: thermo.csv + stress_strain.csv land non-empty
    type: code
    pass_when: |
      Both files exist under
      `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/`
      with non-zero size and ≥ 11 lines each. The header
      of `thermo.csv` matches the schema documented at the
      top of `polyethylene_tension.cpp`
      (`step,stage,T_target_K,T_K,P_bar,...`). The header
      of `stress_strain.csv` matches
      (`step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa`).
      (Closes parent ac-004 CSV portion.)
    status: verified
    last_checked: 2026-05-19

  - id: ac-003
    summary: Resident memory growth ≤ 5% across two consecutive runs
    type: scientific
    pass_when: |
      `RUN.md` records two peak-RSS values (one per
      consecutive run); `(rss_2 - rss_1) / rss_1 ≤ 0.05`.
      This proves the no-allocation-in-launch architectural
      rule (parent §B.2) holds in steady state.
      (Closes parent ac-004 memory portion.)
    status: pending
    observed_2026-05-19: run 1 = 19,333,120 B; run 2 = 19,611,648 B; growth = 1.44% (≤ 5% gate). (Left at pending per /mol:impl rule that only code|runtime auto-flip.)

  - id: ac-004
    summary: RUN.md fully populated with host fingerprint and run metadata
    type: docs
    pass_when: |
      `RUN.md` contains all of: exact `mpirun` command,
      `git rev-parse HEAD` SHA, `uname` output, CPU
      brand string, OMP env vars used, wall time, peak
      RSS (× 2 runs), exit code, and final per-step
      wall time printed by the binary.
    status: pending
    observed_2026-05-19: RUN.md contains exact mpirun command, git SHA ee80e0e, uname output, Apple M4 brand string, OMP env state (unset, with binary's warning), both runs' wall + RSS, exit code 0, and a "Known anomaly" note on the suspicious T_K/E_kin_eV columns. (Left at pending per /mol:impl rule that only code|runtime auto-flip.)
---

# Acceptance criteria

- **ac-001** is the "exit 0" portion of parent ac-004. Run is recorded post-hoc in `RUN.md` so the verification is static (no need to re-run during acceptance).
- **ac-002** verifies both reference CSVs land with the right schema. The schemas are pulled verbatim from the binary's docstring; deviation here means the binary changed.
- **ac-003** is the load-bearing architectural gate: the "no allocation in launch()" rule (parent spec §B.2) is verified by comparing peak RSS across two back-to-back runs. The 5% tolerance allows for OS-level page-cache jitter while still catching genuine leaks.
- **ac-004** documents the run so sub-spec 04 (scaling-report) and downstream readers can reproduce the conditions.

Together these close parent ac-004 in full.
