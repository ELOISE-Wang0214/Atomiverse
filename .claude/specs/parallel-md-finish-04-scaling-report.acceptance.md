---
slug: parallel-md-finish-04-scaling-report
criteria:
  - id: ac-001
    summary: scaling_report.md exists with build-configuration block
    type: docs
    pass_when: |
      `benchmarks/md/reference/polyethylene/scaling_report.md`
      exists. Its body contains a subsection whose heading
      matches regex `(?i)build configuration|compiler` and
      includes the literal substring `clang` (or `g++` /
      `gcc-` depending on host) AND `-O3` AND `-march=native`.
      (Closes parent ac-022.)
    status: pending

  - id: ac-002
    summary: LJ kernel speedup row recorded with current and baseline numbers
    type: performance
    pass_when: |
      `scaling_report.md` contains a row referencing
      `LJ::launch` with three numeric columns: a baseline
      (scalar / pre-B.2) figure quoted from the parent
      spec status table, a current (xsimd) figure pulled
      from `benchmarks/md/reference/polyethylene/profile_1.csv`,
      and a ratio column. The recorded ratio is ≥ 1.2× OR
      the report carries a one-paragraph note explaining
      why measured ratio fell below the historical 1.21×
      claim (regression callout). (Closes parent ac-011,
      refined: uses profile_1.csv LJ::launch row in lieu
      of the never-built CPU `bench_md_step --force-only=lj`.)
    status: pending

  - id: ac-003
    summary: Bonded kernel speedup row recorded
    type: performance
    pass_when: |
      `scaling_report.md` contains a row referencing
      `HarmonicBond::launch` with baseline (pre-B.3, quoted
      as 1.13× from parent spec status table), current
      (from profile_1.csv), and ratio. Ratio ≥ 1.1× OR a
      regression callout. (Closes parent ac-013, refined
      to use profile_1.csv row.)
    status: pending

  - id: ac-004
    summary: Neighbor-list rebuild cadence comparison present
    type: performance
    pass_when: |
      `scaling_report.md` contains a row referencing
      `NeighborList::build_half` showing both an Atomiverse
      rebuild count (from `profile_1.csv` `n_calls` column)
      and a LAMMPS rebuild count (from
      `runs/lammps_n1/log.lammps` "Neighbor list builds"
      summary line, or a documented derivation if the
      summary is absent). The Atomiverse count is annotated
      either ≤ LAMMPS count (target met) or with a one-line
      analysis. (Closes parent ac-015.)
    status: pending

  - id: ac-005
    summary: Single-core LAMMPS ratio row + structural-gap escape hatch
    type: performance
    pass_when: |
      `scaling_report.md` contains a row whose first
      column matches regex `single.core.*PE.tension` with
      both Atomiverse and LAMMPS wall-time values and a
      computed ratio. EITHER the ratio is ≤ 1.2× AND that
      number is clearly displayed, OR the report carries
      a subsection whose heading matches
      `(?i)structural.?gap` and contains ≥ 80 words of
      prose naming the remaining bottleneck and citing
      its `profile_1.csv` row. (Closes parent ac-023.)
    status: pending

  - id: ac-006
    summary: 8-rank LAMMPS ratio row ≤ 2.0× OR documented analysis
    type: performance
    pass_when: |
      `scaling_report.md` contains a row whose first column
      matches regex `8.rank.*PE.tension` with both Atomiverse
      and LAMMPS wall-time values and a computed ratio. The
      ratio is ≤ 2.0× OR — given that parent ac-024 has no
      built-in escape hatch — a dated note from the user
      records explicit approval to defer the strict 2.0×
      gate to Phase 6, in which case the row still reports
      both numbers and the measured ratio. (Refines parent
      ac-024: strict pass when ≤ 2.0×, soft-pass when user
      has signed off on deferring.)
    status: pending

  - id: ac-007
    summary: Hotspot summary lists top-3 components at 1 and 8 ranks
    type: docs
    pass_when: |
      `scaling_report.md` has a subsection whose heading
      matches `(?i)hotspot` and lists the top-3 components
      by `pct_of_step` at each of 1 rank and 8 ranks,
      sourced from `profile_1.csv` / `profile_8.csv`.
    status: pending
---

# Acceptance criteria

- **ac-001** closes parent ac-022 by locating the compiler+flags annotation in the report rather than alongside bench numbers (the CPU `bench_md_step` binary never existed).
- **ac-002 / ac-003** are the *refined* versions of parent ac-011 / ac-013. The refinement points the criteria at `pe_profile_single`'s per-component rows, which already produce the LJ and bonded kernel timings the original criteria wanted. The historical 1.21× / 1.13× baselines from the parent status table are taken as authoritative; the ratio check uses the committed profile_1.csv as the "current" half.
- **ac-004** closes parent ac-015 with explicit pointers to the artifacts that hold both counts.
- **ac-005** closes parent ac-023 — strict ≤ 1.2× OR the structural-gap-analysis escape hatch the parent acceptance already declares.
- **ac-006** is the strictest item in the chain: parent ac-024 demands ratio ≤ 2.0× without an escape hatch. This refinement permits a *user-signed* deferral if the measured ratio exceeds 2.0×, recognizing that perfect single-machine LAMMPS parity is a long-term goal whose closure may require structural work outside Phase 5.
- **ac-007** is the hotspot writeup the parent B.1 acceptance demanded; it's small but locks in that the profile CSVs were actually analyzed, not just committed.
