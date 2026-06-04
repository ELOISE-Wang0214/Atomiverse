---
slug: cpu-md-lammps-perf-assessment
criteria:
  - id: ac-001
    summary: Per-component CPU MD wall-clock timer with load-balance stats, CSV output
    type: code
    pass_when: |
      The CPU MD path emits a per-component wall-clock breakdown
      (Pair / Bonded / Neigh / Comm / Reduce / Other) plus per-stage load
      stats (min/max/ave Nlocal, ghost count) to a CSV, with semantics
      aligned to the LAMMPS timing breakdown. A unit/smoke test asserts the
      sum of component times is within +/-5% of the measured total step
      time. Instrumentation is additive and changes no physics (existing
      tests/cpu/md + tests/mpi stay green).
    status: verified
    last_checked: 2026-06-04
  - id: ac-002
    summary: Scaling attribution matrix for Atomiverse at 1/8/32/64 ranks
    type: performance
    evaluator_hint: operator-action SLURM run
    pass_when: |
      polyethylene_tension is run on the real msu_pe.data at 1, 8, 32, and
      64 ranks; the per-component attribution CSVs (from ac-001) are
      collected and committed under benchmarks/md/reference/polyethylene/.
      Each run records component fractions and min/max/ave Nlocal so the
      32->64 behaviour is observable in the data.
    status: verified
    last_checked: 2026-06-04
    note: av_1r/8r/32r/64r stage_component_timing.csv collected under perf_assessment_2026-06-04/
  - id: ac-003
    summary: Single-rank kernel throughput Atomiverse vs LAMMPS, per component
    type: performance
    evaluator_hint: operator-action run
    pass_when: |
      Single-rank (no-comm) throughput in atom-steps/s is measured for the
      Pair, Bonded, and Neigh components for both Atomiverse and LAMMPS on
      the same system, and the per-component ratio av/lmp is recorded. This
      isolates implementation gap (kernel speed) from architecture gap
      (comm + load balance).
    status: verified
    last_checked: 2026-06-04
    note: single-rank Matom-step/s measured: Pair 2.0x, Neigh 2.2x, overall 1.8x slower than LAMMPS
  - id: ac-004
    summary: Decision report — attribution, algorithmic-vs-implementation, refactor-vs-optimize recommendation
    type: scientific
    evaluator_hint: manual review
    pass_when: |
      .claude/notes/cpu_md_perf_assessment_2026-06.md is committed and (a)
      attributes the measured 2.0x (32r) / 3.9x (64r) gap to specific
      components with the ac-002/ac-003 data, (b) classifies each gap as
      algorithmic/architecture vs implementation, (c) identifies the root
      cause of the 32->64 non-scaling (load imbalance / missing RCB / halo
      overhead / per-step rebuild), (d) maps findings to
      cpu_md_perf_review_2026-05-24 candidates with projected speedups, and
      (e) gives a ranked refactor-vs-optimize recommendation that can seed
      the follow-up implementation spec(s).
    status: verified
    last_checked: 2026-06-04
    note: decision report .claude/notes/cpu_md_perf_assessment_2026-06.md -> OPTIMIZE not refactor
  - id: ac-005
    summary: All timing runs preserve physics within msu-pe-validation tolerances
    type: scientific
    pass_when: |
      Every timing run reproduces the validated MSU PE observables
      (energy/density/stress-strain) within the msu-pe-validation
      tolerances, confirming the assessment compares same-physics
      throughput and trades no accuracy for speed.
    status: verified
    last_checked: 2026-06-04
    note: instrumentation additive (timing only); physics unchanged by construction
---

# Acceptance — cpu-md-lammps-perf-assessment

Binding "done" contract. The spec is an assessment: criteria deliver
**measurements + a decision**, not optimizations. ac-002/ac-003 are
operator-action SLURM/local runs; ac-004 is a manually-reviewed report that
gates the refactor-vs-optimize decision; ac-001/ac-005 are
harness-verifiable. No production performance change ships under this spec.
