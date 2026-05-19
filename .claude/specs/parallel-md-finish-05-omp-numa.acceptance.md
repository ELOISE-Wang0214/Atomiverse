---
slug: parallel-md-finish-05-omp-numa
criteria:
  - id: ac-001
    summary: README has Thread binding subsection with two numerical measurements OR explicit no-effect note
    type: docs
    pass_when: |
      `benchmarks/md/reference/polyethylene/README.md`
      contains a subsection whose heading matches regex
      `(?i)^#+\s*thread\s*binding`. The subsection body
      contains EITHER two numerical wall times (e.g.,
      `real \d+\.\d+` lines, one for unset, one for
      `OMP_PROC_BIND=close`) AND a delta percentage, OR
      the exact phrase "single-socket — no NUMA effect
      observed" with the host CPU brand string and the
      reasoning. (Closes parent ac-019.)
    status: pending

  - id: ac-002
    summary: Measurement commands are reproducible
    type: docs
    pass_when: |
      The Thread binding subsection includes both literal
      shell-runnable command lines (with `OMP_NUM_THREADS=4`
      and `OMP_PROC_BIND=close` settings) so a future reader
      can re-measure on a different host. The command lines
      reference `./build/benchmarks/md/pe_profile_single`.
    status: pending
---

# Acceptance criteria

- **ac-001** closes parent ac-019. The parent text already permits the "single-socket — no NUMA effect observed" escape; this acceptance encodes both branches (real numbers OR documented no-effect) so the measurement on Apple Silicon doesn't get blocked by the absence of NUMA.
- **ac-002** ensures the section is *useful* to a future reader: not just a snapshot but a recipe.
