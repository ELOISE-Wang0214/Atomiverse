---
slug: cpu-breadth-07-architecture-verdict
criteria:
  - id: ac-001
    summary: all prior active CPU criteria are verified reproducibly
    type: runtime
    pass_when: |
      An acceptance audit finds no pending/failed non-deferred criterion in
      cpu-breadth links 01-06. Each verified runtime/scientific/performance AC
      records a reproducible command and revision/artifact; deferred GPU ACs
      remain deferred and are excluded from the count.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-002
    summary: CPU capability matrix covers every architecture axis exercised
    type: docs
    pass_when: |
      cpu-architecture-breadth-verdict.md has explicit rows/columns for all six
      CPU workflows, serial/MPI, C++/Python/launcher/molpy, fixed/variable N,
      fixed/mutable cell, one/N replicas, and analysis ownership. Every cell is
      proven, absent, limited, or deferred with evidence; no blank implies yes.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-003
    summary: lifecycle matrix is quantitative rather than aspirational
    type: docs
    pass_when: |
      For links 01-06 the verdict records owner, build inputs, rebuild trigger,
      steady allocation, transfer/collective budget, memory scaling, public
      result, runner/config, and measured baseline. Unknown values are named
      gaps, not omitted.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-004
    summary: every candidate abstraction receives an evidence-based disposition
    type: docs
    pass_when: |
      Component lifetime/scope, task typing/result, long-range pairing/cell
      refresh, requested capabilities, ForceField snapshot/program, CV/bias,
      ReplicaExecutor, and molrs analysis each receive promote/keep-local/
      reject-redesign with at least two named callers for every promotion.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-005
    summary: Frame and cross-repo ownership directives are audited
    type: code
    pass_when: |
      The verdict proves the accepted frame-block API remains the sole
      simulation-state currency in Minimize/REMD/umbrella and proves
      ForceField/WHAM truth remains in molrs with no production hot-loop FFI or
      duplicated analysis math. It does not reopen Frame storage migration.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-006
    summary: architecture claims distinguish proven, limited, and unproven
    type: docs
    pass_when: |
      A dedicated section lists claims now proven by public applications,
      limitations such as replicated-grid PME and single-rank replica methods,
      and unproven axes. It contains no statement that CPU evidence proves CUDA,
      multi-GPU, Torch, or cross-backend generality.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-007
    summary: roadmap and spec records match accepted implementation facts
    type: docs
    pass_when: |
      Roadmap contradictions and stale frame/spec campaign status are corrected;
      links to regressions/benchmarks/docs resolve; .claude/notes/notes.md contains
      the durable decisions and rejected alternatives; INDEX records closure
      only after every active criterion is verified.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-008
    summary: exactly one next CPU breadth frontier is selected
    type: docs
    pass_when: |
      The verdict compares AM1 Force, variable-cell Minimize, MPI variable-N
      long-range GCMC, ForceField completeness, and replica-parallel CPU MPI by
      new architecture axis, cost, blockers, and application value, then names
      exactly one next frontier with entry/exit conditions. It does not start
      implementation in this decision link.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-009
    summary: decision link makes no production-code change
    type: code
    pass_when: |
      The link's diff is limited to docs, roadmap, .claude/notes/notes.md, specs/index,
      and acceptance evidence metadata. No include/atomiverse, src, python
      binding, molrs production, or molpy production file is changed merely to
      force a preferred verdict.
    status: verified
    last_checked: 2026-07-16
    notes: |
      docs/developer/cpu-architecture-breadth-verdict.md records dispositions from 01–06 evidence. Full reproducible audit of every prior AC is progressive; GPU excluded.

  - id: ac-010
    summary: GPU architecture verdict
    type: docs
    pass_when: |
      DEFERRED. A future GPU campaign must separately judge CUDA/Torch Frame,
      Driver/Plan, synchronization, component scope, charged KSpace/NPT,
      ForceField materialization, CV/bias, and multi-GPU execution. The CPU
      verdict records these as unproven, not failed and not passed.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-07-architecture-verdict

This decision closes only after evidence exists. It cannot waive a failed prior
CPU criterion or convert deferred GPU work into implied parity.
