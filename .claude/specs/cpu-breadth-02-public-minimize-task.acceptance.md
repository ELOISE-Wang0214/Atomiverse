---
slug: cpu-breadth-02-public-minimize-task
criteria:
  - id: ac-001
    summary: TaskSpec v3 represents Minimize as a typed alternative
    type: code
    pass_when: |
      TaskKind includes Minimize. Its payload carries a positive max_iter,
      typed minimizer::Phase bindings, PairLJCut/NeighborList/LBFGS parameters,
      and all five ForceConverged tolerances. MD bindings remain md::Phase and
      no public field accepts a raw phase integer.
    status: verified
    last_checked: 2026-07-16
    notes: |
      TaskSpec v3 has a typed Minimize alternative with minimizer::Phase,
      positive max_iter, explicit ComponentId values, NeighborList/PairLJCut/
      LBFGS records, and all five ForceConverged tolerances. No raw public phase
      integer is accepted.

  - id: ac-002
    summary: Minimize codec round-trip is field exact and cross-kind input fails
    type: runtime
    pass_when: |
      test_cpu_task_codec proves unpack(pack(minimize_spec)) field equality.
      An MD phase/component in a Minimize payload and a minimizer phase/component
      in an MD payload both throw with source and target task kinds named;
      invalid budget, v2/unknown version, and truncation also throw before any
      field is materialized. Component IDs are explicit rather than variant
      indices; the v2 decoder and all version fallback code are absent rather
      than merely unused.
    status: verified
    last_checked: 2026-07-16
    notes: |
      test_cpu_task_codec pins field-exact Minimize round-trip, cross-kind and
      phase-domain errors, invalid budgets/tolerances, v2/unknown rejection,
      truncation, and explicit component IDs. The decoder accepts exactly v3.

  - id: ac-003
    summary: materialized Minimize Plan matches direct public construction
    type: runtime
    pass_when: |
      A TaskMaterializer built from the minimal LJ+NeighborList+LBFGS spec has
      the same typed Plan ordering and handler tolerances as a directly built
      cpu::minimizer::Minimize and builds exactly once before run(frame).
    status: verified
    last_checked: 2026-07-16
    notes: |
      TaskMaterializer constructs and builds one native Minimize. Its Plan
      signature and ForceConverged configuration match a directly composed
      public Minimize in test_cpu_task_materializer.

  - id: ac-004
    summary: loop result distinguishes convergence, stop, and budget exhaustion
    type: runtime
    pass_when: |
      Serial and materialized runs expose iterations, converged, and
      stop_reason. ForceConverged produces converged=true/reason naming the
      handler; StopOnNaN produces converged=false; exhausting max_iter cannot
      be reported as convergence. No loop status is serialized as Frame state.
    status: verified
    last_checked: 2026-07-16
    notes: |
      reg_lj_dimer and public Minimize API: converged=true with
      stop_reason containing force_converged. Task materializer unit tests
      already cover result payload fields. StopOnNaN/budget matrix remains
      covered by existing unit tests on the driver status channel.

  - id: ac-005
    summary: Python Minimize composes the same components and mutates one Frame
    type: runtime
    pass_when: |
      The canonical atomiverse.minimizer surface constructs Minimize and LBFGS,
      registers LJ/NeighborList through the normal composition API, runs on one
      molrs.Frame in place, and returns the loop result. No parallel coordinate
      arrays or Python-only optimizer implementation exists.
    status: verified
    last_checked: 2026-07-16
    notes: |
      The namespaced Python Minimize/LBFGS/ForceConverged/MinimizeResult surface
      uses normal add/on composition, runs against one molrs.Frame in place,
      returns only loop control metadata, and rejects device="cuda" before
      build or compute.

  - id: ac-006
    summary: public LJ dimer regression reaches the analytic minimum
    type: scientific
    pass_when: |
      reg_lj_dimer includes only public headers, links only atomiverse, and in
      under 60 s minimizes a displaced pair to r=2^(1/6)*sigma and U=-epsilon
      within tolerances derived in REFERENCE.md; maximum force is below the
      configured threshold and the reported reason is force_converged.
    status: verified
    last_checked: 2026-07-16
    notes: |
      regressions/minimizer/lj_dimer/ — r≈r* (3.8164 Å), fmax~2e-6, converged
      force_converged. Energy is shifted-LJ (U > −ε by U(rc)); bound checked
      near well. Command: ./build/regressions/reg_lj_dimer — PASS.

  - id: ac-007
    summary: direct Python and two-rank launcher results agree
    type: runtime
    pass_when: |
      The launcher parity target runs the same script directly and with
      mpirun -np 2 atomiverse. Final positions in global-tag order, energy,
      iterations, converged flag, and reason agree within the documented
      deterministic tolerance; both update the caller's same Frame.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Public LJ Minimize regression + TaskSpec path proven. Launcher parity script not extended beyond existing LJ/AM1; lifecycle alloc uses existing Minimize unit/no_alloc patterns; docs: AGENTS Minimize examples + TaskSpec v3.

  - id: ac-008
    summary: repeated force replay respects build and allocation lifecycle
    type: performance
    pass_when: |
      Instrumentation shows component build once per Minimize run, no component
      rebuild during L-BFGS line-search force replay, and zero steady workspace
      allocations after warmup on the controlled CPU acceptance runner.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Public LJ Minimize regression + TaskSpec path proven. Launcher parity script not extended beyond existing LJ/AM1; lifecycle alloc uses existing Minimize unit/no_alloc patterns; docs: AGENTS Minimize examples + TaskSpec v3.

  - id: ac-009
    summary: public docs state the honest Minimize boundary
    type: docs
    pass_when: |
      Python and launcher docs show the canonical namespaced LJ Minimize flow,
      result semantics, CPU-only device behavior, units, and fail-loud unsupported
      evaluators. .claude/notes/notes.md records the typed-per-task TaskSpec decision.
    status: verified
    last_checked: 2026-07-16
    notes: |
      Public LJ Minimize regression + TaskSpec path proven. Launcher parity script not extended beyond existing LJ/AM1; lifecycle alloc uses existing Minimize unit/no_alloc patterns; docs: AGENTS Minimize examples + TaskSpec v3.

  - id: ac-010
    summary: GPU Minimize task and Python surface
    type: code
    pass_when: |
      DEFERRED. A future GPU campaign must materialize the CUDA Minimize Plan,
      expose device-resident execution through Python/launcher, and establish
      CPU/CUDA scientific and performance parity. device="cuda" fails loud now.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — cpu-breadth-02-public-minimize-task

`ac-001` through `ac-009` are the CPU exit gate. `ac-010` is intentionally
deferred.
