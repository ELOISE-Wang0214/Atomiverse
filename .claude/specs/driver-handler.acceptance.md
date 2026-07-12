---
slug: driver-handler
criteria:
  - id: ac-001
    summary: Termination is a Handler, and it is reusable across every workflow
    type: code
    pass_when: |
      LoopHandler exists and is NOT a template — it sees LoopState (iteration, frame,
      context, metrics) and NEVER a phase, which is exactly what lets one handler serve
      MD, Minimize, MC and Path alike. Verdict is {Continue, Converged, Stop}. The
      composition rule Stop > Converged > Continue is implemented and unit-tested. A
      driver accepts SEVERAL handlers.
      A user-written handler (e.g. StopWhenTemperatureExceeds) compiles from OUTSIDE the
      library and halts a run — the thing that is impossible today, because "stop or not"
      is a private `if` inside a driver class.
    status: pending

  - id: ac-002
    summary: The convergence verdict is taken out of LBFGS
    type: code
    pass_when: |
      LBFGS no longer writes ctx.outer_status.converged. It PUBLISHES metrics (max_force,
      rms_force, max/rms displacement, energy change) into LoopMetrics; a
      ForceConverged handler reads them and returns the Verdict. Changing the convergence
      criteria does not touch the optimizer. The existing Minimize tests pass unchanged
      through the new path (behaviour-preserving refactor).
    status: pending

  - id: ac-003
    summary: Status belongs to the loop, not to the Context — the nested-stop bug dies
    type: code
    pass_when: |
      Each loop node owns its own IterationStatus. A unit test builds a driver whose
      inner loop reports Converged and asserts the OUTER loop keeps running — the bug
      that shipped twice (probe 2: LBFGS's unconditional converged write; probe 3: the
      first converged NEB image would have stopped the whole chain) is now structurally
      unwritable, not a rule to remember.
    status: pending

  - id: ac-004
    summary: No int stage anywhere in an API a human writes
    type: code
    pass_when: |
      Backend::run_phase(int) is gone, replaced by a stage-typed entry point; the
      `run_phase(static_cast<int>(p))` calls in path_driver.h are gone with it. A grep for
      `static_cast<int>(` applied to a Phase enum outside the single internal slot_of()
      conversion returns nothing. DefaultBinding still carries an int + PhaseKind — that
      pair IS the erasure (ComponentBase is not a template, so default_bindings() must
      return one type) and it stays constructible ONLY through make_binding(PhaseEnum),
      so no caller ever types an int. The cross-workflow misbind gate still throws.
    status: pending

  - id: ac-005
    summary: A Driver is an executable node — nesting works
    type: code
    pass_when: |
      Driver<PhaseEnum> is itself an executable node (build/launch/destroy, like a
      component), so a Plan slot holds a component OR a sub-driver. Type erasure happens
      at the NODE boundary — a sub-driver keeps its own phase ENUM, distinct from its
      parent's. A unit test nests a driver in a driver, on ONE shared Runtime (no second
      Backend, no rebuilt domain decomposition), and asserts the inner loop's build()
      ran once at the outer build() and allocates nothing per outer iteration.
    status: pending

  - id: ac-006
    summary: MC's accept/reject is an OutcomePolicy, not a termination
    type: code
    pass_when: |
      Commit/rollback after the stages is an OutcomePolicy (Metropolis is one
      implementation), orthogonal to LoopHandler. SamplingDriver as a class is gone; MC
      is a preset carrying TrialBudget + Metropolis.
    status: pending

  - id: ac-007
    summary: The six category bases are DELETED and replaced by presets
    type: code
    pass_when: |
      cpu::EvaluateDriver / DynamicsDriver / OptimizeDriver / ResponseDriver /
      SamplingDriver / PathDriver no longer exist, nor does CPUDriver as a category layer.
      MD / Minimize / MC / Evaluate / Hessian / NEB survive as PRESETS: classes deriving
      Driver<Phase> that configure their default handler + infra auto-injection IN THEIR
      CONSTRUCTOR. They are classes with constructors, NOT make_md() factory functions
      (CLAUDE.md hard rule). Inheritance depth drops from four levels to two. CLAUDE.md's
      "Driver taxonomy" section is rewritten: category is no longer a concept, and the
      two-axis table from probe 3 is folded away (replica count is a nested loop, not an
      axis). The CUDA mirrors follow.
    status: pending

  - id: ac-008
    summary: The regressions pass UNCHANGED — not one line edited
    type: scientific
    pass_when: |
      All five existing regression cases (md/lj_argon, md/eam_cu, md/sw_si,
      md/rigid_water, path/neb_muller_brown) still pass with ZERO edits to their source.
      They are out-of-library consumers of the public API, so if the presets were done
      right they must not be able to tell that the architecture underneath changed. This
      is the strongest acceptance signal in this spec: a refactor that requires editing
      its own end-to-end consumers has changed the API, whatever it claims.
    status: pending

  - id: ac-009
    summary: The verdict is written down and CLAUDE.md no longer teaches categories
    type: code
    pass_when: |
      .claude/NOTES.md records why category was deleted (three probes, one wall: a
      termination semantic welded into a class hierarchy cannot be reused, composed, or
      extended without a new class each time); that Runtime/Driver/Handler is the
      replacement triple; why type erasure sits at the node boundary and not the stage
      boundary; and the ONE thing this spec does not solve — per-instance component state
      under nesting, i.e. the component-clone problem probe 3 found, which nesting makes
      SHARPER rather than dissolving.
    status: pending
---

# Acceptance — driver-handler

The refactor three probes argued for. Its deliverable is a deletion: `category`.

Criteria are written as **requirements**, not mechanisms — the lesson of probes 1–3,
where four criteria across three specs had to be amended or withdrawn because they named
a solution that turned out to be the wrong guess.
