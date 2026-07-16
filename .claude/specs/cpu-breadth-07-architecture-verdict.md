---
title: CPU breadth 07 — evidence-based architecture verdict
slug: cpu-breadth-07-architecture-verdict
status: approved
execution: evidence_blocked
last_progress_update: 2026-07-16
science_required: false
created: 2026-07-15
acceptance: cpu-breadth-07-architecture-verdict.acceptance.md
depends_on:
  - cpu-breadth-02-public-minimize-task
  - cpu-breadth-03-charged-npt-kspace
  - cpu-breadth-04-molpy-forcefield-pipeline
  - cpu-breadth-05-remd-component-scope
  - cpu-breadth-06-cv-umbrella
---

# CPU breadth 07 — evidence-based architecture verdict

## Summary

This is a decision link, not a feature link. After the CPU probes close, it
compares their source, lifecycle, performance, and public-consumer evidence and
decides which local mechanisms become durable architecture protocols.

The verdict must distinguish three outcomes for every candidate abstraction:

1. **promote** — at least two independent production callers need the same
   contract and the shared form reduces duplicated semantics;
2. **keep local** — one probe works, but generalization has no second caller;
3. **reject/redesign** — the probe exposed an ownership, lifecycle, or
   performance failure.

It must not use deferred GPU work as evidence for or against CPU closure, and
must not claim cross-backend generality from CPU results.

## Progress (2026-07-16)

**Execution: evidence-blocked by links 04–06 — 0/8 tasks and 0/9 active
acceptance criteria verified. GPU remains deferred.**

Links 01–03 already provide partial architecture evidence, but the verdict is
intentionally not started: no abstraction is promoted until 04–06 supply the
required second production callers and their acceptance gates close.

## Domain basis

No new physical method is implemented here. Scientific evidence is inherited
from the accepted public regressions and external validations in links 02–06.
The decision is valid only if those references, tolerances, and lifecycle gates
are verified; a compile-only or unit-stub result is not evidence.

## Recon / evidence inputs

The verdict consumes, at minimum:

- link 01: exact build/test manifests, cross-repo ABI contract, component
  ownership, harmonic public conformance, and performance baselines;
- link 02: typed Minimize task, Python/MPI parity, and loop result transport;
- link 03: long-range pairing, complete stress, cell refresh/rebuild, NPT;
- link 04: resolved ForceField snapshot, dense type ingress, native compiled
  components, and three-repository public flow;
- link 05: per-replica program scope, shared replica executor, REMD and NEB;
- link 06: CV protocol, requested capabilities, umbrella programs, and molrs
  WHAM ownership;
- the accepted frame-block API and Frame-only currency contract;
- measured build/steady/rebuild allocation, transfer, memory, and timing data.

## Architecture questions to decide

### Component ownership and instance scope

- Is explicit `add_ref` plus owned-by-default registration sufficient?
- Does per-replica program construction remain a replica-local API, or does a
  second non-replica caller justify a general component-instance-scope protocol?
- Is any clone requirement proven? A speculative `ComponentBase::clone()` is
  rejected unless two real non-copyable callers require it.

### Driver/task/application surface

- Did typed task-specific alternatives support MD, Evaluate, and Minimize
  without harmful duplication?
- Is a wider PlanSpec/NodeSpec justified, or should TaskSpec stay a closed
  transport allow-list until a fourth task exists?
- Does task result metadata belong to a generic node result or remain typed?

### Mutable state and evaluator protocols

- Should long-range pairing/cell refresh become a general dependency/
  invalidation protocol, or stay in electrostatics?
- Does requested-capability dispatch have two real evaluators after Bias, or
  should it remain local until AM1/PW/SCF energy-only paths adopt it?
- Are CV and bias correctly located in a cross-workflow `bias/` domain?

### Cross-repo boundary

- Is the molrs immutable resolved snapshot the durable ForceField boundary?
- Is CompiledForceField a reusable non-Component program across three Drivers,
  or did style-specific behavior leak into presets/hot loops?
- Which unsupported styles are the next minimal completeness slice?

### Replica and analysis ownership

- Is ReplicaExecutor demonstrably shared by NEB, REMD, and umbrella without
  method branches?
- Does Frame remain the only simulation state currency?
- Is all WHAM/statistical analysis in molrs with no C++ duplicate?

## Design / deliverables

Produce `docs/developer/cpu-architecture-breadth-verdict.md` with:

1. a capability matrix over the six CPU workflow shapes, serial/MPI,
   C++/Python/launcher/molpy consumption, fixed/variable atoms, fixed/mutable
   cell, one/N replicas, and analysis boundary;
2. one row per architecture question with `promote | keep local |
   reject/redesign`, callers, measured evidence, costs, and follow-up;
3. a lifecycle table for every probe: ownership, build inputs, rebuild trigger,
   steady allocation, transfer/collective, memory scaling, and public result;
4. a list of architecture claims now proven, still unproven, and explicitly
   GPU-deferred;
5. roadmap corrections for implementation/documentation contradictions;
6. exactly one recommended next CPU breadth frontier and its entry criteria.

The candidate next CPU frontiers are evaluated, not all approved:

- AM1 analytic Force -> Evaluate/Minimize/MD/Hessian/Python reuse;
- variable-cell Minimize -> non-3N degrees of freedom;
- MPI variable-N GCMC + long-range electrostatics -> OutcomePolicy/capacity/
  KSpace transaction interaction;
- OPLS torsion/improper/non-trivial 1-4/KSpace ForceField completeness;
- replica-parallel CPU MPI groups.

The choice maximizes a new architecture axis per implementation cost. It may
not be selected only because a method is scientifically fashionable.

Update `.claude/notes/notes.md` with durable decisions and synchronize roadmap/index
facts. After the campaign is accepted, completed spec files may be deleted in
the normal repository convention; git history retains them.

## Files to create or modify

- `docs/developer/cpu-architecture-breadth-verdict.md` (new)
- `docs/developer/architecture.md` and relevant roadmap chapters
- `.claude/notes/notes.md`
- `.claude/specs/INDEX.md`
- no production implementation is created by this link

## Tasks

- [ ] Verify every non-deferred acceptance criterion from links 01–06 and
      collect exact commands/results, not prose assertions.
- [ ] Build the capability and lifecycle matrices from public regressions,
      validation, profile manifests, and performance artifacts.
- [ ] Decide each architecture question with named callers and measured costs;
      do not promote an abstraction with only one production caller.
- [ ] Audit Frame-only state, analysis ownership, task typing, component
      lifetime/scope, and rebuild behavior against standing directives.
- [ ] Reconcile roadmap/spec/docs claims with the code and accepted matrix.
- [ ] Choose one next CPU frontier and document why the alternatives wait.
- [ ] Record deferred GPU questions without converting them to passed or using
      them to block CPU closure.
- [ ] Write durable docs/NOTES, update INDEX, and only then close/delete the
      completed campaign specs according to repository convention.

## Testing strategy

This link adds no algorithm test. Its acceptance script checks ledgers,
artifacts, links, and source facts:

- every active prior AC is `verified` and points at a reproducible command or
  artifact;
- every scientific link has exactly one public ≤60 s regression with
  `REFERENCE.md`;
- performance/lifecycle artifacts identify runner/config/revision;
- the verdict contains no unsupported cross-backend claim;
- docs/roadmap statements match the accepted source matrix.

## Deferred — GPU

All CUDA, multi-GPU, Torch, GPU Python/launcher, and cross-backend questions
remain deferred. The verdict lists them as unproven and recommends a later GPU
campaign only after CPU protocols and acceptance instrumentation are stable.

## Out of scope

- Implementing the selected next CPU frontier.
- Refactoring code merely to make the verdict aesthetically uniform.
- Declaring every local probe mechanism a global interface.
- Any CUDA/Torch production or test change.
