---
title: CPU breadth 02 — Minimize as the second public task shape
slug: cpu-breadth-02-public-minimize-task
status: approved
execution: in_progress
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-02-public-minimize-task.acceptance.md
depends_on: cpu-breadth-01-acceptance-foundation
---

# CPU breadth 02 — Minimize as the second public task shape

## Summary

Expose CPU geometry optimization through the same four consumption layers as
MD and AM1 single point: public C++, the versioned in-memory `TaskSpec`, Python,
and `mpirun atomiverse script.py`. The minimum delivered method is LJ +
NeighborList + L-BFGS + `ForceConverged`.

This is a probe of the application surface, not a request for a universal
serialized component graph. The current task model hard-codes MD phases and a
closed MD component variant. Adding exactly one second staged loop reveals
which fields are genuinely task-generic and which must remain typed per task.

## Domain basis

For a Lennard-Jones pair

`U(r) = 4 epsilon [(sigma/r)^12 - (sigma/r)^6]`,

the minimum is `r* = 2^(1/6) sigma` and `U(r*) = -epsilon`. The public
regression minimizes a displaced two-atom dimer in a box large enough to avoid
interaction with additional periodic images. It checks distance [Å], energy
[eV], maximum force [eV/Å], and the convergence reason. `REFERENCE.md` cites
the analytic derivative and derives the tolerances from the configured force
threshold and L-BFGS line search.

## Progress / current evidence (2026-07-16)

**Execution: in progress — 7/9 tasks; 6/9 active acceptance criteria
verified. GPU Minimize remains deferred.**

- **Public LJ dimer regression shipped** (`reg_lj_dimer`): analytic r*,
  force_converged; ac-004/ac-006 verified.

- TaskSpec is a strict v3 typed sum. `MinimizeTask` owns typed
  `minimizer::Phase` bindings, a positive iteration budget, explicit component
  IDs, LBFGS parameters, and all five `ForceConverged` tolerances. The v2
  decoder and all variant-order identity are absent.
- Codec tests pin field-exact round-trip, cross-kind/phase rejection, invalid
  budgets and tolerances, unknown/non-current versions, truncation, and explicit
  component IDs.
- `TaskMaterializer` builds one native `cpu::minimizer::Minimize` Plan; a unit
  test compares its typed Plan signature and handler configuration with direct
  public construction. Protocol v2 returns `DispatchResult{Frame, TaskResult}`
  and uses the existing owner-gather path on 1/2/4 ranks.
- Python exposes `atomiverse.minimizer.{Minimize, LBFGS, ForceConverged,
  MinimizeResult}` through ordinary `add`/`on` composition and mutates the
  caller's same molrs Frame. `device="cuda"` fails loud.
- Still blocking closure: the public LJ dimer regression, direct-vs-two-rank
  launcher parity, complete convergence/stop/budget result matrix, repeated
  line-search lifecycle/allocation gate, and the combined acceptance run.

## Architecture questions

1. Can TaskSpec add a new typed phase domain without erasing every task into a
   raw integer phase or inventing a user-visible IR?
2. Can the launcher keep one task-kind-agnostic run path while preserving a
   task-specific result (`iterations`, `converged`, `stop_reason`)?
3. Can the same public Minimize configuration run serial and MPI with identical
   final geometry and verdict?
4. Does the Python surface mirror C++ composition rather than grow a separate
   optimization API or hidden defaults?

## Design

### TaskSpec v3 stays a closed, typed sum

Increment the post-hard-cutover wire version from v2 to v3 and add
`TaskKind::Minimize`. Refactor the in-memory
payload into task-specific alternatives rather than widening the existing MD
record:

- `MDTask` retains `MDComponentEntry{optional<md::Binding>, ...}`;
- `SinglePointTask` retains the AM1 evaluator payload;
- `MinimizeTask` owns `MinimizeComponentEntry` whose explicit binding contains
  `ATV::minimizer::Phase`, an iteration budget, L-BFGS parameters, and the five
  `ForceConverged` tolerances;
- the initial Minimize allow-list reuses PairLJCut and NeighborList parameter
  value types and adds LBFGS. CG and other evaluators remain unrepresentable
  until their Python constructors exist;
- phase encoding still uses only `slot_of` and validates with the matching
  `slot_count<Phase>()`; there is no raw user-authored phase integer.

The codec uses the explicit `ComponentId` mechanism established by frame-block
04; variant order never becomes a wire identity. It rejects a component or
binding from the wrong task kind and names both domains. The v2 codec is
replaced, not retained behind a decoder; no multi-version decoder, version
negotiation, compatibility branch, or translation is permitted. All ranks in
one job run the exact same binary and every non-v3 payload fails before any
task field is materialized.

### Materializer and task result

`TaskMaterializer` adds a `cpu::minimizer::Minimize` alternative and replays
ordinary registrations. It builds once in the constructor and runs through the
existing kind-agnostic `run(frame)` entry.

After execution it exposes a small control result containing iterations,
`converged`, and `stop_reason`. This is loop/control metadata, not a second
simulation-state currency: positions, forces, scalars, box, and topology still
travel only in Frame. The MPI protocol returns the result to rank 0 and gathers
the minimized owner rows through the same domain-decomposition release path as
MD. A stopped/non-converged run is distinguishable from a converged run and is
never reported as success merely because the budget ended.

### Python surface

Add canonical deep paths:

- `atomiverse.minimizer.Minimize`
- `atomiverse.minimizer.LBFGS`
- `atomiverse.driver.MinimizeResult` (or an equivalently named immutable result)

`Minimize.run(frame, max_iter)` mutates the same molrs Frame and returns the
result. Registration uses the same `add`/`on` model as MD; there is no Python
function that accepts parallel x/y/z arrays. `device="cpu"` is the only active
device and any other device fails loud.

### Public and launcher regressions

`regressions/minimizer/lj_dimer` is the one public literature regression. A
separate launcher parity case runs the same Python source directly and with two
MPI ranks, then compares final tag-ordered coordinates, energy, iteration count,
convergence flag, and reason. It does not compare two hidden C++ implementations.

## Files to create or modify

- `include/atomiverse/cpu/task/task_spec.h`, task codec/materializer headers
- `src/cpu/task/` implementations
- task/protocol/worker result plumbing under `cpu/task/` and `app/atomiverse/`
- `python/src/bindings/bind_minimizer.cpp` (new) and binding registration
- `python/src/atomiverse/minimizer/` (new) and public exports/docs
- unit tests for TaskSpec v3, materialization, wrong-domain binding, and results
- `regressions/minimizer/lj_dimer/` (new)
- `regressions/launcher/minimize_lj_parity/` (new)
- CMake registration and `.claude/notes/notes.md`

## Tasks

- [x] Write failing codec tests for Minimize round-trip, wrong-phase-domain
      rejection, invalid budgets/tolerances, unknown version, and truncated
      result payload.
- [ ] Write the public LJ dimer regression and serial-vs-two-rank launcher
      parity case before implementing the new task kind.
- [x] Add typed task-specific payload alternatives and bump the exact wire
      version to v3; keep all phase conversion behind `slot_of` and all
      component encoding behind explicit `ComponentId`.
- [x] Add LBFGS and ForceConverged parameter records matching the actual C++
      constructors field-for-field.
- [x] Materialize and build `cpu::minimizer::Minimize` once, return explicit
      loop result metadata, and gather the final distributed Frame.
- [x] Bind the minimal namespaced Python surface with `device="cpu"` and exact
      Frame-in/Frame-out semantics.
- [ ] Verify build/launch lifecycle: repeated line-search force replay does not
      rebuild components or allocate steady workspaces.
- [ ] Run unit, public regression, Python, serial launcher, two-rank launcher,
      and controlled CPU performance/allocation gates.
- [x] Record what remained task-specific and what became shared in
      `.claude/notes/notes.md`; do not generalize beyond two real task kinds.

## Testing strategy

- **Unit**: codec field equality, typed phase rejection, materializer Plan
  equivalence, result status, non-converged budget exhaustion, and MPI result
  reduction. No full optimization is hidden in default ctest.
- **Regression**: one two-atom LJ public optimization against analytic `r*`
  and `-epsilon`, ≤60 seconds.
- **Launcher acceptance**: direct Python versus two MPI ranks, tag-ordered
  Frame and result parity.
- **Performance**: build and line-search/steady launch reported separately;
  steady LBFGS workspaces allocate zero after build on the controlled runner.

## Deferred — GPU

- `cuda::minimizer::Minimize` Python binding and TaskSpec materialization.
- Device-resident result/status transport and CUDA launcher dispatch.
- CPU/CUDA minimized-geometry and performance parity.

The `device` argument fails loud for deferred devices. No CUDA or Torch
production source is modified in this link.

## Out of scope

- SCF/GPW/PW evaluators in the Minimize wire allow-list.
- Variable-cell optimization, FIRE, transition-state methods, or constraints.
- A generic plugin registry, reflection-based constructor serialization, or
  on-disk Plan format.
- Any decoder, converter, or migration tool for historic TaskSpec versions.
