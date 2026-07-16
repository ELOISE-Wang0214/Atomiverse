# CPU Backend Internals

This page is the developer counterpart to the user tutorials. It is
written for someone who has decided to *modify* the CPU backend — add a
new integrator, swap in a different SCF helper, change the convergence
contract — and needs to understand why the layering is what it is. The
[Architecture](architecture.md) page sets up the whole-project picture
across CPU and CUDA; this page zooms into the CPU half.

The chapter assumes you have read the
[CPU overview tutorial](../tutorial/ch00-architecture.md). The user-level
mental model — Frame, Component, Driver — is repeated here only when
the developer-level rules require it.

## Why There Is No Driver Category

There used to be five "category" base classes here — `EvaluateDriver`,
`DynamicsDriver`, `OptimizeDriver`, `ResponseDriver`, `SamplingDriver`
— one per kind of run. They are **deleted**, and so are their headers.
They differed from one another by a single `if`: *when does the loop
stop?* Evaluate ran one pass, Dynamics ran a budget, Optimize stopped
on convergence. Every new method wanted a sixth, and no user could
write a seventh, because that `if` was private to a driver class.

Stripped to the bone, a Driver is three things: **stages** (a typed
phase enum), **iteration** (walk them repeatedly), and **termination**
(when to stop). Only the third was welded into the hierarchy. Hand it
to a policy object and there is nothing left for a category to be. The
CPU backend now has three collaborating pieces:

| | What it is | How many |
|---|---|---|
| **Runtime** (`cpu::Backend`) | frame binding, plan execution, timing, MPI, domain decomposition | one per run, shared by the whole driver tree |
| **Driver** (`cpu::CPUDriver<P>`) | one **loop node**: stages + plan + handlers | any number, nestable |
| **Handler** (`cpu::LoopHandler`) | asked once per iteration: stop or not? | any number, composable, user-definable |

What survives of the old table is the list of **presets** — the classes
users actually construct. A preset is *not* an architectural layer: it
is a class whose **constructor** attaches the default handler and the
infrastructure a given kind of run wants. It contributes defaults, not
semantics, and every one of them can be overridden.

| Preset | Phase enum | Default policy | Termination |
|---|---|---|---|
| `cpu::Evaluate` | `evaluate::Phase` | — | exactly one pass |
| `cpu::md::MD` | `md::Phase` | — | the step budget |
| `cpu::minimizer::Minimize` | `minimizer::Phase` | handler `ForceConverged` | convergence, or the budget |
| `cpu::vib::Hessian` | `vib::Phase` | — | fixed 6N displacements, then finalize |
| `cpu::mc::MC` | `sampling::Phase` | **outcome** `Metropolis` | the trial budget |
| `cpu::PathDriver` | `path::Phase` | handler `AllReplicasConverged` | every replica converged |

The wrong reflex is to subclass `Driver<>` for every algorithm
("LBFGSDriver", "CGDriver", "FIREDriver"). The right reflex is to
implement the algorithm as a Component and plug it into the existing
`Minimize` driver. The CPU codebase already paid this price once when
the original `cpu::minimizer::LBFGS` *driver class* was renamed to
`Minimize`, and `LBFGSStep` became `LBFGS`. The second-wrong reflex is
to reach for a new driver class because your run stops differently.
That is a `LoopHandler`, and it is about ten lines.

## The Hierarchy in One Diagram

```text
Driver<PhaseEnum, DeviceTag>            facade/driver.h — Plan/Backend ownership
    └─ cpu::CPUDriver<PhaseEnum>        cpu/driver/cpu_driver.h — THE loop node
         ├─ cpu::Evaluate               cpu/driver/evaluate.h
         ├─ cpu::md::MD                 cpu/md/md.h
         ├─ cpu::minimizer::Minimize    cpu/minimizer/minimize.h
         ├─ cpu::vib::Hessian           cpu/vib/hessian.h
         ├─ cpu::mc::MC                 cpu/mc/mc.h
         └─ cpu::PathDriver             cpu/driver/path_driver.h
```

Two levels, not four. `Driver<>` is the underlying primitive — it owns
the `Plan` and the `Backend`. `cpu::CPUDriver<P>` is the loop node, and
provides:

1. A `build()` override that runs `validate_components()` *before* the
   plan is compiled. Each preset overrides `validate_components()` to
   enforce its required role + capability set.
2. `setup_pass(frame)`, `step_once(iter_idx)`, `poll_handlers_(iter)`,
   `teardown_pass()` — the loop primitives each preset's `run(...)`
   composes, plus the continuation machinery (`begin_run_` /
   `flush_live_`) that MD and MC share.
3. **This loop's own** `IterationStatus` and `LoopMetrics` — owned by
   the node, borrowed by the Context (see below).
4. The two policy slots: `add_handler(...)` (plural, termination) and
   `set_outcome(...)` (singular, commit/rollback).
5. The categorical setter API (`add_potential`, `add_potential`, …, `set_thermo`).

The presets are thin: their constructors auto-add infrastructure
components (zeroing, kinetic accumulator, snapshot publisher) and
attach a default policy, and add nothing else. They each own a `run()`
that composes the primitives above — but the termination decision
inside it is delegated to whatever handlers are attached, which is why
adding a stopping rule never means adding a class.

## Termination Is a Handler

A `LoopHandler` ([`cpu/driver/loop_handler.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/driver/loop_handler.h))
is asked once per iteration and answers with an `ATV::Verdict`:
`Continue`, `Converged`, or `Stop`. Handlers are **plural** — several
can be attached, and their verdicts compose by severity, **Stop >
Converged > Continue**, so a run that has gone numerically insane
cannot report success by also happening to satisfy a force tolerance.

The handler sees a `cpu::LoopState` — iteration index, frame, Context,
metrics — and deliberately **never a phase**. That is what lets one
handler serve MD, Minimize, MC and Path alike. Three ship in-tree:
`ForceConverged` (the five-criteria geometry test), `AllReplicasConverged`,
and `StopOnNaN`. A user writes their own from outside the library:

```cpp
class StopWhenTooHot : public cpu::LoopHandler {
  public:
    explicit StopWhenTooHot(double t_max) : t_max_(t_max) {}
    const char* name() const override { return "stop_when_too_hot"; }
    ATV::Verdict after_iteration(const cpu::LoopState& s) override {
        return (temperature(s.frame) > t_max_) ? ATV::Verdict::Stop
                                                      : ATV::Verdict::Continue;
    }
  private:
    double t_max_;
};

sim.add_handler(StopWhenTooHot{2000.0});
```

This is LAMMPS's `fix halt`, and it was impossible to write before.

Under MPI the driver all-reduces the verdict across ranks, once, before
acting on it. Handler authors write `after_iteration` as if the process
were alone; ranks that disagreed about when to stop would deadlock at
the next collective, and that is not a hazard a handler author should
have to know about.

The loop's *other* policy slot is the `OutcomePolicy`
([`cpu/driver/outcome_policy.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/driver/outcome_policy.h)),
set with `set_outcome(...)` and **singular**. It answers a different
question — *does this iteration count?* — with `Commit` or `Rollback`.
Monte Carlo's accept/reject is the reason it exists: rejecting a trial
move does not end the run, it discards one iteration's work and keeps
going. `Metropolis` is one implementation; `MC` attaches it in its
constructor. Swapping in a GCMC or REMD acceptance rule is a
`set_outcome` call, not a new driver class.

## Kind (type) + Role (flags) — not Capability

There is **no `Capability` write-set** and no exclusive `ComponentRole`
enum. Two orthogonal mechanisms remain:

1. **Kind** — a C++ intermediate class (one header under
   `cpu/component/`: `Potential`, `HF`, `Semiempirical`, `PlaneWave`,
   `GPW`, `Bias`, `PathCoupling`, `Observer`, `Recorder`, `Integrator`,
   `Thermostat`, `Barostat`, `Optimizer`, `Constraint`, `Mover`).
   Presets validate with `dynamic_cast` to those types.
2. **Role** — a `RoleFlag` bitset on `role()` (`FrameImmutable`,
   `Deferred`, `ExternalIO`). Unset `Deferred` ⇒ main-loop required.
   See [`facade/component_role.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/facade/component_role.h).

Outputs are Frame keys (`kEnergy`, `kAtom_fx`, …), not capability bits.

- `cpu::scf::SCF` inherits `cpu::HF`. With `with_gradient(true)` it
  writes forces; Evaluate / MD / Minimize all accept any energy/force
  kind (`Potential`, `HF`, `Bias`, …).
- `cpu::md::FixLangevin` inherits `Thermostat` and writes additive
  forces; MD still special-cases Langevin by name when no Potential is
  present (pure Langevin diffusion).

When you write a new component: pick a **kind**, compose **Role** flags
only if you differ from the default `0`, and write the Frame fields you
own.

## The Component Lifecycle, Strictly

The base class
[`cpu::Component`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/component/component.h) defines
exactly three hooks:

```cpp
virtual void build  (const BuildContext&)       {}
virtual void launch (const LaunchContext&) = 0;  // pure
virtual void destroy(const DestroyContext&)     {}
```

The contract is:

1. **`build` runs once, before the first iteration.** All allocation
   goes here: per-thread reduction buffers, cell-list arrays,
   shell-pair lists, FFT plans, anything sized to `frame.capacity()`
   or derived from `frame.simbox`.
2. **`launch` runs once per iteration.** No allocation. Re-use only.
   This is enforced by code review and by regression tests like
   `tests/cpu/test_md_no_alloc.cpp`, which asserts that the inner MD
   loop performs zero `operator new` calls after `build`.
3. **`destroy` runs once, after the last iteration.** RAII handles most
   cleanup automatically (xtensor destructors, `std::vector` going out
   of scope), but if the component opens files, holds GPU handles in a
   PIMPL'd member, or registers with an external library, this is the
   place to close them.

The lifecycle is owned by the **Backend**, not by the component or the
driver. `cpu::Backend::build_components()` walks the compiled plan and
calls each component's `build()`. `cpu::Backend::step(idx)` walks the
plan once per phase per iteration and calls each component's
`launch()`. `cpu::Backend::destroy_components()` walks the plan in
reverse order and calls each component's `destroy()`. A component
neither knows nor cares which iteration is running — it sees only its
own `LaunchContext`.

Two practical consequences:

- **Per-step state lives on the component.** The L-BFGS history vectors
  $(s_i, y_i)$ live in `cpu::minimizer::LBFGS` member arrays. The Verlet
  pair list lives in `cpu::md::NeighborList::pairs_`. The
  predictor-corrector density-matrix history lives in `SCF`'s
  `aspc_alpha_` / `aspc_beta_` members. Drivers and the Frame stay
  stateless across iterations.
- **The Frame is the *only* shared mutable state.** Components
  communicate by reading and writing typed atom columns and typed metadata
  (`kEnergy`, `kKinetic`, `kStress`). There is no driver-side scratch
  buffer that a downstream component is supposed to read. If you find
  yourself wanting one, you have a layering problem; either lift it
  onto the Frame as a fundamental field, or own it inside the
  producing component and surface a `last_*()` accessor for observers
  to capture.

## The Iteration Status Channel

Two small structs carry the loop's verdict and the evidence for it.
`ATV::IterationStatus`
([`facade/iteration_status.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/facade/iteration_status.h))
is the verdict; `ATV::LoopMetrics`
([`facade/loop.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/facade/loop.h))
is the evidence.

```cpp
struct IterationStatus {
    bool converged       = false;
    bool stop_requested  = false;
    std::string reason;
    void reset() noexcept;
};

struct LoopMetrics {          // what the components measured this iteration
    double energy, max_force, rms_force;
    double max_displacement, rms_displacement, energy_change;
    bool   has_displacement;
    double replicas_converged;   // fraction in [0,1]; only a band writes it
};
```

**Status belongs to the loop, not to the Context.** Each loop node owns
its own `IterationStatus` and `LoopMetrics` as members
(`CPUDriver::status_`, `CPUDriver::metrics_`); the Context only
*borrows* the active loop's, through the non-owning `ctx.loop`
(`LoopServices{status, metrics}`). A nested driver re-points `ctx.loop`
at its own on entry and restores the parent's on exit, so **an inner
loop cannot stop its enclosing loop** — there is no pointer from here
to there. That is not a rule to remember; it is the reason a bug that
shipped twice (an optimizer writing `converged` unconditionally; the
first converged NEB image stopping the whole band) is now unwritable.

Ownership rules:

- The **loop node** owns the status. It calls `reset()` at the start of
  every iteration and polls its handlers at the end of every iteration.
- **Components measure, they do not judge.** `LBFGS` and `CG` publish
  `max_force`, `rms_force`, the displacements and ΔE into
  `ctx.loop.metrics`. They do not decide whether those numbers mean
  convergence — a `LoopHandler` does. `LBFGS::converged()` and
  `CG::converged()` are deleted; ask the driver: `opt.converged()` and
  `opt.stop_reason()`.
- A component may **request** a stop (`stop_requested = true` with a
  `reason`), but it **cannot declare convergence** — there is no field
  left for it to write. `converged` is set by the driver, and only when
  a handler returns `Verdict::Converged`.
- **SCF internals never touch any of this.** The SCF cycle's own
  convergence is internal to `cpu::scf::SCFSolver` — it has its own
  `DualConvergenceCheck` helper and its own loop.

There is no per-category table of "who reads which flag" any more,
because there are no categories: every loop node runs the same
`poll_handlers_` reduction, and what differs between an MD run and a
minimization is only which handlers happen to be attached.

## Why SCF Is an Inner Solver, Not a Driver

Neither backend implements SCF as a driver: CPU and CUDA both use the
`scf::SCFSolver` inner-solver shape behind an `SCF` component,
and the CPU plane-wave SCF follows the same shape
(`cpu::pw::PWSolver` + `PWEvaluator`). This was a substantive design
choice, and it is worth explaining because the asymmetry confuses new
developers.

There are two reasons CPU SCF is implemented as
`cpu::scf::SCFSolver` (a public low-level class) rather than as a
driver subclass:

1. **SCF is reused by AIMD and minimization.** `SCF::launch`
   constructs an `SCFSolver` once at `build()` time and calls
   `solver_->run(integrals)` on every MD step or every minimizer step.
   If SCF were itself a `Driver<>` subclass, AIMD would be one driver
   running another driver inside its `force` phase — an awkward
   nested-loop pattern that makes profiling, error reporting, and
   parallel-region coordination harder.
2. **The SCF cycle has no Plan structure visible to users.** No user
   ever wants to insert their own "after diagonalize" component into
   somebody else's SCF loop. The five inner steps — `HamiltonianBuilder`,
   `EigenSolver`, `DensityBuilder`, `DIISMixer`, `DualConvergenceCheck`
   — are an implementation, not a pluggable extension surface.

The SCF inner-step helpers therefore live in
[`include/atomiverse/cpu/scf/internal/`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/scf/internal/),
none of them inherit from `cpu::Component`, and none of them are
registered into any Plan. They are owned by `SCFSolver` as direct
members:

```cpp
class SCFSolver {
    SCFConfig config_;
    RealMatrix X_;                            // S^{-1/2}, built once

    internal::HamiltonianBuilder    hamiltonian_;
    internal::EigenSolver           eigen_;
    internal::DensityBuilder        density_;
    internal::DIISMixer             diis_;
    internal::DualConvergenceCheck  convergence_;

    SCFResult result_;
};
```

`SCFSolver::build(reference)` calls each helper's `prepare(...)` once.
`SCFSolver::run(integrals, guess)` runs the SCF cycle, calling each
helper's `step(...)` per iteration. The internals are public headers
(tests need to include them), but the **facade-public** production
headers — `evaluate.h`, `single_point.h`, `md.h`, `lbfgs.h` — must not
include from `cpu/scf/internal/`. That prevents the `internal::` types
from leaking into user-visible APIs.

If a user wants to run SCF directly without an outer driver, they
construct an `SCFSolver` themselves; if they want SCF as a force
provider in MD or minimization, they wrap it with `SCF`. Two
public surfaces, one inner solver, no nested drivers.

## Hot-Path Data Hygiene

`cpu::Frame`, `cpu::scf::IntegralPack`, and `cpu::scf::ShellList` carry
*only fundamental physical quantities*. Positions, velocities, forces,
masses, atom tags, image flags, integral tensors, shell descriptors —
those are fundamentals. Anything derivable — neighbor maps, atom-to-shell
indices, per-thread reduction buffers, predictor histories, cell-list
geometry — is *not* a fundamental and lives on the component or in the
runtime context.

This rule has bitten developers who tried to lift, e.g., a per-atom
bond-exclusion list onto the Frame. The right place is on the LJ
component (`excluded_by_tag_`, rebuilt from the per-atom topology columns at `build()`
time), not on the Frame. The Frame is the user-visible truth; derived
maps are implementation details that can change with the algorithm.

A correlated rule lives in `cpu::Context`. The active MD pair list and
the spatial-decomposition pointer are kept on the Context, *not* on the
Frame, for the same reason: SCF and minimizer drivers do not own a
neighbor list, so giving the Frame a neighbor-list pointer would be
wrong everywhere outside MD.

## The Box Abstraction

Every CPU file that loops over pairs of atoms, wraps a position into
the cell, or applies the minimum-image convention reads
[`cpu::Box`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/frame/box.h). The full set of
methods is:

```cpp
double length(int d)        const;   // edge length along axis d
double volume()             const;
void   min_image(double& dx, double& dy, double& dz) const;
double wrap_axis(double q, int d) const;
void   scale_diagonal(double sx, double sy, double sz);
void   scale_isotropic(double mu);
void   scale_axis(int d, double mu);
void   set_cell(const Mat3& M);
```

The Box also caches `inv_cell`, per-axis `inv_L`, `half_L`, and an
`ortho` flag. The cache is rebuilt only on `set_cell()` — never inside
`launch()` — so a force kernel reads the cached inverses without paying
for a recomputation each time. Nine call sites that used to inline
their own `static inline minimum_image()` helper have been folded into
this single class.

There is one deliberate exception: the Lennard-Jones xsimd ket-loop
(`src/cpu/potentials/pair_lj_cut.cpp`, vectorized inner loop) inlines the
orthorhombic min-image branch instead of calling through `Box::min_image`. This
is because member loads through `this` defeat the compiler's ability to
hoist the per-axis $L^{-1}$ out of the loop, and the cost was measured
at ~7% per MD step. The xsimd batch overload `min_image_batch<B>` in
[`cpu/frame/box_simd.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/frame/box_simd.h) gives non-LJ
SIMD callers a vectorized path without that landmine.

The takeaway when adding new pair-loop code: prefer `cpu::Box`
helpers; only inline the math when you can show the compiler is failing
to hoist the loads, and document the reason.

## Validation, Phase Compilation, and Plan Sorting

The journey from `driver.add(Component{...})` to "components running in
order" goes through these stages:

1. **Registration.** `driver.add(Component{...})` owns the component, reads its
   `default_bindings()`, and appends `(component, phase, order)` triples to the
   driver's `pending_` vector. `driver.on(p).add(Component{...}, ord)` owns and
   skips the defaults. The corresponding `add_ref` calls borrow an existing
   instance; plain `add(lvalue)` is deleted.
2. **Validation.** `cpu::CPUDriver::build()` calls
   `validate_components()` (overridden by each preset) which walks
   `pending_` and checks role + capability requirements.
3. **Plan compilation.** The base `Driver<>::build()` then asks the
   `PlanBuilder` to sort `pending_` by `(phase, order, insertion_index)`
   and emit a compact `Plan` of `StageSlot`s — one slot per phase, each
   holding a vector of component pointers in execution order.
4. **Backend lifecycle.** `Backend::init()` populates Context fields
   (thread count, MPI size, executor). `Backend::load(plan)` adopts
   ownership of the compiled plan. `Backend::set_frame(frame)` records
   the user's frame. `Backend::build_components()` walks the plan and
   calls each component's `build()`. `Backend::step(idx)` walks the
   plan once per iteration. `Backend::destroy_components()` reverses.

A few sharp edges to know:

- **Phases are strongly typed.** `md::Phase`, `minimizer::Phase`,
  `evaluate::Phase`, and `vib::Phase` are independent enums. A
  component bound by `default_bindings()` to `md::Phase::pre_force`
  cannot accidentally end up in the minimizer's `force` phase — the
  enum types are different, and the driver template is parameterized
  over the phase enum.
- **Order is a signed integer.** Negative orders run before positive
  orders within the same phase. Infrastructure components (zeroers,
  kinetic accumulator, halo exchange) use orders in $[-1000, -10]$;
  user components default to $0$; analysis publishers and recorders
  sit at $+10$ to $+25$. When in doubt, leave order at the default and
  let the constructor of `cpu::md::MD` / `cpu::minimizer::Minimize`
  place the infrastructure for you.
- **Validation runs *before* plan compilation.** `validate_components`
  iterates over `pending_`, not over a sorted plan, because at
  validation time the plan does not exist yet. Components with
  multiple `default_bindings()` entries appear once per binding in
  `pending_` — duplicate-tolerant predicates are expected.

## What Counts as a New Algorithm vs. a New Driver

A short rubric, with examples:

- **New integrator** (e.g. Beeman, Leapfrog with constraints): new
  Component, role `Integrator`, plug into existing `cpu::md::MD`. No
  new driver class.
- **New thermostat** (e.g. CSVR / Bussi): new Component, role
  `Thermostat`, plug into existing `cpu::md::MD`. No new driver class.
- **New optimizer** (e.g. FIRE): new Component, role `Optimizer`, plug
  into existing `cpu::minimizer::Minimize`. No new driver class.
- **New SCF inner step** (e.g. ADIIS replacing DIIS): new
  `cpu::scf::internal::*` helper, owned as a member of
  `cpu::scf::SCFSolver`. No new component, no new driver.
- **New evaluator** (e.g. EAM potential, MLIP): new Component, role
  `Evaluator`, capability `Energy | Force | Stress` (whichever it
  produces). Plugs into `cpu::md::MD` *and* `cpu::minimizer::Minimize`
  *and* `cpu::Evaluate`.
- **New stopping rule** (e.g. "halt above 2000 K", "give up after an
  hour", "stop when the band is flat"): new `cpu::LoopHandler`,
  attached with `add_handler(...)`. No new driver class, and no edit to
  any algorithm.
- **New acceptance rule** (e.g. GCMC insertion, a REMD exchange): new
  `cpu::OutcomePolicy`, attached with `set_outcome(...)`. No new driver
  class.
- **New loop shape** (a genuinely different stage sequence — a new
  phase enum): subclass `cpu::CPUDriver<YourPhase>` and specialize
  `cpu::PhaseTraits` for it. This is rare, and it is the *only* case
  that justifies a new driver class.

When the right answer is "new component", do not invent a new driver
class. The validation rules in the existing presets already enforce the
structural pieces — your new component's role and capability tags will
plug into them automatically.

## Things to Read Next

- `tests/cpu/test_minimize_driver_swap_optimizer.cpp` — the
  canonical hot-swap test for two `Optimizer`-role components.
- `tests/cpu/test_md_no_alloc.cpp` — the regression test for "no
  allocation in `launch()`".
- `tests/cpu/scf/test_aimd_conservation.cpp` — energy / momentum /
  bond-oscillation conservation under composed `MD + SCF`.
- [Extending CPU Subsystems](extending-cpu-subsystems.md) — the
  practical recipe for adding new components in each of MD, SCF, and
  the minimizer.
