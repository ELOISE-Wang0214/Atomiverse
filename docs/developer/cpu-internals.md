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

## Why Driver Categories Exist

User-facing simulation in Atomiverse comes in five flavours, and the
clean way to think about them is: **the algorithm does not pick the
category — the termination rule does.**

| Category | Termination | Concrete CPU class |
|---|---|---|
| `cpu::EvaluateDriver` | exactly one pass | `cpu::Evaluate` |
| `cpu::DynamicsDriver` | full step budget; honours `outer_status.stop_requested` | `cpu::md::MD` |
| `cpu::OptimizeDriver` | early-stop when `outer_status.converged` | `cpu::minimizer::Minimize` |
| `cpu::ResponseDriver` | fixed 6N-displacement budget, then finalize; no early stop | `cpu::vib::Hessian` |
| `cpu::SamplingDriver` | trial-move budget with inline Metropolis accept/reject | `cpu::mc::MC` |

If you ever want to add a new category, the test is whether the new
loop has a *genuinely distinct* termination semantic. NEB qualifies —
the loop terminates when *every image* of a coupled chain has
converged, a coupled-convergence rule none of the existing categories
model. Sampling and Response each earned their row exactly this way
(mid-run rejection logic; a fixed displacement schedule). A new
"molecular mechanics" driver, however, is *not* a new category — it is
a Dynamics driver with a different evaluator at the force phase.

The wrong reflex is to subclass `Driver<>` for every algorithm
("LBFGSDriver", "CGDriver", "FIREDriver"). The right reflex is to
implement the algorithm as a Component and plug it into the existing
`Minimize` driver. The CPU codebase already paid this price once when
the original `cpu::minimizer::LBFGS` *driver class* was renamed to
`Minimize`, and `LBFGSStep` became `LBFGS`.

## The Hierarchy in One Diagram

```text
Driver<PhaseEnum, DeviceTag>            facade/driver.h — Plan/Backend ownership
    └─ cpu::CPUDriver<PhaseEnum>        cpu/driver/cpu_driver.h — CPU backend layer
         ├─ cpu::EvaluateDriver         cpu/driver/evaluate_driver.h
         │     └─ cpu::Evaluate         cpu/driver/evaluate.h
         ├─ cpu::DynamicsDriver         cpu/driver/dynamics_driver.h
         │     └─ cpu::md::MD           cpu/md/md.h
         └─ cpu::OptimizeDriver         cpu/driver/optimize_driver.h
               └─ cpu::minimizer::Minimize  cpu/minimizer/minimize.h
```

`Driver<>` is the underlying primitive — it owns the `Plan` and the
`Backend`. `cpu::CPUDriver<P>` adds CPU-specific glue:

1. A `build()` override that runs `validate_components()` *before* the
   plan is compiled. Subclasses override `validate_components()` to
   enforce their required role + capability set.
2. `setup_pass(frame)`, `step_once(iter_idx)`, `teardown_pass()` —
   protected helpers each category's `run(...)` calls in the right
   pattern.
3. Read/write access to `cpu::Context::outer_status` so the category
   loop can poll `converged` / `stop_requested`.

User-facing CPU drivers (`Evaluate`, `MD`, `Minimize`) are thin: their
constructors auto-add infrastructure components (zeroing, kinetic
accumulator, snapshot publisher) and add nothing else. They do *not*
override `run()`; the category base owns the termination semantics.
This is the load-bearing rule. If you find yourself wanting to override
`run()` in a user-facing class, you almost certainly want a new
*category* base, not a new concrete driver — and you should think hard
about whether the existing three already cover your case.

## Two Tagging Axes, Not One

The single biggest design decision in the component model is splitting
*role* and *capability* into independent axes. The definitions are in
[`facade/component_role.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/facade/component_role.h):

```cpp
enum class ComponentRole : uint8_t {
    None, Evaluator, Integrator, Optimizer,
    Thermostat, Barostat, ConvergenceCheck, Recorder, Observer
};

enum class Capability : uint32_t {
    None          = 0,
    Energy        = 1u << 0,
    Force         = 1u << 1,
    Stress        = 1u << 2,
    Charge        = 1u << 3,
    Dipole        = 1u << 4,
    DensityMatrix = 1u << 5,
    Orbitals      = 1u << 6,
};
```

Why two axes? Because the same component class is, in different
contexts, the right answer to different validation questions.

- `cpu::scf::SCFEvaluator` is always `role() == Evaluator`. Its
  `capabilities()` is `Energy | DensityMatrix | Orbitals` by default,
  and *also* includes `Force` when `with_gradient(true)` is set.
  - In `cpu::Evaluate` (single-point energy), the `Force` bit can be
    off and the driver still validates: it just needs an Evaluator with
    *any* non-empty capability set.
  - In `cpu::md::MD`, the same class with `with_gradient(true)` passes
    Dynamics's "something must produce `Capability::Force`" check.
  - In `cpu::minimizer::Minimize`, the same class passes Optimize's
    "Energy or Force" check.

If role and capability were collapsed into a single tag, this class
would need to be three classes, or each driver would need to special-case
SCF, or both. The split is the principled way out.

A second example: `cpu::md::FixLangevin` declares
`role() == Thermostat` (it is a temperature controller), but its
`capabilities()` is `Capability::Force` because its `launch()` writes
additive friction + random forces into `fx/fy/fz`. That one bit is what
lets a *pure* Langevin run — diffusion in a featureless box, no LJ at
all — satisfy `DynamicsDriver`'s "something writes forces" build-time
check. Roles describe structure; capabilities describe data flow.

When you write a new component, ask both questions independently: what
slot does this fill? what fields does it write?

## The Component Lifecycle, Strictly

The base class
[`cpu::ComponentBase`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/component/component.h) defines
exactly three hooks:

```cpp
virtual void build  (const BuildContext&)       {}
virtual void launch (const LaunchContext&) = 0;  // pure
virtual void destroy(const DestroyContext&)     {}
```

The contract is:

1. **`build` runs once, before the first iteration.** All allocation
   goes here: per-thread reduction buffers, cell-list arrays,
   shell-pair lists, FFT plans, anything sized to `frame.atoms.size()`
   or `frame.box`.
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
  predictor-corrector density-matrix history lives in `SCFEvaluator`'s
  `aspc_alpha_` / `aspc_beta_` members. Drivers and the Frame stay
  stateless across iterations.
- **The Frame is the *only* shared mutable state.** Components
  communicate by reading and writing `frame.atoms.*`, `frame.scalars.energy`,
  `frame.scalars.kinetic`, and `frame.scalars.stress`. There is no driver-side scratch
  buffer that a downstream component is supposed to read. If you find
  yourself wanting one, you have a layering problem; either lift it
  onto the Frame as a fundamental field, or own it inside the
  producing component and surface a `last_*()` accessor for observers
  to capture.

## The Iteration Status Channel

`cpu::Context` carries one piece of structured cross-component state:
`outer_status` of type `ATV::IterationStatus`
([`facade/iteration_status.h`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/facade/iteration_status.h)).

```cpp
struct IterationStatus {
    bool converged       = false;
    bool stop_requested  = false;
    std::string reason;
    void reset() noexcept;
};
```

Ownership rules:

- The **outer category driver** owns the status. It calls `reset()` at
  the start of every iteration and polls it at the end of every
  iteration.
- **`ConvergenceCheck`-tagged components** (in practice, the optimizer
  components — `LBFGS` and `CG`) *write* to it. Their convention is to
  set `outer_status.converged = true` when their five-criterion check
  passes.
- **Anyone** can write `outer_status.stop_requested = true` and a
  `reason` string. NaN detectors, wall-clock guards, and future
  `StopCondition` components are the typical clients.
- **SCF internals never touch `outer_status`.** The SCF cycle's own
  convergence is internal to `cpu::scf::SCFSolver` — it has its own
  `DualConvergenceCheck` helper and its own loop. The outer status
  channel is for the *outer* driver's loop.

The two flags are interpreted differently by category:

| Category | Reads `converged` | Reads `stop_requested` |
|---|---|---|
| Evaluate | (one pass; no loop) | (irrelevant) |
| Dynamics | **ignored** (time evolution doesn't "succeed") | breaks the loop |
| Optimize | breaks the loop | breaks the loop |

That is the smallest possible specification of "what does an outer loop
do?" If you are adding a new category, write down the same row.

## Why SCF Is an Inner Solver, Not a Driver

Neither backend implements SCF as a driver: CPU and CUDA both use the
`scf::SCFSolver` inner-solver shape behind an `SCFEvaluator` component,
and the CPU plane-wave SCF follows the same shape
(`cpu::pw::PWSolver` + `PWEvaluator`). This was a substantive design
choice, and it is worth explaining because the asymmetry confuses new
developers.

There are two reasons CPU SCF is implemented as
`cpu::scf::SCFSolver` (a public low-level class) rather than as a
driver subclass:

1. **SCF is reused by AIMD and minimization.** `SCFEvaluator::launch`
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
none of them inherit from `cpu::ComponentBase`, and none of them are
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
provider in MD or minimization, they wrap it with `SCFEvaluator`. Two
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
component (`excluded_by_tag_`, rebuilt from `frame.atoms (per-atom bonded topology)` at `build()`
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

The journey from `driver.add(component)` to "components running in
order" goes through these stages:

1. **Registration.** `driver.add(c)` reads `c.default_bindings()` and
   appends `(component, phase, order)` triples to the driver's
   `pending_` vector. `driver.on(p).add(c, ord)` skips the defaults and
   adds one explicit triple.
2. **Validation.** `cpu::CPUDriver::build()` calls
   `validate_components()` (overridden by each category base) which
   walks `pending_` and checks role + capability requirements.
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
- **New runtime semantic** (e.g. NEB, REMD, MC): potentially a new
  category base. Must be justified by a *termination rule* that does
  not fit Evaluate / Dynamics / Optimize.

When the right answer is "new component", do not invent a new driver
class. The validation rules in the existing categories already enforce
the structural pieces — your new component's role and capability tags
will plug into them automatically.

## Things to Read Next

- `tests/cpu/test_minimize_driver_swap_optimizer.cpp` — the
  canonical hot-swap test for two `Optimizer`-role components.
- `tests/cpu/test_md_no_alloc.cpp` — the regression test for "no
  allocation in `launch()`".
- `tests/cpu/scf/test_aimd_conservation.cpp` — energy / momentum /
  bond-oscillation conservation under composed `MD + SCFEvaluator`.
- [Extending CPU Subsystems](extending-cpu-subsystems.md) — the
  practical recipe for adding new components in each of MD, SCF, and
  the minimizer.
