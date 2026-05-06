---
title: CPU Driver taxonomy — Evaluate / Dynamics / Optimize (+ SCF inner solver)
status: done
created: 2026-05-06
revised: 2026-05-06
landed: 2026-05-06
supersedes: cpu-scf-driver-alignment, phase2-singlepoint-driver
followups: cuda-driver-taxonomy (planned), phase 5b (IntegralPack/ShellList split + 5 monolith file splits — deferred)
---

# CPU Driver Taxonomy

## Summary

The current CPU driver layer has the wrong axis of variation. Concrete
classes (`cpu::md::MD`, `cpu::minimizer::LBFGS`, `cuda::minimizer::CG`,
`cuda::minimizer::SD`, `cuda::scf::SCF`) name themselves after the
**algorithm** they happen to embed (LBFGS, CG, SD), conflate **outer
task semantics** with **inner solver mechanics**, and inherit a single
template `Driver<PhaseEnum, DeviceTag>` whose `run(frame, n_steps)`
always runs *exactly* `n_steps` iterations regardless of physical
meaning. CPU SCF has no driver at all — it lives as the free function
`scf_solve()`.

This spec separates two distinct conceptual axes the current code
mashes together.

### Axis 1 — outer task semantic (user-facing driver)

What is the user *asking the simulation to do*? For the current CPU
simulation facade, three answers cover the supported workflows:

1. **Evaluate** — run evaluator components once on the current state
   (single-point energy / force / property / electronic structure).
2. **Dynamics** — advance the state for a requested step budget.
3. **Optimize** — iteratively update the state until convergence or
   `max_iter`.

These are the **three current** user-facing driver categories. They
are distinguished by termination semantics, which is a hard contract
of the type. **Future categories** (Sampling for MC/HMC/REMD, Path
for NEB/string method, Response for finite-difference Hessian /
phonon) may be introduced **only** when they have a genuinely
distinct outer-loop semantic — never to wrap an *algorithm* under a
new class name. This rule is enforced at review time by the
binding-rule list below.

### Axis 2 — inner solver (component-internal machinery)

How does a single force / energy / property *get computed*? Some
evaluators (LJ, harmonic bond) compute directly. Others
(`SCFEvaluator`) need an **inner self-consistency loop** to converge
internal variables (density matrix) before they can return a value.
That inner loop is **not** an outer driver category. In ordinary user
workflows it appears through `Evaluate + SCFEvaluator`. For
debugging, low-level testing, and advanced scripting, the inner
solver is **public**:

```cpp
ATV::cpu::scf::SCFSolver solver{config};
solver.build(integrals);
const auto& result = solver.run(integrals);   // density, MO coeffs, eigenvalues
```

`SCFSolver` is a regular class — **not a `Driver<...>` subclass**, no
outer phase enum, no registration into outer Plans. Its internal
step helpers (Hamiltonian build / diagonalize / density build / DIIS
mix / convergence check) live in `namespace cpu::scf::internal::` and
are not included from any facade-public header — but tests may
include them directly to unit-test individual helpers.

### Three concrete pathologies fixed

1. **Wrong subclass axis.** LBFGS / CG / SD / FIRE today are *driver
   classes*. They are *Optimizer policies* that should plug into a
   single `OptimizeDriver`. Same for NoseHoover / Langevin
   (Thermostat policies under `DynamicsDriver`). Subclassing per
   algorithm hard-codes the algorithm into the type system.
2. **Wrong loop semantics.** Evaluate / Dynamics / Optimize have
   genuinely different termination rules. Today they all share
   `Driver::run()` which silently ignores convergence and forces
   wasteful full-budget runs.
3. **Confused SCF placement.** A previous draft introduced
   `cpu::scf::SinglePoint : SCFDriver`. But "single-point" is the
   *outer* task ("compute things at this geometry"), and "SCF" is the
   *inner* solver. They sit on different axes. Single-point is one
   `Evaluate` call; SCF is what happens inside one step of that
   `Evaluate` when the evaluator is a wave function.

### Concrete drivers under the new taxonomy

| New name | Base | Replaces |
|---|---|---|
| `cpu::Evaluate` | `EvaluateDriver` | (new — was `MD.run(frame, 1)` hack or direct `scf_solve()` call) |
| `cpu::md::MD` | `DynamicsDriver` | unchanged identity, new parent |
| `cpu::minimizer::Minimize` | `OptimizeDriver` | `cpu::minimizer::LBFGS` (renamed) |
| **(no class)** | — | `cpu::scf::SinglePoint` — never exists |

Optional **free function** convenience: `cpu::single_point(frame, evaluator)`
may be added as a thin wrapper that constructs a transient `Evaluate`
and runs it. **Not** a class, **not** a driver type.

Single-point energy of a wave function:

```cpp
cpu::Evaluate eval;
auto& scf = eval.add(cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals)});
eval.build();
eval.run(frame);

double energy_eV = frame.energy;
const auto& orbitals = scf.last_result();   // density, MO coeffs, eigenvalues
```

Single-point energy of a force field:

```cpp
cpu::Evaluate eval;
eval.add(cpu::md::LJ{eps, sigma, cutoff});
eval.add(cpu::md::NeighborList<>{cutoff});
eval.run(frame);
```

Same pattern, different evaluator. The user learns one shape.

### Algorithms become components with role tags

| Component | Role | Plugs into |
|---|---|---|
| `cpu::minimizer::LBFGS` (was `LBFGSStep`) | Optimizer | `Minimize` |
| `cpu::minimizer::CG` *(new)* | Optimizer | `Minimize` |
| `cpu::minimizer::FIRE` *(future)* | Optimizer | `Minimize` |
| `cpu::md::VelocityVerlet` | Integrator | `MD` |
| `cpu::md::NoseHoover` | Thermostat | `MD` |
| `cpu::md::FixLangevin` | Thermostat | `MD` |
| `cpu::md::Barostat` | Barostat | `MD` |
| `cpu::md::LJ` | Evaluator (cap: Energy + Force) | `MD`, `Minimize`, `Evaluate` |
| `cpu::scf::SCFEvaluator` | Evaluator (cap: Energy + DensityMatrix + Orbitals + optional Force) | `MD`, `Minimize`, `Evaluate` |

### Two orthogonal tags: role and capability

Roles are **structural** ("what slot in the outer loop does this fill").
Capabilities are **what data the component produces** ("does it write
energy? force? density matrix?"). Splitting them is what lets
`SCFEvaluator` be the same class whether the user wants energy-only,
energy+density, or energy+force.

```cpp
enum class ComponentRole : uint8_t {
    None,
    Evaluator,         // computes things — capabilities() says what
    Integrator,
    Optimizer,
    Thermostat,
    Barostat,
    ConvergenceCheck,
    Recorder,
    Observer,
};

// Bitset of computable quantities. Add new bits as the project grows.
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
// + bitwise operators

class ComponentBase {
public:
    virtual ComponentRole role() const { return ComponentRole::None; }
    virtual Capability    capabilities() const { return Capability::None; }
    // ...
};
```

Validation rules at category-base `build()`:

| Driver | Required (any of) | Required (exactly one of) |
|---|---|---|
| `EvaluateDriver` | ≥ 1 component with `role()==Evaluator` AND `capabilities() != None` | — |
| `DynamicsDriver` | ≥ 1 Evaluator with `Force` capability | exactly 1 Integrator |
| `OptimizeDriver` | ≥ 1 Evaluator with `Energy ∪ Force` capability | exactly 1 Optimizer |

`SCFEvaluator` advertises `Energy | DensityMatrix | Orbitals` by
default and adds `Force` when the user enables gradient
(`SCFEvaluator{config}.with_gradient()` or equivalent — concrete API
decided in Phase 4). This is what makes it valid in both `Evaluate`
(energy-only single-point) and `MD` (AIMD requires `Force`).

### Termination & stop conditions

| Driver | Normal end | Permitted aborts |
|---|---|---|
| `EvaluateDriver` | exactly one pass | error during launch |
| `DynamicsDriver` | exhausted step budget | error during launch; explicit `StopCondition` (future) |
| `OptimizeDriver` | `IterationStatus::converged` set OR `max_iter` reached | error; explicit stop |

Dynamics **does not interpret convergence as success** — that's the
hard contract distinguishing it from Optimize. It may still abort on
genuine errors (NaN detection, integrator blowup), wall-time limits,
user interrupts, or future explicit `StopCondition` mechanisms — but
those mechanisms are out of scope for this spec.

### Iteration status replaces a global atomic flag

Instead of `std::atomic<bool> Context::converged`, the outer status
is a small structured object:

```cpp
struct IterationStatus {
    bool converged = false;          // Optimize success
    bool stop_requested = false;     // explicit early termination
    std::string reason;              // human-readable detail
    void reset() noexcept { converged = false; stop_requested = false; reason.clear(); }
};
```

Only the **outer** category driver owns this object. `OptimizeDriver`
resets it before each iteration, runs the iteration, then breaks if
either bit is set. `DynamicsDriver` resets and breaks only on
`stop_requested` (never on `converged`). `EvaluateDriver` does not
loop, so it ignores everything except errors.

Components that decide convergence (e.g. SCF-internal
`DualConvergenceCheck`, optimizer convergence) **do not** touch the
outer status. SCF internals manage their own internal convergence
inside `SCFSolver::run()`. Outer optimizer ConvergenceCheck components
read frame state and call `context.outer_status.converged = true;`.

This is opt-in (only ConvergenceCheck-tagged components are wired to
the outer status) and avoids the multi-writer / atomic-coupling
hazards of a global flag.

### Phase enums

- New `ATV::evaluate::Phase`: `prepare, evaluate, finalize, COUNT`.
  No `end_of_step` — `Evaluate` has no step. Live in `facade/stage.h`.
- `ATV::md::Phase`: unchanged.
- `ATV::minimizer::Phase`: unchanged.
- `ATV::scf::Phase`: kept temporarily for `cuda::scf::SCF`; marked
  "legacy CUDA-only; CPU SCF is an inner solver" in the header.
  Removed when the CUDA migration follow-up lands.

## Binding rules (apply across the whole spec)

- **Three current outer driver categories: `EvaluateDriver`,
  `DynamicsDriver`, `OptimizeDriver`.** Adding a fourth requires a new
  spec demonstrating a genuinely distinct outer-loop semantic
  (e.g. Sampling = stochastic ensemble walks, Path = chain-of-states,
  Response = perturbation-response sweeps). It must not be added to
  wrap an algorithm.
- **Do not introduce `SinglePointDriver` or `cpu::scf::SinglePoint`.**
  Single-point energy is `Evaluate + SCFEvaluator`. A free function
  `cpu::single_point(frame, evaluator)` is acceptable as syntactic
  sugar; a class is not.
- **SCF is not an outer driver category.** It is represented by
  `cpu::scf::SCFSolver` and used by SCF-capable evaluator components.
  `SCFSolver` is a public low-level API for testing, debugging, and
  advanced workflows. It does not inherit `Driver<...>`, owns no
  outer phase enum, and never registers internal SCF helpers as
  outer-Plan components.
- **Algorithm choices are component classes, not driver subclasses:**
  LBFGS / CG / FIRE are Optimizer components; VelocityVerlet is an
  Integrator; Langevin / NoseHoover are Thermostat components.
- **Evaluators are validated by capability, not by hard-coded role.**
  Energy-only, property-only, and force-providing evaluators must
  all be valid in `Evaluate`. `DynamicsDriver` and `OptimizeDriver`
  require *at least one* evaluator with `Force` (or `Energy`)
  capability — they don't care which class provides it.
- **`OptimizeDriver` is the only current category that treats
  convergence as normal successful termination.** Dynamics may abort
  on errors or explicit `StopCondition`s but does not use convergence
  as its normal success path.
- **Phase enums are owned by the category base.** `EvaluateDriver`
  owns `evaluate::Phase`; `DynamicsDriver` keeps `md::Phase`;
  `OptimizeDriver` keeps `minimizer::Phase`. CUDA-only `scf::Phase`
  survives one release window.
- **`Driver<PhaseEnum, DeviceTag>::run` is removed from the public
  surface.** The three category bases own `run(...)` in their own
  termination-correct shapes. The template stays as a private plan /
  dispatch primitive.
- **Composition still flows through `add(...)` and `on(phase).add(...)`.**
- **No new physics.** PySCF / libcint / LAMMPS remain ground truth.

## Domain basis

No new domain content. All algorithms (Velocity Verlet, Nose-Hoover,
MTK NPT, L-BFGS, DIIS, RHF/UKS) keep their current implementations
and references. The change is purely structural.

## Architectural design

### Category base contracts

Sketches (final API at implementation time):

```cpp
namespace ATV::cpu {

template <typename PhaseEnum>
class CPUDriver : public Driver<PhaseEnum, cpu::device_tag> {
protected:
    void run_one_iteration(cpu::Frame& frame, int iter_index);
    IterationStatus& outer_status();         // resets per-iter
};

class EvaluateDriver : public CPUDriver<evaluate::Phase> {
public:
    void run(cpu::Frame& frame);             // exactly one pass
};

class DynamicsDriver : public CPUDriver<md::Phase> {
public:
    void run(cpu::Frame& frame, int n_steps); // full budget unless stop_requested
};

class OptimizeDriver : public CPUDriver<minimizer::Phase> {
public:
    int run(cpu::Frame& frame, int max_iter); // returns iters actually run
};

}  // namespace ATV::cpu
```

The three `run` overloads share `run_one_iteration` as a tiny core
but the type-system difference is the user-visible contract.

### SCFSolver design

```cpp
namespace ATV::cpu::scf {

/// Public low-level inner solver. Not a Driver subclass.
///
/// Used internally by SCFEvaluator. Also usable directly for unit
/// tests, density-matrix dumps, DIIS comparisons, or initial-guess
/// experiments — but it is not registered into any outer Plan and
/// owns no outer phase enum.
class SCFSolver {
public:
    explicit SCFSolver(SCFConfig config);

    /// Pre-allocate persistent workspaces sized for `reference_shape`.
    /// Must be called before run().
    void build(const IntegralPack& reference_shape);

    /// Execute the self-consistency loop. Performs zero allocations
    /// in steady state (after build()). Returns reference to the
    /// owned SCFResult, valid until the next run() call.
    SCFResult& run(const IntegralPack& integrals,
                   const SCFInitialGuess* guess = nullptr);

    const SCFResult& last_result() const;

private:
    SCFConfig config_;
    internal::HamiltonianBuilder hamiltonian_;
    internal::EigenSolver eigen_;
    internal::DensityBuilder density_;
    internal::DIISMixer diis_;
    internal::DualConvergenceCheck convergence_;
    SCFResult result_;
    // pre-allocated workspaces
};

}  // namespace
```

Internal helpers in `namespace cpu::scf::internal::`:

| Helper | Step |
|---|---|
| `internal::HamiltonianBuilder` | F = Hcore + J − αK + V_xc |
| `internal::EigenSolver` | dsygvd on (F, S) |
| `internal::DensityBuilder` | P = 2 C_occ C_occᵀ (or αβ) |
| `internal::DIISMixer` | Pulay extrapolation |
| `internal::DualConvergenceCheck` | \|ΔE\| < tol ∧ ‖err‖_F < √tol |

These helpers each have a `build()/step()` shape so they pre-allocate
F / J / K / V_xc / DIIS scratches once and rewrite in place. They do
**not** inherit `cpu::ComponentBase` (no `default_bindings`, no
`role()`). They are not included from any header outside
`include/atomiverse/cpu/scf/internal/`. Test code is allowed and
encouraged to include them directly for unit-level testing.

`scf_solve()` is reduced to a ≤ 30-line shim that constructs a
transient `SCFSolver` and runs it, OR is deleted (decision recorded
in NOTES.md at Phase 4).

### Single-point invocation pattern

```cpp
auto integrals = build_integrals(frame, BasisSet::Def2SVP);
cpu::Evaluate eval;
auto& scf = eval.add(cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals)}
                      .with_gradient(false));   // energy/density only
eval.build();
eval.run(frame);

double energy_eV = frame.energy;
const auto& result = scf.last_result();   // density, MOs, eigenvalues
```

For the (rare) case where the caller wants the SCF cycle directly:

```cpp
cpu::scf::SCFSolver solver{config};
solver.build(integrals);
const auto& result = solver.run(integrals);  // public low-level API
```

Mixed-level workflow (PBE/def2-SVP optimization → B3LYP/def2-TZVP
single-point) is two driver constructions:

```cpp
cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{...});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{cfg_A}.with_gradient(true));
opt.build();
opt.run(frame, /*max_iter=*/200);

cpu::Evaluate eval;
eval.add(cpu::scf::SCFEvaluator{cfg_B}.with_gradient(false));
eval.build();
eval.run(frame);
```

No glue class, no `SinglePoint`, no `OptimizeThenSinglePoint`.

### Backwards compatibility (one release window)

- `cpu::minimizer::LBFGS` (the *old* driver name) → `[[deprecated]]`
  alias for `cpu::minimizer::Minimize`. Removed in Phase 6.
- `cpu::minimizer::LBFGSStep` (the old component name) → `[[deprecated]]`
  alias for `cpu::minimizer::LBFGS` (the new component name). Removed
  in Phase 6.
- `cpu::scf::scf_solve()` → ≤ 30-line shim that constructs an
  `SCFSolver` and runs it. Removed in Phase 6 unless an external
  caller surfaces.
- `cpu::md::MD` keeps its name; only its parent changes.

### Out of scope (deliberate)

- **CUDA backend.** `cuda::md::MD`, `cuda::minimizer::CG/SD`,
  `cuda::scf::SCF` keep their current shape. Follow-up
  `cuda-driver-taxonomy.md` migrates them after CPU lands.
- **Future driver categories (Sampling, Path, Response).** Out of
  scope; the binding rules permit them but this spec doesn't add
  any.
- **`StopCondition` mechanism for Dynamics aborts.** Out of scope;
  Phase 1 ships only the `IterationStatus::stop_requested` field;
  no explicit components write it yet.
- **Thermostat / Barostat policy unification.** Stays component-level.
- **Python bindings update.** Tracked in `python-bindings.md`.
- **Performance optimization.** ±5% wall-time gate only.

## Files

### Created

```
include/atomiverse/cpu/cpu_driver.h        (CPUDriver)
include/atomiverse/facade/component_role.h         (ComponentRole + Capability + role()/capabilities())
include/atomiverse/facade/iteration_status.h       (IterationStatus struct)
include/atomiverse/cpu/evaluate_driver.h
include/atomiverse/cpu/dynamics_driver.h
include/atomiverse/cpu/optimize_driver.h
include/atomiverse/cpu/evaluate.h                  (cpu::Evaluate concrete driver)
include/atomiverse/cpu/single_point.h              (free function helper, optional)
include/atomiverse/cpu/scf/scf_solver.h            (public SCFSolver)
include/atomiverse/cpu/scf/internal/hamiltonian_builder.h
include/atomiverse/cpu/scf/internal/eigen_solver.h
include/atomiverse/cpu/scf/internal/density_builder.h
include/atomiverse/cpu/scf/internal/diis_mixer.h
include/atomiverse/cpu/scf/internal/dual_convergence_check.h
include/atomiverse/cpu/scf/shell_list.h
include/atomiverse/cpu/minimizer/cg.h              (new Optimizer component)
src/cpu/{evaluate_driver,dynamics_driver,optimize_driver,evaluate,single_point}.cpp
src/cpu/scf/scf_solver.cpp
src/cpu/scf/internal/{hamiltonian_builder,eigen_solver,density_builder,
                      diis_mixer,dual_convergence_check}.cpp
src/cpu/scf/shell_list.cpp
src/cpu/minimizer/cg.cpp
tests/cpu/test_cpu_driver_role_validation.cpp
tests/cpu/test_evaluate_driver_capabilities.cpp
tests/cpu/test_optimize_driver_early_stop.cpp
tests/cpu/test_dynamics_full_budget.cpp
tests/cpu/test_minimize_driver_swap_optimizer.cpp
tests/cpu/scf/test_scf_solver.cpp
tests/cpu/scf/test_scf_solver_internals.cpp        (direct unit tests of internal helpers)
tests/cpu/scf/test_scf_solver_internal_isolation.cpp (asserts internals don't leak into facade)
tests/cpu/scf/test_scf_evaluator_no_alloc.cpp
tests/cpu/scf/test_integral_pack_no_alias.cpp
```

### Modified

```
include/atomiverse/facade/stage.h        (add evaluate::Phase = {prepare, evaluate, finalize}; mark scf::Phase legacy)
include/atomiverse/facade/driver.h       (Driver<>::run removed from public surface; bases drive it)
include/atomiverse/facade/driver_impl.h  (split termination-aware variants)
include/atomiverse/cpu/component.h       (add role() + capabilities() virtuals)
include/atomiverse/cpu/context.h         (IterationStatus outer_status; accessors)
include/atomiverse/cpu/md/md.h           (parent → DynamicsDriver)
include/atomiverse/cpu/md/{velocity_verlet,nose_hoover,fix_langevin,barostat,
                            lj,...}.h   (role()/capabilities() override)
include/atomiverse/cpu/minimizer/lbfgs.h (split: Minimize driver + LBFGS comp)
include/atomiverse/cpu/minimizer/init.h
include/atomiverse/cpu/scf/scf_evaluator.h   (no-alloc + holds SCFSolver + capabilities incl Force when gradient enabled)
include/atomiverse/cpu/scf/integrals.h   (drop shells/shell_data → ShellList)
include/atomiverse/cpu/scf/gradient.h    (GradientWorkspace + scf_gradient_inplace)
src/cpu/scf/scf_evaluator.cpp            (build/launch refactor; routes through SCFSolver)
src/cpu/scf/scf_solve.cpp                (≤ 30-line shim or deleted)
src/cpu/scf/native_integrals.cpp         (split → ≤ 4 files, each ≤ 600 lines)
src/cpu/scf/rys_roots.cpp                (split → ≤ 2 files)
src/cpu/scf/integral_derivatives.cpp     (split → ≤ 3 files)
src/cpu/scf/rys_eri.cpp                  (split → ≤ 2 files)
src/cpu/scf/xc_functional.cpp            (split → 3 files: lda, gga, dispatch)
src/cpu/driver.cpp                       (instantiate new bases; drop scf::Phase CPU instantiation)
src/cpu/CMakeLists.txt
tests/cpu/CMakeLists.txt
CLAUDE.md
.claude/NOTES.md
.claude/specs/INDEX.md
```

## Tasks

### Phase 1 — Three category bases + role/capability tagging (additive)

- [ ] **1.1** Add `ComponentRole` enum + `Capability` bitflags +
      bitwise operators in `facade/component_role.h`. Add `role()` and
      `capabilities()` virtuals on `cpu::ComponentBase`, defaulting to
      `None` / `None`.
- [ ] **1.2** Add `IterationStatus` in `facade/iteration_status.h`
      and an `outer_status` member + accessor on `cpu::Context`.
- [ ] **1.3** Add `evaluate::Phase = {prepare, evaluate, finalize, COUNT}`
      in `facade/stage.h`. Add `// CUDA-only legacy` comment on
      `scf::Phase`.
- [ ] **1.4** Add `cpu::CPUDriver<P>` template + the three
      concrete bases (`EvaluateDriver`, `DynamicsDriver`,
      `OptimizeDriver`). Each implements category-correct termination
      and validates required role + capability sets at `build()`.
- [ ] **1.5** Test: `tests/cpu/test_cpu_driver_role_validation.cpp`
      asserts:
      - `DynamicsDriver` without an Integrator → throw at build().
      - `DynamicsDriver` without any `Force`-capable evaluator → throw.
      - `OptimizeDriver` without an Optimizer → throw.
      - `OptimizeDriver` without any Energy/Force-capable evaluator → throw.
      - `EvaluateDriver` without any role==Evaluator AND capabilities!=None → throw.
      - All three with the right roles/caps build() successfully.
- [ ] **1.6** Test: `tests/cpu/test_evaluate_driver_capabilities.cpp`
      registers a stub energy-only evaluator (`role==Evaluator`,
      `capabilities==Energy`) on `EvaluateDriver` and asserts
      `build()` accepts it (no `Force` requirement).
- [ ] **1.7** Test: `tests/cpu/test_optimize_driver_early_stop.cpp` —
      stub optimizer + stub ConvergenceCheck that flips
      `outer_status.converged` on iter 3 → `run` returns 3.
- [ ] **1.8** Test: `tests/cpu/test_dynamics_full_budget.cpp` —
      stub ConvergenceCheck that always sets `converged=true` is
      ignored by Dynamics. A second case sets `stop_requested=true`
      and Dynamics breaks early (with reason recorded).

### Phase 2 — Migrate existing CPU drivers (parent change + algorithm-as-component)

- [ ] **2.1** `cpu::md::MD` parent: `Driver<md::Phase, cpu::device_tag>`
      → `DynamicsDriver`. No behavior change. Aspirin AIMD step within
      ±2% wall time.
- [ ] **2.2** Tag MD components:
      - `VelocityVerlet → Integrator`, `caps = None`.
      - `NoseHoover / FixLangevin → Thermostat`.
      - `Barostat / NoseHooverNPT → Barostat`.
      - `LJ → Evaluator`, `caps = Energy | Force | Stress`.
      - `HarmonicBond / HarmonicAngle / MultiHarmonicDihedral → Evaluator`,
        `caps = Energy | Force`.
      - `ZeroForces / ZeroEnergy / etc.` keep `role()==None`.
- [ ] **2.3** Rename `cpu::minimizer::LBFGS` (driver) →
      `cpu::minimizer::Minimize`. Add `[[deprecated]] using LBFGS = Minimize`.
- [ ] **2.4** Rename `cpu::minimizer::LBFGSStep` (component) →
      `cpu::minimizer::LBFGS`. Tag `role()==Optimizer`. Add
      `[[deprecated]] using LBFGSStep = LBFGS`.
- [ ] **2.5** Add `cpu::minimizer::CG` Optimizer component
      (Polak-Ribière, simplest correct implementation).
- [ ] **2.6** Test: `tests/cpu/test_minimize_driver_swap_optimizer.cpp` —
      same `Minimize` accepts either `LBFGS` or `CG`; both converge on
      H₂O / STO-3G to the same minimum within 1e-3 Ha.
- [ ] **2.7** Update existing tests to the new construction pattern;
      verify all minimizer / MD / AIMD tests stay green.

### Phase 3 — `cpu::scf::SCFSolver` (public low-level inner solver)

- [ ] **3.1** Carve `internal::HamiltonianBuilder`, `EigenSolver`,
      `DensityBuilder`, `DIISMixer`, `DualConvergenceCheck` out of the
      current `scf_solve.cpp` body, into `cpu/scf/internal/`. **Not**
      ComponentBase subclasses; **no** `default_bindings`. Each helper
      has its own `build()/step()` with pre-allocated workspaces.
- [ ] **3.2** Implement `cpu::scf::SCFSolver` — public class, public
      surface `{ctor, build, run, last_result}`. Internally drives
      the SCF cycle using the helpers from 3.1.
- [ ] **3.3** Test: `tests/cpu/scf/test_scf_solver.cpp` —
      `SCFSolver{default_config}.run(integrals)` reproduces legacy
      `scf_solve(integrals, default_config)` energy + density
      bit-for-bit (1e-12 Ha) on aspirin and ≥ 5 of the 12-benchmark
      molecules.
- [ ] **3.4** Test: `tests/cpu/scf/test_scf_solver_internals.cpp` —
      direct unit tests of `internal::DIISMixer`,
      `internal::DensityBuilder`, `internal::DualConvergenceCheck`
      via their public-by-namespace headers. (These are the tests the
      previous draft would have made impossible.)
- [ ] **3.5** Test: `tests/cpu/scf/test_scf_solver_internal_isolation.cpp` —
      asserts:
      - No internal helper inherits `cpu::ComponentBase`.
      - No internal helper overrides `default_bindings`.
      - No header outside `include/atomiverse/cpu/scf/internal/`
        includes any `internal/*.h` (grep gate).
- [ ] **3.6** Reduce `scf_solve()` to a ≤ 30-line shim. Update
      `rhf_solve` / `uhf_solve` / `uks_solve` similarly. Decision (delete
      shim vs keep) recorded in `.claude/NOTES.md`.

### Phase 4 — `SCFEvaluator` no-alloc + capabilities-based registration

- [ ] **4.1** Test (RED): `tests/cpu/scf/test_scf_evaluator_no_alloc.cpp`
      mirrors `tests/cpu/test_md_no_alloc.cpp`. Asserts steady-state
      launches allocate zero times. Skips on macOS.
- [ ] **4.2** Add `rebuild_integrals_inplace(frame, basis, charge, pack)`.
- [ ] **4.3** Add `GradientWorkspace` + `scf_gradient_inplace(frame,
      ints, result, config, ws, out)`.
- [ ] **4.4** Refactor `SCFEvaluator`: hoist `IntegralPack`,
      `GradientWorkspace`, density-guess buffers, ASPC scratch, and
      the `SCFSolver` member into `build()`. `launch()` becomes
      in-place updates + `solver_.run(integrals_)`. Tag
      `role()==Evaluator`. `capabilities()` returns
      `Energy | DensityMatrix | Orbitals` always, plus `Force` when
      `with_gradient(true)` was set (default true for backwards-compat;
      tests cover both).
- [ ] **4.5** Test 4.1 turns GREEN. All AIMD / SCF regression tests pass.
- [ ] **4.6** Bench: aspirin AIMD step within ±5% of pre-Phase-4 baseline.

### Phase 5 — `cpu::Evaluate` driver + IntegralPack split + file granularity

- [ ] **5.1** Implement `cpu::Evaluate : EvaluateDriver`. Default
      infrastructure: `ZeroForces`, `ZeroEnergy` (and stress if
      `VirialStress` is added). One-pass `run(frame)` API.
- [ ] **5.2** Optional helper: `cpu::single_point(frame, evaluator)` as a
      free function (not class) that constructs a transient
      `Evaluate`. Kept thin; documented as syntactic sugar.
- [ ] **5.3** Test: `tests/cpu/test_evaluate_driver.cpp` —
      - Single-shot LJ force eval matches a 0-step `cpu::md::MD`
        invocation bit-for-bit.
      - Single-shot SCF energy via `Evaluate + SCFEvaluator(.with_gradient(false))`
        matches `SCFSolver.run(integrals)` to 1e-12 Ha on aspirin and
        5 benchmark molecules.
      - Energy-only evaluator (no Force capability) is accepted by
        `Evaluate` and rejected by `Dynamics`.
- [ ] **5.4** Test: `tests/cpu/scf/test_integral_pack_no_alias.cpp` (RED).
- [ ] **5.5** Introduce `cpu::scf::ShellList` (header + impl + free
      function `build_shell_list(atoms, basis)`).
- [ ] **5.6** Move `shells` + `shell_data` out of `IntegralPack` into
      `ShellList`. Update consumers. Test 5.4 GREEN.
- [ ] **5.7** Split `native_integrals.cpp` (1710) → ≤ 4 files.
- [ ] **5.8** Split `rys_roots.cpp` (1455) → ≤ 2 files.
- [ ] **5.9** Split `integral_derivatives.cpp` (1146) → ≤ 3 files.
- [ ] **5.10** Split `rys_eri.cpp` (1043) → ≤ 2 files.
- [ ] **5.11** Split `xc_functional.cpp` (639) → 3 files
       (`xc_lda.cpp`, `xc_gga.cpp`, `xc_dispatch.cpp`).

### Phase 6 — Cleanup, docs, deprecation removal

- [ ] **6.1** Remove `[[deprecated]] using LBFGS = Minimize` and
      `using LBFGSStep = LBFGS`. Remove `scf_solve()` shim if no
      remaining caller. Update affected docstrings.
- [ ] **6.2** Update `CLAUDE.md`:
      - "Driver taxonomy" subsection: three current outer categories +
        SCF-as-inner-solver. Explicit note that future categories are
        permitted only with distinct outer-loop semantics.
      - Updated Drivers table.
      - A paragraph: "Single-point energy is `Evaluate + SCFEvaluator`.
        There is no `SinglePoint` class. SCF iteration is encapsulated
        by `SCFSolver`, which is a public low-level API used internally
        by `SCFEvaluator`."
      - `ComponentRole` + `Capability` in Conventions; LBFGS/CG/FIRE
        as Optimizer components, etc.
- [ ] **6.3** Append `.claude/NOTES.md` entry: three-category
      taxonomy, SCF-as-inner-solver, deprecation timeline, CUDA
      follow-up.
- [ ] **6.4** Confirm `cpu-scf-driver-alignment.md` and
      `phase2-singlepoint-driver.md` already marked superseded.
- [ ] **6.5** Open follow-up `cuda-driver-taxonomy.md` (status: planned).

## Testing strategy

- **Role + capability validation**
  (`test_cpu_driver_role_validation.cpp`,
  `test_evaluate_driver_capabilities.cpp`).
- **Termination semantics** (`test_optimize_driver_early_stop.cpp`,
  `test_dynamics_full_budget.cpp`).
- **Evaluate one-pass parity** (`test_evaluate_driver.cpp`).
- **No-alloc** (`test_scf_evaluator_no_alloc.cpp`) — Linux/glibc only.
- **Inner solver parity** (`test_scf_solver.cpp`) — bit-for-bit vs
  `scf_solve` on aspirin + 5 benchmarks.
- **Inner solver unit tests** (`test_scf_solver_internals.cpp`) — direct
  tests of DIISMixer, DensityBuilder, DualConvergenceCheck.
- **Inner-solver isolation** (`test_scf_solver_internal_isolation.cpp`)
  — no `ComponentBase` inheritance, no facade-public include of
  `internal/*.h`. (Tests can include freely.)
- **Pointer discipline** (`test_integral_pack_no_alias.cpp`).
- **Optimizer hot-swap** (`test_minimize_driver_swap_optimizer.cpp`).
- **Regression**: every existing test stays green.
- **Domain**: PySCF parity for `Evaluate + SCFEvaluator` at PBE/def2-SVP
  and HF/STO-3G; libcint parity for analytic gradient.
- **Architecture grep gates** (acceptance):
  - No `class\s+\w+\s*:\s*public\s+Driver<` outside the three category
    base files.
  - No `SCFDriver`, `SinglePoint`, `LangevinDriver`, `NPTDriver`,
    `REMDDriver`, `Workflow`, `Pipeline`, `Composite`,
    `OptimizeThenSinglePoint` class names anywhere in CPU code.

## Risks

- **Scope is wide.** Six phases; each independently shippable.
- **`evaluate::Phase` over-design.** The smallest sensible shape is
  two phases (`evaluate, finalize`); spec proposes three for
  ZeroForces/VirialStress/Recorder placement. Drop to two at
  implementation time if three feels excessive.
- **SCFSolver vs scf_solve performance.** Phase 4.6 bench gate ≤ 5%
  regression; worst case, inline aggressively.
- **`IterationStatus` mis-write hazard.** Mitigation: only
  ConvergenceCheck-tagged components have a documented contract that
  reads/writes `outer_status`; SCF internals use their own internal
  state, never the outer one. A test in Phase 1.8 covers the
  separation.
- **CUDA divergence.** CUDA keeps current shape until the follow-up
  spec.
- **Mental model risk.** "SCFSolver vs SCFEvaluator" is a meaningful
  distinction (public inner solver vs outer evaluator component).
  Mitigation: CLAUDE.md and the SCFSolver header open with a
  paragraph naming both and saying which is which.
- **Capability bitflag growth.** New capabilities (Hessian, NACME,
  RDM-2) will arrive over time. Reserved bits up to 32; if it ever
  overflows, switch to `std::bitset<N>`. Out of scope today.

## Exit gate

All Phase 1-6 tasks ticked + `cpu-driver-taxonomy.acceptance.md`
criteria all green + CLAUDE.md updated + predecessors marked
superseded + INDEX entry updated + user sign-off recorded.
