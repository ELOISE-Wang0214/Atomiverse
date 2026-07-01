# Chapter 15: Cross-Method Workflows — Swap the Evaluator

Every chapter of this book has quietly repeated one sentence: *change method
by swapping the evaluator, keep the driver.* This closing chapter stops
repeating it and instead runs it, end to end. We take one molecular system
and drive it through single-point energy, geometry optimization, and
molecular dynamics — first with a cheap classical force field, then with a
quantum-mechanical SCF force provider — changing, in each case, exactly one
line. The point is not any single result; it is that the architecture makes
"cheap classical" and "expensive quantum" the *same code with one line
different*, which is what lets you develop a workflow on fast forces and
promote it to *ab initio* accuracy when you are ready.

The six questions:

1. **The concept.** What "method" and "cost vs accuracy" mean once the loop
   is factored out.
2. **The Atomiverse representation.** The `Evaluator` component as the single
   pluggable degree of freedom.
3. **The C++ API.** `cpu::single_point`, `cpu::Evaluate`, `cpu::md::MD`,
   `cpu::minimizer::Minimize`, and the `SCFEvaluator` knobs.
4. **The expected result.** Matching structure, diverging cost, between
   classical and quantum runs.
5. **The validation.** A golden end-to-end run and the composition tests that
   certify the evaluator coexists with every driver.
6. **The recurring pattern.** This *is* the recurring pattern — the spine of
   the whole engine.

## The Concept: Method Is a Force Provider, Not a Program

In most simulation software, "run a classical MD" and "run an *ab initio*
MD" are different programs with different input decks, different data
structures, and different mental models. Atomiverse rejects that split. It
observes that a time integrator, a geometry optimizer, and a single-point
evaluator do not actually *care* where forces come from — they only care that
*something* fills `frame.atoms.fx/fy/fz` and `frame.scalars.energy` at the
right phase. So "method" is factored down to a single replaceable object: the
**Evaluator**.

Once method is a force provider, **cost versus accuracy** becomes a property
of that one object, cleanly separated from everything else:

- A **classical** evaluator (`PairLJCut` + bonded terms) costs $O(N)$ per
  step with a small constant. It is accurate only insofar as its
  parameters were fit to the chemistry you are studying.
- A **quantum** evaluator (`SCFEvaluator`) costs a full self-consistent-field
  solve per step — orders of magnitude more — but is *transferable*: it needs
  no per-system parameters and describes bond breaking, polarization, and
  charge transfer that no fixed force field captures.

The architecture's payoff is that you choose that trade-off with a single
constructor, at the last possible moment, without touching the workflow
around it.

## The Representation: One Pluggable Degree of Freedom

The three CPU driver categories share the `force` phase as a common socket:

```text
Evaluate     : prepare → evaluate → finalize
Dynamics (MD): initial_integrate → pre_force → force → post_force → final_integrate → end_of_step
Optimize     : search_direction → line_search → force → convergence_check → end_of_step
Response     : setup → displace → force → accumulate → diagonalize → thermochem → end_of_step
```

Any component with `role() == Evaluator` and `Capability::Force` (or, for
single points, `Capability::Energy`) drops into that socket on *any* of the
four drivers. The driver validates the socket is filled and sequences it; it
never inspects what the evaluator *is*. That single indirection is the whole
trick, and it is why this book could teach RHF, UHF, DFT, hybrids, plane
waves, and classical force fields as *variations on one loop* rather than as
separate subsystems.

## The C++ API: One System, Every Driver, Two Evaluators

Take a molecular system in a `cpu::Frame` and walk it through the three
workflows. Keep the evaluator abstract for a moment by defining two:

```cpp
using namespace ATV;

// Classical force provider (parameters are illustrative)
auto make_classical = []{ return cpu::PairLJCut{0.0104, 3.4, 6.8}; };

// Quantum force provider: RHF/STO-3G, forces on by default
auto make_quantum = [](const cpu::scf::IntegralPack& ints){
    return cpu::scf::SCFEvaluator{cpu::scf::rhf_config(ints)};
};
```

### Single-point energy — `cpu::Evaluate` / `cpu::single_point`

The `cpu::single_point` helper wraps the one-pass `cpu::Evaluate` driver
around a single evaluator. For an energy-only quantum single point, turn the
gradient off so the evaluator drops `Capability::Force` and skips the gradient
work:

```cpp
auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
auto scf  = cpu::scf::SCFEvaluator{cpu::scf::rhf_config(ints)};
scf.with_gradient(false);                 // energy only

cpu::single_point(frame, scf);            // frame.scalars.energy [eV] populated
double e_ha = scf.last_result().energy;   // and in Hartree on the result
```

For a classical single point, drive `cpu::Evaluate` directly with the
categorical setters (a pair potential needs a neighbor list):

```cpp
cpu::Evaluate eval;
eval.set_pair(make_classical());
eval.set_neighbor(cpu::md::NeighborList<>{6.8, 0.3});
eval.build();
eval.run(frame);                          // frame.scalars.energy [eV]
```

### Geometry optimization — `cpu::minimizer::Minimize`

The optimizer driver takes an `Optimizer` component (LBFGS/CG) plus *any*
force evaluator at the `force` phase. Swapping classical for quantum is the
one-line change:

```cpp
cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{5, 0.01, 0.05, 10, 0.03, 0.1, 0.05, 1e-3});

// Classical relaxation:
opt.on(minimizer::Phase::force).add(make_classical());
//   ── OR ── ab initio geometry optimization (one line different):
opt.on(minimizer::Phase::force).add(make_quantum(ints));

opt.build();
int iters = opt.run(frame, /*max_iter=*/100);
```

### Molecular dynamics — `cpu::md::MD`

Identically, classical MD and AIMD differ only in the force-phase component.
The integrator, thermostat, and observers are untouched:

```cpp
cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{/*dt=*/0.5});
sim.set_thermostat(cpu::md::NoseHoover{300.0, 100.0, 0.5});   // NVT

// Classical MD:
sim.set_pair(make_classical());
sim.set_neighbor(cpu::md::NeighborList<>{6.8, 0.3});
//   ── OR ── AIMD (swap the force provider):
sim.on(md::Phase::force).add(make_quantum(ints).predictor(3));

sim.set_thermo(100);
sim.build();
sim.run(frame, n_steps);
```

The two `SCFEvaluator` knobs that appear here are the workflow-tuning levers
this chapter is really about:

- **`with_gradient(bool)`** — on by default (writes forces, advertises
  `Capability::Force`); turn it off for energy-only single points, where it
  also drops the force capability so the evaluator is legal in a pure
  `Evaluate` context.
- **`predictor(order)`** — the ASPC density predictor-corrector. For MD it
  extrapolates the density matrix from previous steps, cutting the SCF
  iteration count by 30–50% on production-scale systems (see
  [Tutorial 6](06-aimd.md)). It is a *dynamics* optimization — meaningful only
  when the geometry changes smoothly step to step, which is exactly the MD
  case.

## The Expected Result: Same Structure, Divergent Cost

Run the *same* molecule through classical and quantum optimization and you
should see the two hallmarks of this architecture:

- **The workflows are structurally identical.** Both optimizations converge
  by the same five-criterion rule, return the same `int` iteration count
  contract, and leave the relaxed geometry on the same `cpu::Frame`. Both MD
  runs report the same `step / PE / KE / T` thermo lines through the same
  `cpu::Thermo` component. Nothing downstream of the force phase can tell the
  two apart.
- **The cost is worlds apart, and the accuracy tracks it.** The classical
  step is microseconds; the quantum step is a full SCF. In exchange, the
  quantum run gets the physics right where the force field's parameters were
  never fit — the canonical example being bond dissociation, where a
  fixed-form classical potential cannot break a bond and an SCF (with an
  appropriate method) can.

The disciplined way to report a cross-method comparison is to name each
number's class. "Classical and quantum optima differ by X Å" is a
**physical-approximation** statement about the two *methods*. "The AIMD total
energy drifts by less than Y eV/atom/ps under NVE" is a
**numerical-tolerance** statement about the *integrator*, independent of which
evaluator produced the forces. And "a fresh single point on a mid-trajectory
snapshot reproduces the propagation force to ~1e-5 eV/Å" is the internal
consistency check that the evaluator carries no illegal step-to-step state.

## The Validation: A Golden Run and Composition Tests

The cross-method promise is validated at two levels.

**Composition.** The `SCFEvaluator` is tested to coexist, without conflict,
with every ensemble and driver it is supposed to plug into — Nosé-Hoover,
Langevin, and Berendsen-barostat MD; L-BFGS and CG minimization; and the
single-pass evaluate driver. These composition tests are what certify that
"swap the evaluator" is literally true and not just aspirational: the same
evaluator instance produces correct forces regardless of the loop it is
embedded in. A specific state-leakage regression confirms that at any point
in an ASPC-warmstarted run, a stand-alone single point on the snapshot
geometry agrees with the propagation force to ~1e-5 eV/Å — catching the one
failure mode (cached geometry-dependent state) that a pluggable evaluator
could introduce.

**End-to-end golden run.** The designated cross-method validation is a full
optimize-then-propagate workflow run once as a golden reference, checking that
the pipeline — build integrals, minimize, run dynamics, read observables —
produces the committed energies and geometry within tolerance. This is the
integration test that the *whole* architecture, not any one component, does
what the manual claims.

Underneath, each half of the swap stands on its own validated foundation: the
classical forces against LAMMPS ([Tutorial 7](07-classical-md.md)), the
quantum forces against PySCF (Chapters [2](02-rhf.md), [6](ch06-hybrid-dft.md),
[14](ch14-vibrations.md)). Cross-method validation does not re-derive those; it
certifies that *composing* them changes nothing about either.

## The Recurring Pattern: This Was the Pattern All Along

There is no new pattern to introduce here, because this chapter *is* the
pattern the whole book has been building toward. Every earlier "how the
pattern recurs" section pointed at this one:

```text
cpu::Evaluate            + PairLJCut     = classical single point
cpu::Evaluate            + SCFEvaluator  = quantum single point
cpu::minimizer::Minimize + PairLJCut     = classical relaxation
cpu::minimizer::Minimize + SCFEvaluator  = ab initio geometry optimization
cpu::md::MD              + PairLJCut     = classical molecular dynamics
cpu::md::MD              + SCFEvaluator  = ab initio molecular dynamics
cpu::vib::Hessian        + PairLJCut     = force-field vibrational analysis
cpu::vib::Hessian        + SCFEvaluator  = ab initio vibrational analysis
```

Four drivers down the left, two evaluators across the top, eight distinct
computational-chemistry workflows in the grid — and every cell is *real*,
buildable code that shares one `cpu::Frame`, one component contract, and one
set of units. Adding a ninth workflow does not mean writing a ninth program;
it means writing one new evaluator (a new functional, a new force field, a
machine-learned potential) and dropping it into the socket, or adding one new
driver category (a Monte-Carlo sampler, a nudged-elastic-band path) that
reuses every evaluator that already exists.

That factorization — *one loop, many physics* — is the whole thesis of
Atomiverse, and the reason a single manual could teach classical mechanics
and quantum chemistry as one subject. You have now seen it hold across every
method the engine supports. The next method you reach for will be one line
different from the one you already know.

## Validation Reference

The evaluator-composition tests (SCFEvaluator coexisting with NVE/NVT/NPT
integrators and both minimizers, plus the state-leakage single-point-vs-
propagation check) live in
[`tests/cpu/scf/test_scf_evaluator.cpp`](../../tests/cpu/scf/test_scf_evaluator.cpp)
and
[`tests/cpu/scf/test_aimd_conservation.cpp`](../../tests/cpu/scf/test_aimd_conservation.cpp);
the optimizer hot-swap (LBFGS ↔ CG on one `Minimize` driver) in
[`tests/cpu/test_minimize_driver_swap_optimizer.cpp`](../../tests/cpu/test_minimize_driver_swap_optimizer.cpp).
The end-to-end golden cross-method run is the designated integration
validation for the "swap the evaluator" workflow.
