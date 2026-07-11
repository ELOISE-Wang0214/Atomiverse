# Chapter 12: Geometry Minimization with L-BFGS or CG

The `cpu::minimizer::Minimize` driver answers a different question from
the MD chapters. Instead of "how does the system evolve in time?" it
asks "where is the nearest local minimum of the potential energy
surface?" The wiring is closely related — same `cpu::Frame`, same
`SCFEvaluator` at the force phase — but the outer loop is no longer a
time integrator. It is a quasi-Newton or conjugate-gradient *search*,
and it terminates when the gradient and the displacements between
successive geometries fall below tolerance.

This chapter walks through the full geometry-optimization recipe with
an aspirin (acetylsalicylic acid) molecule at RHF / STO-3G. The example
to read alongside this chapter is
[`examples/05_geom_opt.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/05_geom_opt.cpp).

## What the Minimize Driver Actually Does

`cpu::minimizer::Minimize` is a thin subclass of
[`cpu::OptimizeDriver`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/driver/optimize_driver.h).
The base class owns the loop:

```text
loop until converged or max_iter reached:
    outer_status.reset()                 // clear convergence flag
    run one pass of the compiled plan    // search_direction → line_search →
                                          // force → convergence_check → end_of_step
    if outer_status.converged:  break    // normal success
    if outer_status.stop_requested: break // explicit early stop
```

The driver does not know about L-BFGS, conjugate-gradient, FIRE, or any
specific algorithm. All it knows is to break the loop when an
`outer_status.converged` flag has been set. The flag is set by whichever
component is bound at `convergence_check` — by convention, the optimizer
component itself.

Build-time validation at `Minimize::build()` checks that:

1. At least one component has `role() == ComponentRole::Optimizer`.
2. At least one component has `role() == ComponentRole::Evaluator` and
   advertises `Capability::Energy` *or* `Capability::Force`.

If either is missing, `build()` throws with an explanatory message.
Unlike MD, the optimizer side has only one composition rule — pick
*one* optimizer component. The two shipped optimizers, `LBFGS` and
`CG`, are mutually exclusive in the same way Nose-Hoover and Langevin
are.

## The Phase Order

```text
search_direction → line_search → force → convergence_check → end_of_step
```

The `Minimize` driver auto-adds two infrastructure components in its
constructor:

- `cpu::minimizer::ZeroForces` at `force @ -1000`
- `cpu::minimizer::ZeroEnergy` at `force @ -999`
- `cpu::SnapshotPublisher` at `end_of_step @ +10`

Note that `force` is the same logical phase as in MD — same evaluator
components plug into both drivers. The `cpu::minimizer::ZeroForces`
component is a separate class from `cpu::md::ZeroForces` only because
its default binding is `minimizer::Phase::force`, not `md::Phase::pre_force`.

The optimizer component is bound at *two* phases: `search_direction` (it
computes a step) and `convergence_check` (it publishes the verdict).
That is why a single `LBFGS{...}` registration is enough for the entire
algorithm.

## L-BFGS, the Default Optimizer

`cpu::minimizer::LBFGS` is the limited-memory BFGS quasi-Newton
optimizer (Nocedal, *Math. Comp.* 35, 773 (1980)). At every iteration
it stores up to $m$ pairs of position and gradient differences,

$$s_k = x_{k+1} - x_k, \qquad y_k = g_{k+1} - g_k,$$

and uses Nocedal's two-loop recursion to compute a search direction
that approximates the Newton step

$$d_k = -H_k g_k$$

without ever materializing $H_k$. Step zero falls back to steepest
descent. From step one onward the L-BFGS direction is used.

The constructor takes:

```cpp
cpu::minimizer::LBFGS{
    /*memory=*/   5,         // history size m
    /*alpha0=*/   0.01,      // initial step length [Å²/eV]
    /*ftol=*/     0.05,      // max-force tol      [eV/Å]
    /*max_ls=*/   10,        // backtracking attempts (only used when line search is on)
    /*rms_ftol=*/ 0.03,      // rms-force tol      [eV/Å]
    /*max_dtol=*/ 0.002,     // max-disp tol       [Å]
    /*rms_dtol=*/ 0.001,     // rms-disp tol       [Å]
    /*etol=*/     1e-5       // |ΔE| tol           [eV]
};
```

Convergence requires *all five* criteria — max force, RMS force, max
displacement, RMS displacement, and energy change — to be below their
tolerances simultaneously. That is stricter than a single-criterion
test (e.g. "max force < 0.05 eV/Å") and avoids the failure mode where a
flat saddle with a tiny gradient passes but the geometry is still
walking around. If you find yourself loosening tolerances to make a run
succeed, prefer to loosen them all uniformly rather than disable any.

By default `LBFGS` uses a fixed step length `alpha0`. You can opt into
backtracking line search:

```cpp
auto step = cpu::minimizer::LBFGS{m, alpha0, ftol, max_ls,
                                   rms_ftol, max_dtol, rms_dtol, etol}
              .line_search(/*c1=*/1e-4, /*contract=*/0.5,
                           /*min_alpha=*/1e-6);
opt.add(std::move(step));
```

When the line search is on, every step that fails the Armijo condition

$$E_{k+1} - E_k \le c_1 \alpha_k g_k \cdot d_k$$

rolls back: positions return to $x_k$, $\alpha_k$ is contracted by
`contract`, and the trial is re-taken from $x_k$ along the same
direction. The $(s, y)$ history is updated only on accepted trials.
`min_alpha` floors the contraction so the loop cannot run forever on a
pathological evaluation.

## Conjugate Gradient as an Alternative

`cpu::minimizer::CG` is a Polak-Ribière+ conjugate-gradient optimizer.
It uses the same convergence channel and the same five-criterion rule
as `LBFGS`. The only thing that changes is the direction:

$$\beta_k = \max\!\left(0,\; \frac{g_k \cdot (g_k - g_{k-1})}{g_{k-1} \cdot g_{k-1}}\right), \qquad d_k = -g_k + \beta_k d_{k-1}.$$

The non-negativity projection on $\beta$ already restarts the algorithm
when curvature information becomes ill-conditioned, so there is no
periodic-restart parameter.

To swap CG in for L-BFGS, change one line:

```cpp
opt.add(cpu::minimizer::CG{/*alpha0=*/0.01, /*ftol=*/0.05,
                           /*rms_ftol=*/0.03, /*max_dtol=*/0.002,
                           /*rms_dtol=*/0.001, /*etol=*/1e-5});
```

The rest of the driver setup is unchanged. The hot-swap is exercised in
[`tests/cpu/test_minimize_driver_swap_optimizer.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/tests/cpu/test_minimize_driver_swap_optimizer.cpp),
which builds the same molecule twice and checks that both optimizers
report convergence under the same tolerances.

In practice, L-BFGS converges in fewer iterations than CG on
well-conditioned problems (rigid molecules near minimum), and the
two are roughly comparable when the surface is flat or the gradient is
noisy. CG's smaller per-step memory makes it the right default for
very large systems where storing $2m$ vectors of length $3N$ would be
expensive.

## Putting It Together: Aspirin at RHF/STO-3G

Aspirin is 21 atoms (9 C + 8 H + 4 O), 73 STO-3G basis functions, 94
electrons. It is small enough that an SCF + analytic gradient finishes
in a fraction of a second per step, but large enough that you actually
see L-BFGS doing useful work.

The setup, in full:

```cpp
auto frame = make_aspirin();              // see examples/05_geom_opt.cpp

auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
auto config    = cpu::scf::rhf_config(integrals);   // RHF: no XC grid, fast
                                                    // (default_config => KS/PBE)
config.tol     = 1e-6;                    // Hartree
config.max_iter = 200;
config.force_native_two_electron_deriv = true;

cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{
    5,        // memory
    0.01,     // alpha0  [Å²/eV]
    0.05,     // ftol    [eV/Å]
    10,       // max_ls
    0.03,     // rms_ftol
    0.1,      // max_dtol  [Å]   — slightly looser; fixed-step L-BFGS takes ~0.05–0.2 Å steps
    0.05,     // rms_dtol  [Å]
    1e-3      // etol      [eV]
});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.build();

int iters = opt.run(frame, /*max_iter=*/100);
```

The single line that does *all* the chemistry is
`opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config})`.
The evaluator's `launch()` rebuilds the integrals at the current
geometry, runs the SCF cycle, computes the analytic nuclear gradient
($\partial E / \partial R_A$), and writes the negative gradient into
`frame.atoms.fx/fy/fz`. The optimizer then reads those forces from the
frame in its `convergence_check` pass, decides the new geometry, and
writes the updated coordinates back.

`Minimize::run` returns the number of iterations actually executed
(≤ `max_iter`). On a converged run that number is your iteration count.
On a non-converged run that number equals `max_iter`, and you can
inspect `LBFGS::trace()` for the per-step record.

## What Loose Tolerances Actually Cost You

The aspirin example above uses `max_dtol = 0.1 Å` — quite loose by
quantum-chemistry standards. The reason is that the L-BFGS line search
is opt-in, and a fixed-step optimizer near a minimum routinely takes
0.05 Å — 0.2 Å steps when the gradient is small but non-vanishing.
Setting `max_dtol = 0.001 Å` on a fixed-step run will keep the loop
running until `max_iter`, even when the energy and forces are clearly
at minimum.

Two clean cures, in order of effort:

1. **Match displacement tolerance to step size.** If your `alpha0` is
   $0.01\,\text{Å}^2/\text{eV}$ and the typical $|g|$ near minimum is
   $0.05\,\text{eV/Å}$, the predicted step is ~$5 \times 10^{-4}\,\text{Å}$,
   so `max_dtol = 1\text{e-}3` is appropriate. Loosen it only if you
   are far from this regime.
2. **Turn on the line search.** With Armijo backtracking the optimizer
   refuses to take steps that increase the energy, and the typical step
   length adapts to the local curvature. Tight displacement tolerances
   become physically meaningful.

The five-criterion rule means you should report all five thresholds
when you publish a geometry, not just the force tolerance.

## Recorders for the Optimization Trace

Three recorder components live under
[`include/atomiverse/cpu/minimizer/`](https://github.com/MolCrafts/Atomiverse/blob/dev/include/atomiverse/cpu/minimizer/)
and are designed specifically for optimization runs.

`OptLog` writes a tab-separated row per step with `step`, `E [eV]`,
`max_force [eV/Å]`, `rms_force [eV/Å]`, and `step_norm [Å]`. It
defaults to binding at `end_of_step @ 0` and works against any force
evaluator (classical or SCF). Useful as the canonical
"is-this-converging" log:

```cpp
opt.on(minimizer::Phase::end_of_step)
   .add(cpu::minimizer::OptLog{"aspirin_opt.tsv", /*also_stdout=*/true});
```

`FinalGeometry` writes a single-frame XYZ file exactly once, at the end
of the run. After the run it holds the last-visited geometry — converged
or not. Cheap to enable and ideal for "I just want the final
coordinates":

```cpp
opt.on(minimizer::Phase::end_of_step)
   .add(cpu::minimizer::FinalGeometry{"aspirin_final.xyz"});
```

`TrajectoryWriter` appends a frame per step to an XYZ trajectory.
Useful for visualizing the optimization path or re-running an
animation. Throttle long runs with the `stride` parameter:

```cpp
opt.on(minimizer::Phase::end_of_step)
   .add(cpu::minimizer::TrajectoryWriter{"aspirin_traj.xyz", /*stride=*/1});
```

All three are MPI-aware (rank-0-only writes) and respect the same
build/launch/destroy lifecycle as every other component.

## When Optimization Diverges

Geometry optimization with quantum forces fails for two qualitatively
different reasons. Distinguishing them quickly is what separates a
fifteen-minute fix from a two-day debug.

**The optimizer is the problem** when the energy is converging
monotonically and the forces are shrinking, but the loop hits
`max_iter`. Symptoms in `LBFGS::trace()`: `max_force` plateaus near
the threshold and oscillates; `step_norm` is large because the fixed
step keeps overshooting. Cure: turn on the line search, or tighten
`alpha0`.

**The SCF is the problem** when the energy jumps non-monotonically,
the forces blow up, or the atoms move by tens of Ångström in one step.
Symptoms: SCFEvaluator's `last_result().converged` is false; the
trajectory shows atoms passing through each other. Cure: this is not
an optimizer bug. Either the SCF tolerance is set too loose
(`config.tol > 1e-6` produces noisy gradients near minimum), the
initial guess is poor (try MINAO via `default_config`), or you are
running the wrong method for the chemistry (e.g. RHF on a singlet
biradical, or a small basis on an anion). The
[SCF tutorials](ch02-rhf.md) cover these cases.

If forces, energy, and SCF are all healthy but the geometry will not
move, you may be at a saddle point. None of the shipped optimizers
follow modes; they all minimize. Use a transition-state search package
or perturb the starting geometry along a guessed mode and minimize from
there.

## What's Next

The exact same recipe — `Minimize` driver, `SCFEvaluator` at force
phase, optimizer component at `search_direction` — also drives
classical relaxation. Replace the SCFEvaluator with `cpu::PairLJCut` and
`cpu::BondHarmonic`, and you have a force-field minimization.
The optimizer does not know the difference, and the convergence
tolerances stay in the same units.

If you came here from [Chapter 15: AIMD](ch15-cross-method.md), the symmetry is
worth pausing on. The same `SCFEvaluator{config}` instance plugs into
either driver:

```cpp
// Geometry optimization (Optimize loop)
cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{...});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{cfg});

// Ab initio molecular dynamics (Dynamics loop)
cpu::md::MD                sim;
sim.add(cpu::md::VelocityVerlet{dt});
sim.add(cpu::md::NoseHoover{T, tau, dt});
sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{cfg}.predictor(3));
```

Two drivers, one evaluator, two ensembles. That is the architecture
paying off.
