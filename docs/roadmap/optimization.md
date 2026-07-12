# Geometry Optimization and Reaction Paths

Legend and caveats: [Roadmap index](index.md).

The Optimize driver (`cpu::minimizer::Minimize`, `cuda::minimizer::Minimize`)
hosts the optimizer as a hot-swappable Component, so any force-providing
evaluator — LJ, bonded, SCF, GPW — can be relaxed by any minimizer. What is
missing is breadth of algorithm (no FIRE, no quasi-Newton beyond L-BFGS) and,
more consequentially, **any notion of relaxing the cell**.

Transition states and reaction paths are entirely absent. NEB in particular
would need a new Driver category (a coupled chain of images terminating when
every image converges), which the architecture anticipates but does not
implement.

## 1. Geometry Optimization

| Method                                     | CPU | CUDA |
| ------------------------------------------ | :-: | :--: |
| Steepest descent                           |  ○  |   ✓  |
| Fletcher–Reeves conjugate gradient|  ⊘  |   ⊘  |
| Polak–Ribière conjugate gradient           |  ✓  |   ✓  |
| Nonlinear conjugate gradient               |  ✓  |   ✓  |
| BFGS|  ⊘  |   ⊘  |
| L-BFGS                                     |  ✓  |   ✓  |
| Trust-region optimization                  |  ○  |   ○  |
| Newton–Raphson optimization|  ⊘  |   ⊘  |
| Rational-function optimization             |  ○  |   ○  |
| Partitioned rational-function optimization |  ○  |   ○  |
| Augmented-Hessian optimization             |  ○  |   ○  |
| FIRE minimization                          |  ○  |   ○  |
| Quick-min minimization|  ⊘  |   ⊘  |
| Damped-dynamics minimization|  ⊘  |   ⊘  |
| Cartesian-coordinate optimization          |  ✓  |   ✓  |
| Internal-coordinate optimization           |  ○  |   ○  |
| Redundant-internal-coordinate optimization |  ○  |   ○  |
| Constrained geometry optimization          |  ○  |   ○  |
| Variable-cell optimization                 |  ○  |   ○  |
| Joint atomic and cell optimization         |  ○  |   ○  |

Notes:

- **Steepest descent is CUDA-only** as a distinct component (`cuda/minimizer/sd.cuh`).
  On CPU, SD exists only as the degenerate β = 0 limit of CG — there is no SD
  component to select.
- The CG flavour is **Polak–Ribière⁺** (restart by β clamp) on both backends;
  Fletcher–Reeves is not implemented. It is a nonlinear CG, hence both rows ✓.
- Line search on CPU offers fixed-step, Armijo backtracking, and **strong Wolfe**
  (Nocedal & Wright bracket + zoom). The CUDA L-BFGS uses Armijo backtracking.
- Convergence on CUDA checks five criteria (max/RMS force, max/RMS displacement,
  ΔE); the CUDA SD component deliberately runs its full budget without one.
- **Every minimizer moves atoms only.** No cell degrees of freedom, no
  constraints, no internal coordinates. Combined with the missing periodic stress
  tensor (see [Periodic DFT](dft.md#11-periodic-response-and-analysis)),
  variable-cell relaxation is doubly
  blocked: no stress to relax against, and no optimizer that would accept it.

## 2. Transition-State and Reaction-Path Methods

| Method                                        | CPU | CUDA |
| --------------------------------------------- | :-: | :--: |
| Eigenvector-following transition-state search |  ○  |   ○  |
| Partitioned-RFO transition-state search       |  ○  |   ○  |
| Dimer method                                  |  ○  |   ○  |
| Nudged elastic band|  ⊘  |   ⊘  |
| Climbing-image nudged elastic band            |  ○  |   ○  |
| Adaptive nudged elastic band                  |  ○  |   ○  |
| String method|  ⊘  |   ⊘  |
| Growing string method                         |  ○  |   ○  |
| Freezing string method                        |  ○  |   ○  |
| Intrinsic reaction coordinate                 |  ○  |   ○  |
| Minimum-energy path optimization              |  ○  |   ○  |
| Minimum-energy crossing point search          |  ○  |   ○  |
| Conical-intersection optimization             |  ○  |   ○  |
| Basin hopping                                 |  ○  |   ○  |
| Minima hopping                                |  ○  |   ○  |

## Superseded (⊘)

| Method | Dominated by |
| ------ | ------------ |
| BFGS | L-BFGS (implemented) — same superlinear convergence, O(mN) memory instead of O(N²). The dense inverse Hessian buys nothing at simulation sizes. |
| Fletcher–Reeves CG | Polak–Ribière⁺ (implemented) — same cost, strictly better restart behaviour on non-quadratic surfaces. |
| Newton–Raphson | Trust-region / RFO — raw Newton steps toward *any* stationary point, including saddles, and diverge far from the minimum. The regularized variants are the ones worth having. |
| Quick-min, damped dynamics | FIRE — the same damped-MD idea with the adaptive timestep that makes it competitive. |
| Nudged elastic band | Climbing-image NEB — plain NEB does not converge on the saddle itself, which is the point of running it. |
| String method | Growing / freezing string — same path, without needing a good initial guess for the whole chain. |

**Steepest descent stays ○ on CPU** rather than ⊘: it is dominated by CG for
minimization proper, but it is the robust first-resort on a badly broken start
geometry, and it is already implemented on CUDA.
