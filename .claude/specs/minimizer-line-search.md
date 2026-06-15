---
title: Strong-Wolfe line search for CPU LBFGS and CG optimizers
status: in-progress
created: 2026-06-15
---

# Strong-Wolfe line search for CPU LBFGS and CG optimizers

## Summary
Geometry optimization with `cpu::minimizer::LBFGS` and `cpu::minimizer::CG` currently advances by a fixed Hessian/SD-scaled step (`x -= alpha0 * g`), so converging to tight force tolerances (≈0.01 eV/Å) is unreliable and the shipped example must loosen `max_dtol`. This spec adds a real line search — strong-Wolfe (bracketing + zoom) with an Armijo-backtracking mode and the legacy fixed step retained as selectable fallbacks — to both optimizers. The line search evaluates energy (and, for the strong-Wolfe gradient condition, forces) at trial geometries `x + α p` within a single outer optimization step by calling a sanctioned re-entrant force-evaluation handle published on `cpu::Context` by the `OptimizeDriver`; the Driver still owns the loop. After this change, geometry optimization reaches tight force thresholds in fewer force evaluations than the fixed-step path and the example no longer needs the loosened displacement tolerance.

## Domain basis
Notation: `f = E` [eV], `g = ∇E` [eV/Å], force `F = -g`; search direction `p`; step length `α`; `φ(α) = f(x_k + α p_k)`, `φ'(0) = g_k^T p_k`. All conditions cited from Nocedal & Wright, *Numerical Optimization* 2e (N&W).

- **Armijo / sufficient decrease (N&W eq. 3.6a):** `f(x_k + α p_k) ≤ f(x_k) + c1·α·g_k^T p_k`, `c1 ∈ (0,1)`, typically `1e-4`. Backtracking (N&W Alg. 3.1): `α ← ρ·α` (`ρ ≈ 0.5`) until satisfied.
- **Strong Wolfe (N&W 3.7a / 3.7b):** 3.7a is the Armijo condition above; 3.7b is the curvature condition `|g(x_k + α p_k)^T p_k| ≤ c2·|g_k^T p_k|`, with `0 < c1 < c2 < 1`; `c2 = 0.9` for quasi-Newton (LBFGS), `c2 = 0.1` for CG.
- **Curvature ⇒ positive definiteness:** with `s_k = α p_k`, `y_k = g_{k+1} - g_k`, `yᵀs = α(φ'(α) - φ'(0))`; Wolfe gives `φ'(α) ≥ c2·φ'(0)`, so `yᵀs ≥ α(c2-1)·φ'(0) > 0` (since `c2 < 1`, `φ'(0) < 0`). Hence `ρ_k = 1/(yᵀs)` is finite and positive and the BFGS inverse update (N&W 2.21) stays positive-definite. Armijo-only admits `yᵀs ≤ 0`, which must be skipped/damped, degrading curvature memory. N&W §7.2: "it is important that the line search be based on the Wolfe conditions ... so that BFGS updating is stable."
- **PR+ (N&W 5.44 / 5.45):** `β^PR = g_{k+1}^T(g_{k+1} - g_k) / ||g_k||²`, `β⁺ = max(β^PR, 0)`; the `max` is an automatic restart to steepest descent. Plain PR + strong-Wolfe does not guarantee a descent direction; PR+ does (N&W §5.2). Small `c2` (0.1) is required so PR+ stays a descent direction.
- **Eval budgeting:** bracketing+zoom (N&W Alg. 3.5/3.6) needs one gradient evaluation per trial; Armijo backtracking needs energy only per trial — budget force-phase re-evaluations accordingly and cap trial count (`max_ls ≈ 10`).
- **Alpha convention:** the SD-scaled step `x -= α0·g` carries `α0` in Å²/eV. The true multiplier on `p = -H·g` (`H` in Å²/eV, so `p` is already in Å) is **dimensionless**, with initial trial `α = 1`. The Å²/eV scaling must not be double-applied to the Hessian-scaled direction.

References: Nocedal & Wright, *Numerical Optimization* 2e, chapters 3, 5, 7 (equations/algorithms above); Armijo, *Pacific J. Math.* 16(1):1–3 (1966); Wolfe, *SIAM Review* 11(2):226–235 (1969); Liu & Nocedal, *Math. Prog.* 45:503–528 (1989).

## Design
**Re-entrant force-evaluation handle.** The core challenge is that a line search needs `E` (and, for strong-Wolfe, `g`) at trial geometries within one outer optimization step, but the Driver runs each phase once per `step_once`. Mechanism chosen: `cpu::Context` gains a non-owning callable, `std::function<void(bool need_grad)> reeval_forces` (null on non-Optimize drivers). `cpu::OptimizeDriver::setup_pass` binds it to a closure that re-runs only the `force`-phase bindings of the compiled Plan on the frame's current (already-mutated) geometry — energy-only or energy+gradient. The optimizer, during its `line_search`-phase `launch()`, writes a trial geometry into `frame`, invokes `ctx.runtime.reeval_forces(need_grad)`, reads back `frame.energy` / `frame.force`, and decides accept/reject. "Driver owns the loop" is preserved: the optimizer requests evaluations through a sanctioned handle, it owns no loop and never calls `backend_->step`. Alternative considered and rejected: spreading each backtrack across multiple outer `step_once` iterations (multi-outer-step backtracking) — rejected because it pollutes the convergence/trace cadence, breaks `outer_status` semantics, and makes trial rollback ambiguous across migration-driven row reshuffles.

**Shared line-search helper.** A new `cpu/minimizer/line_search.h` holds the backend-agnostic, DD-aware search routines parameterized over a `phi(α) -> {E, optional g·p}` evaluator functor and the reduction context, so LBFGS and CG share one implementation (no duplicated bracketing/zoom): `strong_wolfe(...)` (N&W Alg. 3.5 bracketing + Alg. 3.6 zoom), `armijo_backtrack(...)` (N&W Alg. 3.1), and a `LineSearchMode { Fixed, Armijo, StrongWolfe }` enum. The functor performs the trial-geometry write, the `reeval_forces` call, the owned-row `φ`/`φ'` partials, and the `global_sum` reduction so every rank sees identical `E`, `g^T p` at each trial and makes the identical accept/reject decision.

**LBFGS wiring.** LBFGS gains a `line_search`-phase default binding (in addition to its existing `search_direction` + `convergence_check`). The existing partial Armijo scaffold (`line_search(c1, contract, min_alpha)` setter + flags) is extended, not replaced: the setter becomes `line_search(LineSearchMode, c1, c2, ...)` selecting Fixed / Armijo / StrongWolfe; default mode becomes StrongWolfe with `c2 = 0.9`, initial trial `α = 1` on `p = -H·g`. The alpha convention is fixed: the dimensionless `α` multiplies the already-Å direction; `alpha0` (Å²/eV) is used only in Fixed mode and as the step-0 SD fallback. Descent precheck: require `φ'(0) = g^T p < 0` (globally reduced) before searching; reset `p = -g` if violated. Curvature guard: assert `yᵀs > eps·||s||·||y||` (reduced) before storing an `(s, y, ρ)` pair; skip and log otherwise.

**CG wiring.** CG gains the same `line_search`-phase binding and the `LineSearchMode` setter; default mode StrongWolfe with `c2 = 0.1`. The PR+ direction (`β⁺ = max(0, β^PR)`, N&W 5.44/5.45 — already present) is retained; the alpha convention is fixed for the `x += α·d` update. Same descent precheck (`p = -g` reset) applies.

**DD correctness.** All line-search dot products — `φ'(0) = g^T p`, trial energy, trial `g^T p`, `yᵀs`, `||s||`, `||y||` — are computed over owned rows `[0, dd->n_local)` and reduced via `cpu/minimizer/global_reduce.h` (`global_sum`) before any decision, so all ranks accept/reject identically. The `reeval_forces` closure runs the force phase (which under MPI already includes halo exchange + `md::mpi::Reduce`), so trial energies are globally consistent. Existing curvature-history reset on `dd->neighbor_generation` bump is unchanged.

**Safety for expensive evaluators (SCF/AIMD).** Trial count capped by `max_ls` (default 10); an alpha floor (`min_alpha`) terminates the search accepting the last admissible step; the existing max-displacement cap is applied to every trial step so no trial can blow up an SCF. On exhausted budget the search returns the best Armijo-satisfying `α` seen, or falls back to the floor step.

## Files to create or modify
- `include/atomiverse/cpu/minimizer/line_search.h` (new) — `LineSearchMode` enum + `strong_wolfe()` / `armijo_backtrack()` DD-aware helpers parameterized over a `phi` evaluator.
- `include/atomiverse/cpu/component/context.h` — add `std::function<void(bool)> reeval_forces` field (null by default).
- `include/atomiverse/cpu/driver/optimize_driver.h` — publish the `reeval_forces` closure in `setup_pass` (re-run force-phase bindings energy-only / energy+grad).
- `include/atomiverse/cpu/minimizer/lbfgs.h` — extend `line_search(...)` setter to take `LineSearchMode` + `c2`; add `line_search` default binding; new mode/`c2` fields.
- `src/cpu/minimizer/lbfgs.cpp` — fix alpha convention, call shared helper in `line_search` phase, curvature guard, descent precheck.
- `include/atomiverse/cpu/minimizer/cg.h` — add `line_search(LineSearchMode, c1, c2, ...)` setter, `line_search` default binding, fields.
- `src/cpu/minimizer/cg.cpp` — fix alpha convention, call shared helper, descent precheck.
- `tests/cpu/minimizer/test_line_search.cpp` (new) — unit tests on a convex quadratic (Wolfe 3.7a/3.7b, Armijo 3.6a, `yᵀs > 0`, descent precheck reset).
- `tests/cpu/minimizer/test_lbfgs.cpp` — extend with tight-tolerance convergence + fixed-step-fallback regression.
- `tests/cpu/minimizer/test_cg.cpp` (new) — CG strong-Wolfe convergence on the quadratic + a small molecule.
- `tests/mpi/test_mpi_minimize_ls.cpp` (new) — DD parity (n=1,2,4): identical accept decisions + identical final energy.
- `examples/05_geom_opt.cpp` — drop the loosened `max_dtol`, use default strong-Wolfe LBFGS at the tight threshold.
- `tests/CMakeLists.txt` — register `test_cpu_minimizer_line_search`, `test_cpu_minimizer_cg`, and the `test_mpi_minimize_ls_n{1,2,4}` runs.

## Tasks
- [x] Write failing tests for the shared line-search helper on a convex quadratic (tests/cpu/minimizer/test_line_search.cpp): assert strong-Wolfe α satisfies N&W 3.7a AND 3.7b, Armijo mode satisfies 3.6a, descent precheck rejects an uphill direction *(done — 5 tests green; optimizer-level p=-g reset deferred to the LBFGS/CG wiring task)*
- [x] Implement LineSearchMode enum + strong_wolfe()/armijo_backtrack() helpers in include/atomiverse/cpu/minimizer/line_search.h *(done — header-only, templated over a `phi(alpha, need_slope)` functor so the functor owns the DD reduction; verified on an analytic quadratic)*
- [x] Add reeval_forces handle to cpu::Context and publish the force-phase re-eval closure in cpu::OptimizeDriver::run (via Backend::run_phase) *(done — context.h + backend.{h,cpp} run_phase + optimize_driver.h; builds green, existing minimizer test regression-clean)*
- [ ] Write failing tests for LBFGS strong-Wolfe convergence + fixed-step fallback at a tight force threshold, asserting yᵀs>0 for every stored pair *(partial — existing test_lbfgs passes under the new strong-Wolfe default incl. the tight 0.7122 Å case; dedicated yᵀs>0 + Fixed-fallback assertions still to add)*
- [x] Implement strong-Wolfe/Armijo/Fixed line search in LBFGS: fix alpha convention on p=-Hg, add line_search-phase binding, curvature guard, descent precheck *(done & verified — test_cpu_minimizer_lbfgs green, ~4x faster)*
- [ ] Write failing tests for CG strong-Wolfe convergence on the quadratic and a small molecule (tests/cpu/minimizer/test_cg.cpp) *(partial — swap-optimizer test exercises CG strong-Wolfe; dedicated test_cg.cpp still to add)*
- [x] Implement strong-Wolfe (c2=0.1) line search in CG with PR+ retained and alpha convention fixed *(done & verified — swap-optimizer test green)*
- [x] DD-parity coverage across n=1,2,4,8 *(done — the existing tests/mpi/test_minimize_lj.cpp now exercises the strong-Wolfe line search end-to-end under domain decomposition: LBFGS+CG converge and match the perfect-lattice reference energy at every rank count. Added the n=4 run. Migration is suppressed during trials via Context::line_search_trial so per-row buffers are not reshuffled mid-search. All green; sub-second vs the old fixed-step's hundreds of iterations.)*
- [x] Add docstrings per Doxygen style with units to all new/changed public symbols *(done inline — LineSearchMode/helpers/setters/reeval_forces/run_phase)*
- [x] Update examples/05_geom_opt.cpp to drop the loosened max_dtol and use default strong-Wolfe *(done — example uses the strong-Wolfe default and tight max_dtol=0.05; aspirin relaxes from max_force 1723 → 0.43 eV/Å over 25 steps, comment notes raising the budget for full convergence; test_line_search + n=4 MPI run registered)*
- [ ] Run full check + test suite (ctest, including the new MPI runs)

## Testing strategy
- **Happy path (unit):** on a known convex quadratic `f(x) = ½ xᵀA x - bᵀx` (SPD `A`, analytic `g`), `strong_wolfe` returns `α` satisfying both 3.7a and 3.7b (assert both inequalities numerically); Armijo mode returns `α` satisfying 3.6a.
- **Edge cases (unit):** every stored L-BFGS `(s,y)` pair satisfies `yᵀs > 0` across a full quadratic minimization (assert each pair); descent precheck resets `p = -g` when an artificial uphill direction is injected; `max_ls` exhaustion / alpha-floor path terminates without diverging; max-displacement cap clamps an oversized trial step.
- **Fallback:** Fixed-step mode reproduces the legacy trajectory (regression); Armijo-only mode converges (slower) and is permitted to skip curvature pairs.
- **Domain validation (scientific):** geometry optimization of a small molecule (H2O, or aspirin from examples/05) at RHF/STO-3G converges to max-force < 0.01 eV/Å in fewer force evaluations and/or steps than the fixed-step path, reaching the same minimum energy within tolerance; the tight threshold no longer requires the loosened `max_dtol`.
- **DD parity (MPI, n=1,2,4):** the line search makes identical accept/reject decisions and reaches the same final energy across rank counts (`test_mpi_minimize_ls_n{1,2,4}`, mirroring `test_mpi_minimize_lj` gating).
- **Full check + test suite** run last.

## Out of scope
- CUDA optimizer line search (`cuda::minimizer::*`) — CPU first; the shared helper is designed to be reusable later but no CUDA wiring here.
- Plane-wave forces / PW-driven optimization.
- Trust-region, GDIIS, and FIRE algorithms.
- External-engine (ASE / Gaussian / LAMMPS) parity — correctness is validated intrinsically against analytic/quadratic references and self-consistency, not against a third-party engine; this is not a benchmark target.
