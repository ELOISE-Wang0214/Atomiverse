# maze-symplectic-charge-dynamics

- status: code-complete (all 9 acceptance criteria verified 2026-07-27; dual
  gate green: 4/4 ctest incl. pre-existing const-pot suites,
  run_reg_maze-symplectic-charge-dynamics PASS, reg_maze-charge-dynamics
  foothold twin PASS under the new default)
- loop: metalwalls-method-loop (wave 5)
- oracle: metalwalls (no binary installed → literature pin, Coretti et al. 2020; perf_exempt)
- map: `.claude/notes/metalwalls-capability-map.md` §1 "MaZe mass-zero charge dynamics"

## Problem

`ConstPotMaZe` (wave 4) is an adiabatic foothold: it is a thin owner of
`ConstPotBO` running CG + Kolafa ASPC each `pre_force`. The ASPC predictor is
dissipative and a loose CG tolerance produces a systematic energy drift in
NVE constant-potential MD. The true MaZe algorithm (Coretti, Scalfi, Bacon,
Rotenberg, Vuilleumier, Salanne, Ciccotti, Bonella, *J. Chem. Phys.* **152**,
194701 (2020)) evolves electrode charges as mass-zero constrained degrees of
freedom: Verlet extrapolation of the charges provides the (time-reversible)
initial guess and a SHAKE-like iteration applies corrections **along the
constraint gradients** until the constant-potential constraint
\(\sigma_i = (A Q - (\Psi - b))_i - \chi = 0\) (with \(\chi\) the uniform
potential shift enforcing electroneutrality) is satisfied to machine
precision. For the linear constant-potential constraint the SHAKE system is
\(A^2 \gamma = -\sigma(\tilde Q)\), \(Q_{n+1} = \tilde Q + A\gamma\) with
\(\tilde Q = 2Q_n - Q_{n-1}\). Solved tightly, the resulting ionic dynamics
is time-reversible and shows no secular energy drift — the paper's headline
result vs dissipative predictor/loose-SCF schemes.

## Design (arch law)

Upgrade `ConstPotMaZe` in place (reuse, never duplicate — Component on
`cpu::md::MD` at `pre_force`, never a Driver category, never inside KSpace*):

- New public members (ctor args/members — **no** Config bag):
  - `bool symplectic = true` — true mass-zero constrained SHAKE path
    (Coretti et al.). `false` keeps the wave-4 CG+ASPC adiabatic foothold
    twin (unchanged behaviour, used as the comparison twin).
  - `double maze_tol = 1e-10` — constraint residual tolerance
    \(\lVert\sigma\rVert_\infty\) [V].
  - `int max_maze_iters = 2000` — cap on constraint-CG iterations.
- New accessors: `constraint_residual()` [V], `last_maze_iters()`,
  existing `max_dQ()` [e] kept.
- Electroneutrality: neutrality constraint gradient is the ones vector;
  solved simultaneously with the potential constraints via the cached
  \(u = A^{-1}\mathbf{1}\) direction (uniform potential shift \(\chi\)),
  identical KKT solution to `ConstPotBO`'s projection.
- Implementation reuses `ConstPotBO` internals (matrix/RHS build, charge
  write-back, electrode forces) via friendship — no duplicated elastance
  code, no new public surface on `ConstPotBO`.
- Bootstrap: first launch does a cold tight solve (no history); Verlet
  extrapolation from the second launch on.

## Acceptance criteria

1. **Exactness**: on a static parallel-plate frame, `symplectic` charges
   match the dense-Cholesky `ConstPotBO` solution to ≤ 1e-8 e per site and
   `constraint_residual() ≤ maze_tol` after every launch. [unit]
2. **Electroneutrality**: `neutrality_residual() ≤ 1e-10` e every step,
   including global electroneutrality mode. [unit]
3. **Verlet warm start**: on a smoothly moving electrolyte trajectory the
   mean constraint-CG iteration count with history is strictly lower than
   cold-start (history cleared) solves of the same sequence. [unit]
4. **History correctness**: after two launches on a static frame the Verlet
   guess equals the converged charges (zero-iteration re-solve;
   `last_maze_iters() == 0`). [unit]
5. **Time-reversibility (regression)**: NVE constant-potential MD
   (parallel plate + `PairGaussianCoul` + `FixFreeze` walls +
   `VelocityVerlet`), forward `n` steps, velocity flip, backward `n` steps:
   tight symplectic MaZe returns ion positions and electrode charges to
   their initial values orders of magnitude closer than the loose
   (`cg_tol = 1e-3`) CG+ASPC foothold twin; hard-coded pins with margins.
6. **Energy conservation (regression)**: over the same NVE trajectory the
   conserved quantity \(E = KE + U_{\mathrm{pair}} - \sum_i \Psi_i Q_i\)
   drifts less with symplectic MaZe than with the loose CG+ASPC twin, and
   below an absolute hard-coded pin.
7. **Regression packaging**: `regressions/md/maze-symplectic-charge-dynamics/`
   out-of-library consumer of public headers only, ≤ 60 s, targets
   `reg_maze-symplectic-charge-dynamics` + `run_reg_maze-symplectic-charge-dynamics`
   registered; `REFERENCE.md` states `oracle: metalwalls`,
   `oracle_version`, full Coretti 2020 citation, units.
8. **MD-only build**: everything builds and passes with the binding MD-only
   configure (HF/PW/GPW/CUDA/TORCH all OFF).
9. Foothold twin unchanged: existing `reg_maze-charge-dynamics` still passes.

## Out of scope

- Nonlinear constraints (PIM dipole MaZe), MPI-distributed constraint solve,
  CUDA port, MetalWalls-binary parity (no binary installed — literature pin).
