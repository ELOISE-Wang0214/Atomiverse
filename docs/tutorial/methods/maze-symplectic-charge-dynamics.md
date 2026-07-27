# True mass-zero constrained (MaZe) electrode charge dynamics

*(Technical-report draft; narrative house style TBD — rewrite later.)*

| Field | Value |
| ----- | ----- |
| Slug | `maze-symplectic-charge-dynamics` |
| Regression | `regressions/md/maze-symplectic-charge-dynamics/` |
| Public API | `cpu::md::ConstPotMaZe` members `symplectic` (default `true`), `maze_tol` [V], `max_maze_iters`; accessors `constraint_residual()`, `last_maze_iters()`, `max_dQ()` (Component on `cpu::md::MD`, default `pre_force`) |
| Roadmap / queue | `docs/developer/metalwalls-method-queue.md` wave 5; map §1 MaZe ~ → ✓ |
| Oracle | **MetalWalls** physics — literature pin (Coretti et al. 2020; no MW binary → `perf_exempt`) |

## Summary

Wave 4 landed `ConstPotMaZe` as an **adiabatic foothold**: CG + Kolafa ASPC
on `ConstPotBO` each `pre_force`. This wave replaces its default path with
the **true mass-zero constrained (SHAKE-like) charge dynamics** of Coretti
et al., *J. Chem. Phys.* **152**, 194701 (2020): electrode charges are
auxiliary degrees of freedom whose Verlet extrapolation
\(\tilde Q = 2Q_n - Q_{n-1}\) provides the (time-reversible) initial state,
and SHAKE corrections **along the constraint gradients**
(\(Q_{n+1} = \tilde Q + A\gamma\), CG on \(A^2\gamma = -\sigma(\tilde Q)\))
drive the constant-potential constraint
\(\sigma_i = (AQ - (\Psi - b))_i - \chi\) to `maze_tol` every step
(\(\chi\) = uniform potential shift of the simultaneous electroneutrality
constraint, applied along the cached \(u = A^{-1}\mathbf 1\) direction).
The CG residual of the SHAKE system **is** the constraint violation, so the
stopping rule bounds \(\lVert\sigma\rVert_\infty\) directly.

`symplectic = false` keeps the wave-4 CG+ASPC adiabatic twin (bitwise the
old behaviour) for comparison studies. Still a **Component** on
`cpu::md::MD` — no new Driver category, no charge solve inside `KSpace*`.

## Why it matters (the paper's result, reproduced by the regression)

- Solved tightly, the electrode charges are a function of the ionic
  positions only — the NVE ionic dynamics is **time-reversible** and the
  constant-potential grand energy \(KE + U - \sum_i\Psi_iQ_i\) shows **no
  secular drift**.
- A loose-tolerance CG+ASPC adiabatic twin (production-style
  `cg_tol = 1e-3`) is dissipative: the velocity-flip palindrome does not
  retrace and the grand energy drifts. The regression pins both contrasts
  on an identical 10-atom parallel-plate NVE trajectory.

## Composition

```cpp
cpu::md::ConstPotMaZe maze(/*eta [1/Å]=*/2.0);
maze.add_electrode(/*type=*/1, /*Psi [V]=*/0.0);
maze.add_electrode(/*type=*/2, /*Psi [V]=*/1.0);
maze.maze_tol = 1e-10;          // ‖σ‖∞ constraint tolerance [V]

cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{0.25});
sim.add_potential(cpu::PairGaussianCoul{6.0, {0.0, 2.0, 2.0, 0.0}});
sim.set_neighbor(cpu::md::NeighborList<>{6.0, 0.5}); // cutoff < L/2
sim.add_ref(maze);              // pre_force −50 default binding
sim.add(cpu::md::FixFreeze{1});
sim.add(cpu::md::FixFreeze{2});
sim.build();
sim.run(frame, 400);
maze.constraint_residual();     // ≤ maze_tol [V]
maze.last_maze_iters();         // 0 when the Verlet guess already satisfies σ
```

The regression itself drives the components with an explicit velocity-Verlet
loop instead of the `MD` driver: its conservation twin needs the pair energy
to see **every** min-image pair (the elastance/RHS build has no cutoff),
which a production NeighborList (cutoff < L/2) cannot cover in a small cell.

## Tests

- Unit: `tests/cpu/dynamics/test_const_pot_maze_symplectic.cpp` — matches
  dense BO per site, neutrality (electrode-only and global), zero-iteration
  re-solve on a static frame, Verlet warm start beats cold start, foothold
  twin unchanged, parameter validation.
- E2E: `reg_maze-symplectic-charge-dynamics` — constraint ≤ 1e-10 V every
  step, palindrome reversibility, grand-energy conservation vs the loose
  CG+ASPC twin. Pins in `REFERENCE.md`.
