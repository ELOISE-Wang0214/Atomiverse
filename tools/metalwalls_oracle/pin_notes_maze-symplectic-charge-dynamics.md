# Pin notes — maze-symplectic-charge-dynamics

- Wave: 5 (metalwalls-method-loop), 2026-07-27
- Oracle priority: **metalwalls** (MaZe constant-potential charge dynamics is
  MetalWalls-only physics — no LAMMPS twin: LAMMPS ELECTRODE has conp/conq
  matrix/CG solvers but no mass-zero constrained SHAKE charge dynamics).
- No MetalWalls binary is installed on this host → **literature pin** with
  full citation, `perf_exempt` per the binding oracle-priority rule.

## Primary citation

Coretti, Scalfi, Bacon, Rotenberg, Vuilleumier, Salanne, Ciccotti, Bonella,
*Mass-zero constrained molecular dynamics for electrode charges in
simulations of electrochemical systems*, **J. Chem. Phys. 152, 194701
(2020)**, DOI 10.1063/5.0007192.

Supporting: Marin-Laflèche et al., *MetalWalls: A classical molecular
dynamics software dedicated to the simulation of electrochemical systems*,
JOSS **5**(53), 2373 (2020), DOI 10.21105/joss.02373 (MaZe is one of the MW
charge solvers; grand-energy bookkeeping eq. 6).

## What is pinned (algorithmic invariants, not fitted numbers)

The paper's testable statements for the linear constant-potential constraint:

1. **Constraint at solver precision every step.** SHAKE for the linear
   constraint reduces to an iterative linear solve driven to (near) machine
   precision; pinned as \(\lVert AQ-(\Psi-b)-\chi\mathbf 1\rVert_\infty
   \le 10^{-10}\) V each step.
2. **Time reversibility.** With the charge solution independent of history
   (tight solve; Verlet-extrapolated guess only affects iteration count),
   an NVE velocity-flip palindrome retraces the ionic trajectory. Pinned:
   return error ≤ 1e-6 Å and < 0.1 × the loose (cg_tol 1e-3) CG+ASPC
   adiabatic twin's return error on the identical trajectory.
3. **No secular energy drift.** Grand energy \(KE+U-\sum\Psi_iQ_i\)
   conserved; pinned ≤ 1e-6 eV over 2×400 steps and below the loose twin.

The exactness of the converged charges themselves (MaZe = dense BO solution,
same KKT system) is pinned in the **unit** tests
(`tests/cpu/dynamics/test_const_pot_maze_symplectic.cpp`) against the dense
Cholesky `ConstPotBO` path, which earlier waves pinned against the
MetalWalls/LAMMPS-ELECTRODE plate-capacitor algebra
(`pin_notes.md`, `reg_const-pot-bo-parallel-plate`).

## Measured values at pin time (2026-07-27, this host, Release build)

Identical 10-atom parallel-plate NVE palindrome (2 × 400 steps, dt 0.25 fs):

| Metric | MaZe symplectic | loose (cg_tol 1e-3) CG+ASPC twin |
| --- | ---: | ---: |
| Grand-energy drift [eV] | 7.57e-08 | 1.76e-06 |
| Palindrome position return [Å] | 1.4e-14 | 1.07e-05 |
| Palindrome charge return [e] | 9.0e-13 | 5.09e-05 |
| Max constraint residual [V] | 9.99e-11 | — |

The MaZe path is ~9 orders of magnitude more reversible and 23× better
conserving than the loose adiabatic twin — the Coretti 2020 contrast.
Fair wall-clock vs MetalWalls not possible (no binary) →
`perf_exempt=true`.
