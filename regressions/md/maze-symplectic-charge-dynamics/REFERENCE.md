# maze-symplectic-charge-dynamics

oracle: metalwalls
oracle_version: literature — Coretti et al., J. Chem. Phys. 152, 194701 (2020), DOI 10.1063/5.0007192 (no MetalWalls binary installed; perf_exempt)

## Physical statement (the literature pin)

Coretti, Scalfi, Bacon, Rotenberg, Vuilleumier, Salanne, Ciccotti, Bonella,
*"Mass-zero constrained molecular dynamics for electrode charges in
simulations of electrochemical systems"*, **J. Chem. Phys. 152, 194701
(2020)**, establishes for constant-potential electrode charges evolved as
mass-zero constrained degrees of freedom (SHAKE-like solve of the
constant-potential constraint, Verlet-extrapolated initial guess, corrections
along the constraint gradients):

1. the constraint \(\sigma_i = \partial U/\partial Q_i - \Psi_i\) is
   satisfied at solver precision **every step** (their §II.C: for the linear
   constant-potential constraint SHAKE reduces to a linear solve driven to
   machine precision);
2. the resulting ionic dynamics is **time-reversible** — unlike dissipative
   predictor / loose-SCF adiabatic propagation (their §II.B, §III);
3. NVE total energy (the constant-potential grand energy
   \(E = KE + U - \sum_i \Psi_i Q_i\); MetalWalls JOSS **5**, 2373 (2020)
   eq. 6) shows **no secular drift**, while loose-tolerance adiabatic
   propagation drifts (their Fig. 2 discussion).

MetalWalls implements this MaZe path as one of its charge solvers
(Marin-Laflèche et al., JOSS **5**, 2373 (2020), DOI 10.21105/joss.02373).
No MetalWalls binary is available offline in this environment, so the pins
are the algorithmic invariants above (exact statements, not fitted numbers)
plus internal twin discrimination — see `tools/metalwalls_oracle/
pin_notes_maze-symplectic-charge-dynamics.md`.

## Geometry

Parallel-plate capacitor, orthorhombic box \(6 \times 6 \times 20\) Å,
min-image electrostatics (`PairGaussianCoul`, cutoff 12 Å > max min-image
distance 10.9 Å — every pair interacts exactly once, identical kernels to
the elastance build):

| Region | type | N | z [Å] | η [1/Å] | Ψ [V] | mass [amu] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Bottom plate | 1 | 4 | 4 | 2.0 | 0.0 | static |
| Top plate | 2 | 4 | 16 | 2.0 | +1.0 | static |
| Ion + | 3 | 1 | 8.5 | point | — | 40.0 |
| Ion − | 3 | 1 | 11.5 | point | — | 40.0 |

Ions ±0.5 e at (2.1, 2.1, 8.5) and (4.2, 4.2, 11.5) Å — lateral min-image
separation 2.97 Å, offsets chosen away from the L/2 min-image force kink —
zero initial velocity. Explicit NVE velocity Verlet (public Component API;
the driver NeighborList cannot cover the no-cutoff min-image elastance
twin), dt = 0.25 fs, 400 steps forward, velocity flip, 400 steps back
(palindrome).

## Hard-coded pins (asserted in the regression source)

| Pin | Value | Source |
| --- | --- | --- |
| Constraint residual \(\lVert\sigma\rVert_\infty\), every step | ≤ 1e-10 V | Coretti 2020 §II.C (solver-precision constraint) |
| Electrode electroneutrality \(|\sum Q|\), every step | ≤ 1e-9 e | constraint of the KKT system |
| Palindrome ion-position return error (MaZe) | ≤ 1e-6 Å | Coretti 2020 time-reversibility |
| Palindrome electrode-charge return error (MaZe) | ≤ 1e-6 e | Coretti 2020 time-reversibility |
| MaZe return error < 0.1 × loose CG+ASPC twin | ratio | reversible vs dissipative propagation |
| NVE grand-energy drift (MaZe, after 5-step settle) | ≤ 1e-6 eV | Coretti 2020 no-secular-drift |
| MaZe drift < loose (cg_tol 1e-3) CG+ASPC twin drift | ratio | Coretti 2020 vs adiabatic loose solve |
| \(Q_{top} > 0\) at ΔΨ = +1 V | polarity | electrostatics sanity |
| Ion min separation > 1.8 Å | validity | bound smooth trajectory, no collapse |

Units: eV, Å, fs, amu, V, e (Atomiverse public surface).

## Measured at pin time (2026-07-27, Release, GCC 14.3, x86_64 epyc)

| Metric | MaZe symplectic | loose (1e-3) CG+ASPC twin |
| --- | ---: | ---: |
| Grand-energy drift [eV] | 7.57e-08 | 1.76e-06 |
| Palindrome position return [Å] | 1.4e-14 | 1.07e-05 |
| Palindrome charge return [e] | 9.0e-13 | 5.09e-05 |
| Max constraint residual [V] | 9.99e-11 | — |
| Max \(|\sum Q|\) [e] | 5.6e-17 | 5.6e-17 |
| \(Q_{top}\) first solve [e] | +0.135822 | +0.135821 |
| min ion separation [Å] | 3.798 | 3.798 |

Pins above carry 10×+ margins over these measurements.

## Runtime

Two 800-step NVE runs of a 10-atom cell (8×8 elastance): < 1 s (≪ 60 s).
