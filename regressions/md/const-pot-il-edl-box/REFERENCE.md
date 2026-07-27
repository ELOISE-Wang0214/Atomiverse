# const-pot-il-edl-box

oracle: metalwalls
oracle_version: MW-form BO algebra, offline numpy pin (no `mw` binary on host);
  scipy 1.18.0 / numpy 2.5.1, x86_64 login node, 2026-07-27
notes: Wave-5 production EDL box — planar graphene-like constant-potential
  walls + non-polarizable single-site coarse-grained IL, short NVT.

## System

- Box 17.04 × 17.216585027234639 × 40 Å (orthorhombic, 3D PBC).
- Two honeycomb C walls (a_CC = 1.42 Å, 4-atom rectangular cell, 4×7 cells →
  112 sites/wall): z = 12 Å (type 2, Ψ = −1 V), z = 28 Å (type 3, Ψ = +1 V),
  frozen (`FixFreeze`). Interior (IL) gap 16 Å < periodic return path 24 Å,
  so the EDL response dominates.
- 12 + 12 single-site CG IL ions (±1 e, non-polarizable, unscaled charges) on
  a deterministic 3×3(−centre)×3 grid, checkerboard parity plus one charge
  swap per outer layer (5+/3− facing the −1 V wall, 3+/5− facing the +1 V
  wall — a nascent-EDL polarization so the step-0 solve shows screening).
- LJ: cation σ = 4.38 Å, ε = 2.56 kJ/mol = 0.026532 eV, m = 139.22 amu
  (BMIM+-like); anion σ = 5.06 Å, ε = 4.71 kJ/mol = 0.048815 eV,
  m = 144.96 amu (PF6−-like); wall C σ = 3.37 Å, ε = 0.23 kJ/mol =
  0.0023837 eV. Ion–ion cross terms by Lorentz–Berthelot; ion–wall terms as
  explicit LB entries; the wall types carry **no LJ diagonal**, so the
  frozen honeycomb (C–C 1.42 Å ≪ σ_CC) has zero intra-wall LJ (the LAMMPS
  `neigh_modify exclude` analogue). Sizes/masses follow the coarse-grained
  imidazolium family of Roy & Maroncelli, J. Phys. Chem. B 114, 12629 (2010)
  and Merlet et al., J. Phys. Chem. C 116, 7687 (2012) (each ion reduced to
  its dominant interaction site; this regression pins solver algebra and
  physics invariants, not FF-sensitive observables).
- Electrostatics: `CoulLong{α=0.3 Å⁻¹, r_c=8 Å}` + `KSpacePME{α=0.3,
  spacing 0.5 Å, order 5}`; LJ cutoff 10 Å; NeighborList 10 Å + 0.5 skin.
- Const-pot: `cpu::md::ConstPotBO`, η = 1.805132 Å⁻¹ (MetalWalls graphite
  Gaussian width 0.955234 bohr⁻¹; MW example inputs), CG solve
  (`use_cg`, tol 1e-10) + ASPC order 2, electrode-only electroneutrality,
  `xy_image_range = 0` (min-image E1 foothold, documented).
- NVT: VelocityVerlet dt = 1 fs, Nose–Hoover T = 400 K, τ = 100 fs,
  600 steps, `InitVelocity{400 K, seed 11}`.

## Oracle

MetalWalls-form constant-potential Born–Oppenheimer solve (Scalfi et al.,
J. Chem. Phys. 153, 174704 (2020); Marin-Laflèche et al., JOSS 5, 2373 (2020),
10.21105/joss.02373). No MetalWalls binary is installed on this host and
LAMMPS ELECTRODE is not built into the available `lmp`, so the pin is an
**independent offline numpy implementation of the same algebra** on the exact
deterministic initial geometry:
`tools/metalwalls_oracle/inputs/metalwalls/const-pot-il-edl-box/pin_il_edl.py`
(same kernels as `ConstPotBO` with `xy_image_range=0`: Gaussian self
k_e·η·√(2/π), Gaussian–Gaussian erf((η/√2)r)/r, Gaussian–point erf(η·r)/r,
min-image, electrode-only neutrality projection; direct `numpy.linalg.solve`
vs ATV CG at tol 1e-12).

## Pins (hard-coded in the regression, tol 1e-6 e)

| solve | Q_left [e] | Q_right [e] |
| ----- | ---------- | ----------- |
| empty capacitor (IL charges zeroed) | −4.965149459545e−01 | +4.965149459545e−01 |
| with IL (nascent-EDL config) | −1.552770528510e+00 | +1.552770528510e+00 |

Derived assertions: EDL screening |Q|_IL / |Q|_empty ≈ 3.13 > 1; polarity
Q_left < 0 < Q_right; electrode neutrality ≤ 1e-8 e.

## Run-time invariants (trajectory-dependent, therefore not hard-pinned)

- finite E, K (> 0), charges each 100-step chunk; |Q_L + Q_R| ≤ 1e-6 e;
- polarity preserved after 600 fs; walls immobile to 1e-10 Å;
- Ψ(z) (PotentialProfile, 80 bins over [0, 40 Å)): Ψ(z_R) − Ψ(z_L) > 0
  (positive applied ΔΨ = 2 V);
- ρ_s(z) (SpeciesDensityProfile, 32 bins over [12, 28 Å)): all 12+12 ions
  confined between the walls at the final step; run-mean profile
  normalization Σ ρ_s·A·dz = 12 ± 1e-6 per species; first-layer (4 Å)
  cation/anion occupancies printed as EDL layering diagnostics.

perf: **perf_exempt** — MW-only physics (Gaussian constant-potential
electrodes + IL EDL box); no MetalWalls twin binary, no LAMMPS ELECTRODE
package in the available `lmp`.
