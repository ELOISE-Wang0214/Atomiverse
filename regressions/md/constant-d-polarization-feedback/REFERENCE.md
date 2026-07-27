# constant-d-polarization-feedback

oracle: metalwalls
oracle_version: literature — Zhang & Sprik, Phys. Rev. B 93, 144201 (2016);
  Cox & Sprik, J. Chem. Phys. 151, 064506 (2019). No MetalWalls binary is
  installed on this machine and mainline LAMMPS has no constant-D ensemble,
  so the pins are the closed-form constant-D expressions from the finite-field
  papers (the same E = D − 4πP electric-enthalpy functional MetalWalls
  implements for its finite-field route), evaluated with the project's
  Coulomb constant. perf_exempt: no fair twin binary.

## Physics pinned

Constant-D ensemble (Zhang–Sprik), z axis, point charges, eV/Å/e/fs units
with k_e = 1/(4πε₀) = 14.399645 eV·Å/e² (`units::coulomb`):

- itinerant polarization  P_z = Σ q_i z_i^unwrapped / V   [e/Å²]
- Maxwell field           E_z = D_z − 4π k_e P_z          [V/Å]
- electric enthalpy       U_D = V E_z² / (8π k_e)         [eV]
- force                   F_{i,z} = −∂U_D/∂z_i = q_i E_z  [eV/Å]

## Pins (derivation: tools/metalwalls_oracle/inputs/literature/constant_d_pins.py)

Geometry: ±1 e pair on z, separation d = 2 Å, cubic L = 10 Å (V = 1000 Å³).

Gate A (static, D_z = 0.5 V/Å):

| quantity | value | units |
| --- | --- | --- |
| P_z | 2.0e-3 | e/Å² |
| E_z | 1.3809744842959204e-1 | V/Å |
| U_D | 5.2696244279044818e-2 | eV |
| F_z(+q) | 1.3809744842959204e-1 | eV/Å |
| foothold twin E_z | 0.5 (= D_z, unscreened) | V/Å |

Gate B (dynamics, D_z = 0): itinerant-polarization harmonic well
U = 2π k_e q² d²/V (Cox & Sprik) → plasma-like oscillation of the pair,

- ω = sqrt(4π k_e q² / (μ V) · force_to_accel) = 5.9091697154893236e-2 fs⁻¹
  (period 106.329 fs), μ = 0.5 amu
- d(t) = d0 cos(ωt), d0 = 2 Å; pinned at t = 26.6 fs and t = 53.2 fs
  (tol 5e-5 Å; VV dt = 0.1 fs discretization error ≤ 4.6e-6 Å)
- energy conservation KE + U_D = U0 = 3.6190255157040796e-1 eV
  (checked at the half-period turning point, tol 1e-6 eV)
- self-consistency E_z(t) = −4π k_e P_z(t) at every probe

Gate C (vacuum-foothold twin, feedback off, D = 0): zero force — atoms
never move, E_z stays D_z. This is the queue observable "E_z response vs
vacuum-foothold twin": feedback screens/pushes, the foothold does not.

## Citations

- S. J. Cox and M. Sprik, "Finite field formalism for bulk electrolyte
  solutions", J. Chem. Phys. 151, 064506 (2019); theory follows
  Zhang & Sprik. (Primary constant-D electrolyte formalism.)
- C. Zhang and M. Sprik, "Finite field methods for the supercell modeling
  of charged insulator/electrolyte interfaces", Phys. Rev. B 93, 144201
  (2016). (Constant-D electric enthalpy F = E_PBC + V(D − 4πP)²/8π.)
- M. Stengel, N. A. Spaldin, D. Vanderbilt, "Electric displacement as the
  fundamental variable in electronic-structure calculations",
  Nat. Phys. 5, 304 (2009). (Origin of the constant-D functional.)

Runtime: ~1064 MD steps of a 2-atom system — well under 60 s.
