# Spec: dftb2-3ob-organic

## Goal

Ship **production 3ob-style static SK tables** for organic H/C/N/O under
SCC-DFTB2, beyond the in-tree mio-style exponential sketch. Pin total energy
against offline **DFTB+** (family oracle).

## Method / oracle

| Field | Value |
| ----- | ----- |
| method | `dftb2` (SCC-DFTB2) |
| SK set | DFTB.org **3ob-3-1** (organic CHNO), static C++ tables |
| oracle | **dftb+** |
| oracle_version | 24.1 |
| system | methane CH4 (small organic) |
| quantity | total energy [eV] |

## Public API

- `cpu::semiempirical::DFTB2` gains public member
  `DftbSkSet sk_set` with values `MioSketch` (default, existing path) and
  `ThreeOb` (new static 3ob tables).
- No `*Config` bag, no `make_*` / `build_*` factory.
- Units: energy eV, positions Å, forces eV/Å.

## Data

- New header `basis_data/dftb_3ob_tables.h`: transcribed radial H/S + repulsive
  spline for all ordered pairs among H/C/N/O (CC-BY-SA 4.0, DFTB.org 3ob-3-1).
- No runtime `.skf` parser.

## Implementation notes

- 3ob path: linear SK interpolation on Bohr grid; spline repulsive; Ohno–Klopman
  γ; energy `Tr(P H0) + E_scc + E_rep` (DFTB2 double-counting form).
- Hamiltonian twin for pin: DFTB+ `SCC=Yes`, **no** ThirdOrder (pure DFTB2 with
  3ob SK files). 3ob was trained for DFTB3; this slug pins **SK fidelity** under
  DFTB2, not full DFTB3 product (see `dftb3-water-mio`).

## Acceptance

1. Unit test: 3ob tables cover H/C/N/O; CH4 energy finite; unsupported Z throws.
2. Regression `reg_dftb2-3ob-organic`: CH4 Evaluate+DFTB2(`ThreeOb`) energy
   within stated tol of DFTB+ pin in `REFERENCE.md`.
3. Existing mio sketch path remains default; `h2o_dftb2_mio` not broken.

## Out of scope

- Analytic DFTB forces (FD ok); DFTB3 third-order; spin; periodic; GFN-xTB.
