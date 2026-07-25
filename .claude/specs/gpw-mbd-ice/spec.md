# Spec — `gpw-mbd-ice`

## Goal

Ship **Many-body dispersion (MBD@rsSCS)** as an additive free-atom (volume
ratio = 1) correction on Γ-point GPW, matching the offline **libMBD** pin on the
ice-relevant water-dimer geometry (same coordinates as `gpw-ts-ice` /
`gpw-pbe-d3bj-ice`).

## Method

- Ambrosetti, Reilly, DiStasio, Tkatchenko, *J. Chem. Phys.* **140**, 18A508
  (2014) — range-separated self-consistent screening + long-range MBD.
- libMBD: Hermann et al., *J. Chem. Phys.* **159**, 174802 (2023).
- Free-atom α₀ / C₆ / R_vdW from in-tree `basis_data/ts_tables.h` (same as TS).
- PBE damping β = 0.83, Fermi steepness a = 6 (libMBD defaults).
- Frequency quadrature: Gauss–Legendre n = 15 (libMBD pure-Python path).
- Isolated (no lattice) is the dual hard pin; optional orthorhombic periodic
  with Ewald + k-mesh when `cell` is set.

## Public API

- `cpu::gpw::GpwDft::mbd = true` (and `GPWEvaluator::mbd`)
- Mutually exclusive with `d3` / `d3bj` / `d4` / `ts`
- Energy on `GpwResult::e_disp` [Hartree], folded into `e_total`
- Forces: energy-only v1 (same policy as D4); no silent wrong forces

## Out of scope

- Hirshfeld density volumes (volume ratio ≠ 1)
- Analytic MBD forces / stress
- vdW-DF nonlocal correlation (`gpw-vdw-df-ice`)
- Live CP2K in ctest (CP2K stock `PAIR_POTENTIAL` has no MBD; libMBD is the
  offline oracle, same as TS)

## Acceptance

See `acceptance.md`.
