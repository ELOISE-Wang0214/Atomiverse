# H2O SCC-DFTB2 (static mio-style SK)

## Method

Self-consistent-charge DFTB2 (Elstner et al., *Phys. Rev. B* **58**, 7260
(1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)).

## Data

In-tree static tables `basis_data/dftb_mio_tables.h` for H/C/N/O only.
**Not** a full DFTB.org mio-1-1 `.skf` dump — hoppings are mio-style
exponentials for the product path (no runtime parser).

## Geometry

Water, slightly off-equilibrium (Å):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

Cubic box 20 Å (isolated).

## Gates

- Finite, non-zero total energy [eV] on `frame.meta[kEnergy]`
- Finite analytic forces [eV/Å] with |F|² > 0
- One Cartesian component matches central FD to 0.5 eV/Å
