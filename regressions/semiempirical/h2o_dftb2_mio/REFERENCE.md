# H2O SCC-DFTB2 (static mio-style SK)

```text
oracle: literature
oracle_version: Elstner et al. PRB 58, 7260 (1998) — framework gate (not DFTB+ E pin)
```

## Method

Self-consistent-charge DFTB2 (Elstner et al., *Phys. Rev. B* **58**, 7260
(1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)).

## Data

In-tree static tables for H/C/N/O: H onsite/U + H–H from transcribed mio-1-1
(`dftb_mio_hh_sk.h`); other pairs exponential sketches (`dftb_mio_tables.h`).
No runtime `.skf` parser. Full multi-element mio E identity is a separate pin
(`dftb2-dftbplus-energy-pin` for H₂).

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
- One Cartesian component matches central FD to **1.5 eV/Å** (fixed-Δq HF vs
  re-SCC residual after mio-1-1 H–H / H U)
