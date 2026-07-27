# H2O SCC-DFTB3 (static mio-style SK)

```text
oracle: dftb+
oracle_version: formula-Gaus2011 (DFTB+ module boundary; no host DFTB+ binary)
```

## Method

Self-consistent-charge **DFTB3** (third-order charge fluctuation):

- Gaus et al., *J. Chem. Theory Comput.* **7**, 931 (2011),
  DOI [10.1021/ct100684s](https://doi.org/10.1021/ct100684s)
- Second-order SCC backbone: Elstner et al., *Phys. Rev. B* **58**, 7260
  (1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)

Family-locked numerical oracle package is **DFTB+** (never OpenMOPAC / xtb).

## Hamiltonian twin note

In-tree static tables `basis_data/dftb_mio_tables.h` for H/C/N/O only.
**Not** a full DFTB.org mio-1-1 `.skf` dump — hoppings / \(U^{\mathrm{d}}\)
are mio-style product-path sketches (no runtime parser). Absolute energy
identity with DFTB+ + production mio-1-1 is **deferred** to SK-fidelity
slugs (`dftb2-3ob-organic`, energy-pin). This regression pins the
**ATV DFTB3 Hamiltonian twin** energy at the stated geometry (deterministic
offline value measured under SE-only build), with formula lineage from
DFTB+/Gaus.

## Geometry

Water, slightly off-equilibrium (Å):

| Atom | x    | y    | z   |
|------|------|------|-----|
| O    | 0.00 | 0.00 | 0.0 |
| H    | 0.96 | 0.00 | 0.0 |
| H    |−0.24 | 0.93 | 0.0 |

Cubic box 20 Å (isolated).

## Parameters

- `scc_tol = 1e-10`, `max_scc = 50`, forces off (energy gate)
- Units: energy **eV**, positions **Å**

## Pin (ATV Hamiltonian twin)

| Quantity | Value | Unit | Tol |
| -------- | ----- | ---- | --- |
| \(E_{\mathrm{tot}}\) | **−84.561719022489** | eV | \(10^{-6}\) |

Measured under SE-only Release build, `scc_tol=1e-10`, forces off
(2026-07-27). Components at pin geometry (informational):
\(E^{(2)}\approx 26.14\) eV, \(E^{(3)}\approx 6.26\) eV (sketch SK + \(U^{\mathrm{d}}\)).

## Gates

- Finite, non-zero total energy [eV] on `frame.meta[kEnergy]`
- Hard-coded energy pin vs table above (tol \(10^{-6}\) eV)
- \(|E_{\mathrm{DFTB3}} - E_{\mathrm{DFTB2}}| > 10^{-6}\) eV (third-order active)

## Perf

`perf_exempt: true` — host has no DFTB+ binary for fair wall ratio; formula
oracle only. Revisit when offline DFTB+ pin infrastructure lands under
`tools/tb_oracle/inputs/dftbplus/`.
