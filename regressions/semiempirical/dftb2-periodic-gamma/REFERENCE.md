# dftb2-periodic-gamma — Γ-point periodic SCC-DFTB2

```text
oracle: literature
oracle_version: ATV-mio-style-sketch (Hamiltonian twin; not DFTB.org mio-1-1)
family_oracle: dftb+
```

**Note:** No offline DFTB+ binary was available on the build host. Absolute
energy is pinned to the ATV static mio-style SK Hamiltonian at the geometry
below (same hygiene as `h2o_dftb2_mio`). Family numerical oracle for this
method remains **DFTB+** once bit-identical SK tables ship
(`dftb2-dftbplus-energy-pin` / `dftb2-3ob-organic`).

## Method

Self-consistent-charge DFTB2 (Elstner et al., *Phys. Rev. B* **58**, 7260
(1998), DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260))
at the **Γ-point**: lattice-summed two-center \(H^0\) and \(S\), screened
Ohno–Klopman-like \(\gamma\), and exponential repulsive within the static
pair cutoff from `basis_data/dftb_mio_tables.h`. No k-mesh; no Ewald
long-range Coulomb (γ is real-space truncated at \(R_\mathrm{cut}\)).

## Data

In-tree static tables `basis_data/dftb_mio_tables.h` for H/C/N/O only.
**Not** a full DFTB.org mio-1-1 `.skf` dump.

## Geometry

1-D carbon dimer crystal (two C atoms per orthorhombic cell):

| Atom | x [Å] | y [Å] | z [Å] |
|------|-------|-------|-------|
| C    | 0.00  | 0.00  | 0.00  |
| C    | 1.40  | 0.00  | 0.00  |

Cell: orthorhombic \(L_x = 2.8\) Å, \(L_y = L_z = 12.0\) Å.

With \(r(\mathrm{C–C})=1.4\) Å and \(L_x=2.8\) Å, each carbon has
nearest-neighbor contacts along ±x from the partner and its periodic
images; lattice-sum H/S differs from a single min-image molecular dimer.

## SCC settings

- `scc_tol = 1e-10`
- `max_scc = 50` (default)
- `with_forces = false` (energy gate)

## Pin

| Quantity | Value | Unit | Tol |
|----------|-------|------|-----|
| Total energy \(E\) | **−113.5976598145** | eV | ±0.05 eV |

Also required: \(|E - E_\mathrm{mol}| > 0.05\) eV where \(E_\mathrm{mol}\) is
the same dimer in a 20 Å cubic box (molecular limit).

## Unmoved by the 2026-07-28 SCC repair

`E = −113.5976598145 eV` (Γ-PBC) and `−90.9958822966 eV` (isolated dimer) are
byte-unchanged across the charge-mixing, double-counting and γ-gradient fixes.
The homonuclear C chain has Δq ≡ 0 by symmetry, so E² = 0 and the spurious
Σ_A V_A q⁰_A term vanishes. Like `spin-dftb2-o2` this is a control, not a check
of the Δq-dependent machinery.

## Perf

```text
perf_exempt: true
perf_exempt_reason: no offline DFTB+ binary on host; fair wall deferred
```
