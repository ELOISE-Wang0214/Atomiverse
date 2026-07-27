# Reference — H2 non-SCC DFTB0 energy

```text
oracle: dftb+
oracle_version: methodology / Hamiltonian twin (DFTB+ non-SCC path);
  offline DFTB+ binary unavailable on build host (2026-07-27).
  Numerical pins are the ATV DFTB0 closed form under the static
  mio-style tables in basis_data/dftb_mio_tables.h (not bit-identical
  to DFTB.org mio-1-1 .skf). True DFTB+ binary energy identity deferred
  to dftb2-dftbplus-energy-pin after SK fidelity.
```

## Method

Non-SCC DFTB / **DFTB0** (band structure of H0 + pair repulsive; no SCC):

1. Parameter set → static SK tables (`dftb_mio_tables.h`)
2. H0 / S construction (two-centre exponential / Gaussian radial forms)
3. Generalized eigen \(H_0 C = S C \varepsilon\)
4. Closed-shell occupation \(n_\mathrm{occ} = n_e/2\)
5. \(E = \sum_i 2\varepsilon_i + E_\mathrm{rep}\),  
   \(E_\mathrm{rep}=\sum_{A<B} A\,e^{-B R}\)

Family reference package: **DFTB+** non-SCC / DFTB0 electronic path
(Elstner et al. *Phys. Rev. B* **58**, 7260 (1998) lineage; DFTB+ docs:
non-SCC vs SCC-DFTB2). **Not** OpenMOPAC; **not** xtb.

## Data

In-tree static tables `include/atomiverse/cpu/semiempirical/basis_data/dftb_mio_tables.h`
for H/C/N/O. **Not** a full DFTB.org mio-1-1 dump.

## Geometry

Gas-phase H2, bond along \(x\), \(R = 0.74\) Å (and \(R = 3.0\) Å for
non-flatness). Cubic box 30 Å (isolation). Charge 0, singlet.

| Atom | x    | y | z |
|------|------|---|---|
| H    | 0.00 | 0 | 0 |
| H    | 0.74 | 0 | 0 |

## Pins

| Quantity | Value | Unit | Tol |
| -------- | ----- | ---- | --- |
| \(E(R=0.74\,\text{Å})\) | **−18.1376035790** | eV | \(1\times10^{-4}\) eV |
| \(|E(0.74)-E(3.0)|\) | \(> 10^{-6}\) | eV | pipeline not flat |
| Force skeleton at \(R=0.90\) Å | \(\|F\|^2 > 0\) | eV²/Å² | finite, non-zero |

## Perf

Offline DFTB+ binary was not available on the host; fair wall
\(T_\mathrm{atv}/T_\mathrm{oracle}\) not measured → **perf_exempt**.
