# Reference — CH4 SCC-DFTB2 with static 3ob-3-1 SK

## Oracle

```text
oracle: dftb+
oracle_version: 24.1
```

Offline pin only — ctest never launches DFTB+.

## Method

- **Hamiltonian:** SCC-DFTB2 (Elstner et al., *Phys. Rev. B* **58**, 7260
  (1998); DOI [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260)).
- **SK set:** DFTB.org **3ob-3-1** organic H/C/N/O (Gaus, Goez, Elstner,
  *J. Chem. Theory Comput.* **9**, 338 (2013); DOI
  [10.1021/ct300849w](https://doi.org/10.1021/ct300849w)). Static tables in
  `basis_data/dftb_3ob_tables.h` (CC-BY-SA 4.0 transcription; no runtime `.skf`).
- **Pin Hamiltonian twin:** DFTB+ `SCC = Yes`, **no** `ThirdOrder` (pure DFTB2
  with 3ob SK files). 3ob was trained for DFTB3; this case pins **SK fidelity**
  under DFTB2.

## Geometry

Methane, tetrahedral, \(r(\mathrm{C{-}H}) = 1.09\) Å:

| Atom | x [Å] | y [Å] | z [Å] |
|------|-------|-------|-------|
| C | 0 | 0 | 0 |
| H | \(s\) | \(s\) | \(s\) |
| H | \(s\) | \(-s\) | \(-s\) |
| H | \(-s\) | \(s\) | \(-s\) |
| H | \(-s\) | \(-s\) | \(s\) |

with \(s = 1.09 / \sqrt{3}\). Isolated (cubic box 20 Å, molecular path).

## Pin

| Quantity | Value | Units |
|----------|-------|-------|
| Total energy (DFTB+ 24.1) | **−3.2291103259** | Ha |
| Total energy (public ATV) | **−87.86856831** | eV |

Conversion: \(E[\mathrm{eV}] = E[\mathrm{Ha}] \times 27.211386245988\).

DFTB+ detailed components (documentation only):

| Term | Ha | eV |
|------|-----|-----|
| Electronic | −3.1685644194 | −86.2210 |
| SCC | 0.0011039499 | 0.0300 |
| Repulsive | −0.0605459065 | −1.6475 |
| **Total** | **−3.2291103259** | **−87.8686** |

## Gate

- Public API: `cpu::Evaluate` + `cpu::semiempirical::DFTB2` with
  `sk_set = DftbSkSet::ThreeOb`, `with_forces = false`.
- Assert `|E_ATV − (−87.86856831)| ≤ 0.05` eV (linear SK interp + Ohno–Klopman
  γ vs DFTB+ poly-interp / γ implementation residual).

## Perf (optional, not ctest)

DFTB+ wall for this molecule: ~0.12 s (single core, release). ATV expected
well under 10× (small dense eigenproblem).
