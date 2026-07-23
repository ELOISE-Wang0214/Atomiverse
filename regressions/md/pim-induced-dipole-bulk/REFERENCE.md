# Reference — bulk PIM self-consistent induced dipoles (E4)

oracle: metalwalls  
oracle_version: MetalWalls JOSS 5 2373 (2020) bulk PIM multipole form; Thole Chem. Phys. 59 341 (1981) exponential damping; offline analytic pin (no MetalWalls binary on agent host)  
perf_exempt: true  
perf_exempt_reason: No MetalWalls binary on PATH for a fair bulk-PIM wall twin; pin is the published multipole + Thole formula (same path as gaussian-electrode-charge literature pins).

## Physics

Classical Polarizable Ion Model (bulk, no electrodes):

\[
\boldsymbol{\mu}_i=\alpha_i\mathbf{E}_i^{\mathrm{tot}},\quad
\mathbf{E}_i^{\mathrm{tot}}
  =\mathbf{E}_i^{\mathrm{perm}}
  +\sum_{j\neq i}\mathbf{T}_{ij}\cdot\boldsymbol{\mu}_j
\]

\[
E_{\mathrm{pol}}
  =-\tfrac12\sum_i\boldsymbol{\mu}_i\cdot\mathbf{E}_i^{\mathrm{perm}}
  \quad\text{[eV via `units::coulomb`]}
\]

Permanent field from Frame point charges (min-image, cutoff). Mutual induction
uses the dipole tensor with **Thole** exponential screening
(\(u=a r/(\alpha_i\alpha_j)^{1/6}\), \(\lambda_3,\lambda_5,\lambda_7\)).
Optional **Tang–Toennies** \(f_n(br)\) multiplies the same kernels (`tt_b>0`).
Forces: Hellmann–Feynman at SCF (charge–dipole + dipole–dipole); permanent
charge–charge is **not** included (use `PairCoulCut` / PME separately).

This is the **MetalWalls bulk PIM** electrostatic polarization piece (map E4).
LAMMPS is **not** the oracle here (LAMMPS polarizable path is Drude /
permanent multipoles, not this SCF bulk PIM).

## Geometry — offline pin

Cubic box \(L=50\) Å (min-image ≡ free space at this spacing), cutoff 20 Å,
`thole_a = 0.39`, `tt_b = 0`.

| Site | type | \(q\) [e] | \(\alpha\) [Å³] | \((x,y,z)\) [Å] |
| ---- | ---- | --------- | --------------- | --------------- |
| A | 1 | +1 | 1.0 | (0, 0, 0) |
| B | 2 | −1 | 1.5 | (3, 0, 0) |

### Pin values (hard-coded in regression)

| Observable | Value |
| ---------- | ----- |
| \(E_{\mathrm{pol}}\) [eV] | **−0.12040946414707571** |
| \(F_{A,x}\) [eV/Å] | **+0.040341072868572643** |
| \(F_{B,x}\) [eV/Å] | **−0.040341072868572643** |
| \(F_{*,y}\), \(F_{*,z}\) | 0 |
| \(\mu_{A,x}\) [e·Å] | 0.082882545215708481 |
| \(\mu_{B,x}\) [e·Å] | 0.12342697961661356 |
| SCC iterations | 8 (tol \(10^{-12}\) e·Å) |

**Tol:** \(|E - E_{\mathrm{pin}}| < 10^{-9}\) eV;  
\(|F_x - F_{\mathrm{pin}}| < 10^{-8}\) eV/Å; transverse \(|F| < 10^{-12}\).

**Dual gate:** energy **and** forces (action–reaction; FD of \(E_{\mathrm{pol}}\)
matches analytic \(F\) in unit tests within \(5\times 10^{-6}\) eV/Å).

**API:** `cpu::Evaluate` + `cpu::PimInducedDipole{20.0, {0,1.0,1.5}, 0.39}`.

## Provenance

- Thole, *Chem. Phys.* **59**, 341 (1981) — exponential damping.  
- Applequist isotropic induced-dipole SCC.  
- MetalWalls JOSS **5**, 2373 (2020) — classical PIM at electrodes/bulk.  
- Unit closed-form undamped two-site SCC + one-pol force \(F=2k\alpha/r^5\).  
- Map: MetalWalls capability **E4** (`pim-induced-dipole-bulk`).

## Fair wall (perf gate)

No MetalWalls binary available on the agent host → **perf_exempt**.
When a MW twin appears, record \(T_{\mathrm{atv}}/T_{\mathrm{MW}}\) here and
require \(\le 10\).
