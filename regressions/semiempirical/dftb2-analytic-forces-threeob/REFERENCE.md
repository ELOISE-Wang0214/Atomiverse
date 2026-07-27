# Reference — DFTB2 ThreeOb analytic nuclear forces (H2O)

oracle: dftb+
oracle_version: formula-level (Hourahine et al. JCP 152, 124101 (2020); Gaus
  et al. JCTC 9, 338 (2013) 3ob-3-1; no DFTB+ binary required on pin host)

## Method

SCC-DFTB2 on static **3ob-3-1** organic SK (`DftbSkSet::ThreeOb`):

\[
E = \mathrm{Tr}(P H_0) + \tfrac12\sum_{AB}\Delta q_A\gamma_{AB}\Delta q_B + E_{\mathrm{rep}}
\]

Nuclear forces are the **closed-form** gradient at the SCC density
(Hellmann–Feynman + Pulay energy-weighted density + Hamiltonian shift on
\(\partial S\) + fixed-\(\Delta q\) \(\gamma\) geometry terms + repulsive
spline derivative) — **not** a re-SCC finite difference.
This is the DFTB+ force module boundary applied to production 3ob tables.

## Hamiltonian / tables

- ATV class: `cpu::semiempirical::DFTB2`
- SK / Hubbard / repulsive: static `basis_data/dftb_3ob_tables.h` (DFTB.org
  3ob-3-1 H/C/N/O; Gaus et al. 2013). Same tables as `dftb2-3ob-organic`.

## Geometry (Å)

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| H | 1 | 1.05 | 0.00 | 0.00 |
| H | 1 | −0.26 | 1.00 | 0.00 |

Isolated molecule, cubic box 20 Å. `scc_tol = 1e-10`, `max_scc = 80`.

## Pins (hard-coded gates)

1. **Energy:** finite, \(|E| > 10^{-3}\) eV.
2. **Forces:** all components finite; \(|F|^2 > 10^{-8}\) at off-eq geometry.
3. **Translation:** \(|\sum_i F_i| \le 10^{-6}\) eV/Å per Cartesian.
4. **FD identity:** analytic \(F\) vs central **re-SCC** FD of the same \(E\)
   (\(h=10^{-4}\) Å) on \(F_{O,x}\) and \(F_{H1,y}\) within **1.0 eV/Å**.

Offline pin note:
`tools/tb_oracle/inputs/dftbplus/dftb2_analytic_forces_threeob_note.md`.

## Perf

`perf_exempt: true` — no DFTB+ fair wall on pin host. Analytic path replaces
\(6N\) full re-SCC force evaluations with one post-SCC contraction.
