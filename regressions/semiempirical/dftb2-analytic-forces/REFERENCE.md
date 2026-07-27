# Reference — DFTB2 analytic nuclear forces (H2O)

oracle: dftb+
oracle_version: formula-level (Hourahine et al. JCP 152, 124101 (2020); no DFTB+ binary on pin host)

## Method

SCC-DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260 (1998)):

\[
E = \sum_i n_i\varepsilon_i - \tfrac12\sum_{AB}\Delta q_A\gamma_{AB}\Delta q_B + E_{\mathrm{rep}}
\]

Nuclear forces are the **closed-form** gradient at the SCC density (Hellmann–Feynman
+ Pulay energy-weighted density + Hamiltonian shift on \(\partial S\) + fixed-\(\Delta q\)
\(-\partial E_{\mathrm{SCC}}\) + repulsive) — **not** a re-SCC finite difference.
This is the DFTB+ force module boundary.

## Hamiltonian / tables

- ATV class: `cpu::semiempirical::DFTB2`
- SK / Hubbard / repulsive: static `basis_data/dftb_mio_tables.h` (mio-**style**
  sketch for H/C/N/O). **Not** bit-identical to DFTB.org mio-1-1.

## Geometry (Å)

| atom | Z | x | y | z |
| ---- | - | - | - | - |
| O | 8 | 0.00 | 0.00 | 0.00 |
| H | 1 | 1.05 | 0.00 | 0.00 |
| H | 1 | −0.26 | 1.00 | 0.00 |

Isolated molecule, cubic box 20 Å. `scc_tol = 1e-10`.

## Pins (hard-coded gates)

1. **Energy:** finite, \(|E| > 10^{-3}\) eV.
2. **Forces:** all components finite; \(|F|^2 > 10^{-8}\) at off-eq geometry.
3. **Translation:** \(|\sum_i F_i| \le 10^{-6}\) eV/Å per Cartesian.
4. **FD identity:** analytic \(F\) vs central **re-SCC** FD of the same \(E\)
   (\(h=10^{-4}\) Å) on \(F_{O,x}\) and \(F_{H1,y}\) within **1.0 eV/Å**.
   (H2 matches to \(5\times10^{-3}\) eV/Å in unit tests; multi-element residual is
   Mulliken charge response beyond fixed-\(\Delta q\) HF on sketch SK.)

Offline pin note: `tools/tb_oracle/inputs/dftbplus/dftb2_analytic_forces_note.md`.

## Perf

`perf_exempt: true` — no DFTB+ binary / fair wall on pin host. Analytic path
replaces \(6N\) full re-SCC force evaluations with one post-SCC contraction.
