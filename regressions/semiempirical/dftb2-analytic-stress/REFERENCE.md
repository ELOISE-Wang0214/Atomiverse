# Reference — DFTB2 analytic cell stress (H2O)

oracle: dftb+
oracle_version: formula-level (Hourahine et al. JCP 152, 124101 (2020); no DFTB+ binary on pin host)

## Method

SCC-DFTB2 (Elstner et al., Phys. Rev. B **58**, 7260 (1998)):

\[
E = \sum_i n_i\varepsilon_i - \tfrac12\sum_{AB}\Delta q_A\gamma_{AB}\Delta q_B + E_{\mathrm{rep}}
\]

Intensive cell stress (mechanical / DFTB+ / PW–GPW convention):

\[
\sigma_{\alpha\beta} = -\frac{1}{\Omega}\frac{\partial E}{\partial\varepsilon_{\alpha\beta}}
\]

closed-form at the SCC density via the **pair virial** of the same analytic
nuclear force terms (Hellmann–Feynman + Pulay + fixed-\(\Delta q\) SCC \(\gamma\)
+ repulsive) — **not** a multi-SCC strain finite difference on the product path.

Voigt order on Frame `kStress`: `(xx, yy, zz, xy, xz, yz)` in **eV/Å³**.

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

Isolated molecule, cubic box 20 Å. `scc_tol = 1e-10`. `with_stress = true`.

## Pins (hard-coded gates)

1. **Energy:** finite, \(|E| > 10^{-3}\) eV.
2. **Stress:** all six Voigt components finite.
3. **Strain-FD identity:** analytic \(\sigma_{xx}\), \(\sigma_{yy}\) vs multi-SCC
   affine strain-FD of the same \(E\) (\(h=2\times10^{-4}\)) within **5×10⁻³ eV/Å³**
   (fixed-\(\Delta q\) vs re-SCC residual on multi-element sketch SK; H2 unit test
   is tighter).

Offline pin note: `tools/tb_oracle/inputs/dftbplus/dftb2_analytic_stress_note.md`.

## Perf

`perf_exempt: true` — no DFTB+ binary / fair wall on pin host. Analytic path
replaces 12 multi-SCC energy evaluations (6 Voigt × ±h) with one post-SCC
pair contraction.
