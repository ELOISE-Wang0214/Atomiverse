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

Measured at merge onto `dev` (2026-07-28), on dev's mio-1-1 H–H tables:

```text
E = -44.60971806 eV
sigma = (1.76174259e-03, 1.66840261e-03, 0.0,
        -4.62267469e-04, -1.7e-18, 1.2e-18) eV/A^3
sigma_xx: analytic 1.76174259e-03  strain-FD 1.89593569e-03  |d| 1.342e-04
sigma_yy: analytic 1.66840261e-03  strain-FD 1.77276608e-03  |d| 1.044e-04
```

\(\sigma_{zz}=0\) exactly and the \(xz,yz\) shears are ~1e-18 — correct for a
planar molecule in the \(xy\) plane.

## Gaps

1. **The 1.3e-4 eV/Å³ residual is the force defect, not strain-FD noise.**
   This stress is by construction the pair virial of the *same* bond terms as
   `analytic_forces_restricted`, so it inherits that function's disagreement
   with FD — 0.23 / 0.75 eV/Å on this very geometry (see the Gaps section of
   `regressions/semiempirical/dftb2-analytic-forces/REFERENCE.md`, where an A/B
   rules out a sign slip or missing factor). Relative error here is ~7%. Do not
   read the passing gate as "analytic stress verified": it is verified only to
   the accuracy of a gradient that is itself known to be wrong for polar
   systems. Tighten this gate to ~1e-8 once the gradient is fixed.
2. **D3 dispersion is absent from the stress.** `DFTB2::d3` contributes to the
   energy and forces but not to `kStress` — neither the analytic pair virial nor
   the open-shell strain-FD path (which differentiates `scc_total`, excluding
   dispersion) sees it. Setting `d3 = true` with `with_stress = true` silently
   yields an electronic+repulsive-only stress. A D3 virial is owed.
3. **No periodic case.** Both this gate and the geometry are a molecule in a
   20 Å box, so the stress is never checked against the Γ-point lattice-summed
   energy path, where the min-image force/virial expression is not the gradient
   of the lattice-summed energy (`dftb2-periodic-forces` remains ○).

## Perf

`perf_exempt: true` — no DFTB+ binary / fair wall on pin host. Analytic path
replaces 12 multi-SCC energy evaluations (6 Voigt × ±h) with one post-SCC
pair contraction.
