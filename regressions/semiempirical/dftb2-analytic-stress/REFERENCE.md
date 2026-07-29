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
3. **Strain-FD identity:** analytic \(\sigma_{xx}\), \(\sigma_{yy}\) vs
   multi-SCC affine strain-FD of the same \(E\) (\(h=2\times10^{-4}\)) within
   **1×10⁻⁸ eV/Å³**. The stress is the pair virial of the same bond terms as
   the analytic force, and a converged SCC-DFTB2 gradient is stationary, so
   this is FD-exact to the strain-FD truncation error.

Offline pin note: `tools/tb_oracle/inputs/dftbplus/dftb2_analytic_stress_note.md`.

Measured after the gradient + double-counting fixes (2026-07-28):

```text
E = -91.70909683 eV
sigma = (-1.93986571e-03, -1.71915651e-03, 0.0,
          1.02900333e-03, -7.8e-20, 3.6e-20) eV/A^3
sigma_xx: analytic -1.93986571e-03  strain-FD -1.93986582e-03  |d| 1.147e-10
sigma_yy: analytic -1.71915651e-03  strain-FD -1.71915660e-03  |d| 8.937e-11
```

\(\sigma_{zz}=0\) exactly and the \(xz,yz\) shears are ~1e-20 — correct for a
planar molecule in the \(xy\) plane.

## Closed — the 1.3e-4 eV/Å³ residual was the force defect (fixed 2026-07-28)

The gate used to read **5×10⁻³ eV/Å³** and pass at 1.3e-4 / 1.0e-4, i.e. ~7%
relative, described as a "fixed-\(\Delta q\) vs re-SCC residual". It was not:
this stress is by construction the pair virial of the same bond terms as
`analytic_forces_restricted`, and it inherited that function's wrong-sign
\(\gamma\) gradient. With the sign corrected (and the SCC actually converged,
and the energy on the \(\mathrm{Tr}(\rho H^0) + E^{(2)} + E_{\mathrm{rep}}\)
route) the residual is **1.1e-10 / 8.9e-11 eV/Å³** — six orders of magnitude
tighter, and the gate is tightened to 1e-8 to match.

Both \(E\) and \(\sigma\) moved (the stress even changed sign): the old
numbers were a limit-cycle snapshot over a double-counting route that was wrong
at second order, differentiated by a gradient with a sign error. Neither the old
nor the new value is an oracle number — no DFTB+ binary exists on this host
(`module spider dftb` → nothing) and no DFTB+ identity is claimed.

## Gaps (still open)

1. **D3 dispersion is absent from the stress.** `DFTB2::d3` contributes to the
   energy and forces but not to `kStress` — the analytic pair virial does not
   see it. Setting `d3 = true` with `with_stress = true` silently yields an
   electronic+repulsive-only stress. A D3 virial is owed.
2. **No periodic case, and the self-image virial is missing.** This gate's
   geometry is a molecule in a 20 Å box, so only the \(a<b\) pair virial is
   exercised. The pair walk (`src/cpu/semiempirical/dftb2.cpp:1325`,
   `for (int b = a + 1; ...)`) skips self-images (\(a=a\), \(T\neq0\)). Those
   exert zero nuclear force — so forces are unaffected and
   `reg_dftb2-periodic-forces` is green — but they **do** carry a virial once a
   lattice image falls inside the SK / γ / repulsive cutoff, so under dense
   Γ-PBC (e.g. the 2.8 Å C-chain cell) the analytic stress is incomplete.
   A periodic strain-FD gate plus the self-image virial are owed together.

**Closed 2026-07-29:** the open-shell branch no longer falls back to multi-SCC
strain FD — `analytic_forces_restricted` / `analytic_forces_unrestricted` were
folded into one spin-resolved `analytic_forces_stress` that also emits the
virial, and this case's numbers are unchanged
(`regressions/semiempirical/dftb2-spin-analytic-stress`).

## Perf

`perf_exempt: true` — no DFTB+ binary / fair wall on pin host. Analytic path
replaces 12 multi-SCC energy evaluations (6 Voigt × ±h) with one post-SCC
pair contraction.
