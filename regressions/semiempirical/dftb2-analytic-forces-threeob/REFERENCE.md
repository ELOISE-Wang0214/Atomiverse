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

Measured at merge onto `dev` (2026-07-28):

```text
E = -110.03663082 eV
F_O_x  analytic 1.52154493  re-SCC FD 2.36539390  |d| 8.438e-01 eV/A
F_H1_y analytic 0.26705128  re-SCC FD 0.03167127  |d| 2.354e-01 eV/A
```

## Gaps — the gate passes, the numbers do not agree

Gate 4's tolerance is **1.0 eV/Å** and the residuals are 0.84 and 0.24 eV/Å;
on `F_H1_y` that is an **8x relative** disagreement (0.267 vs 0.032). This is
not "closed-form vs re-SCC FD at stationarity" — a converged SCC-DFTB2 gradient
is stationary and therefore FD-exact.

It is the **same defect** already recorded for the MioSketch path in
`regressions/semiempirical/dftb2-analytic-forces/REFERENCE.md`, now reproduced
independently on the 3ob tables. The fixed-Δq expression

```text
dE/dR |_gamma  =  -dq_A dq_B  d(gamma_AB)/dR
```

now exists in three places — `analytic_forces_restricted`,
`analytic_forces_unrestricted` and `analytic_forces_3ob` — each with its own
γ implementation (`gamma_ab` / `gamma_ab_ha`), and every polar case is wrong
while every Δq = 0 case is exact:

| case | Δq | analytic vs re-SCC FD |
| ---- | -- | --------------------- |
| O₂ triplet (`reg_dftb2-spin-analytic-forces`) | 0 | 3.0e-07 eV/Å |
| H₂ (unit test) | 0 | 5e-03 eV/Å |
| H₂O MioSketch (`reg_dftb2-analytic-forces`) | ≈ ±1.9 e | 0.23 / 0.75 eV/Å |
| H₂O 3ob (this case) | ≠ 0 | **0.84 / 0.24** eV/Å |

Three independent implementations agreeing on the same wrong answer points at a
shared error in the **formula**, not a transcription slip in one of them. Route:
`/mol:fix` on the DFTB2 charge-dependent gradient term, covering all three
functions at once; tighten every gate to the Δq = 0 class (≤1e-5 eV/Å)
afterwards. Do not widen these tolerances.

Also owed: ThreeOb now has an analytic SK derivative path
(`fill_sk_block_deriv_3ob`, `rep_3ob_dr`), so its `with_stress` fallback to
multi-SCC affine strain-FD is no longer necessary — a closed-form 3ob pair
virial is derivable and not implemented.

## Perf

`perf_exempt: true` — no DFTB+ fair wall on pin host. Analytic path replaces
\(6N\) full re-SCC force evaluations with one post-SCC contraction.
