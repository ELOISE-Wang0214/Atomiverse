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

## Gaps — the 1.0 eV/Å tolerance is masking a defect (found 2026-07-28)

Gate 4 passes at \(|d| = 0.23\) eV/Å (\(F_{O,x}\)) and \(0.75\) eV/Å
(\(F_{H1,y}\)) against a tolerance of **1.0**, i.e. the gate is ~4x looser than
the error it is meant to catch. The standing explanation — "Mulliken charge
response beyond fixed-\(\Delta q\) HF" — does not hold: a converged SCC-DFTB2
gradient is **stationary**, so the analytic force is FD-exact and needs no
coupled-perturbed correction (this is why DFTB+ reproduces its own FD to
~1e-8).

Evidence that the charge-dependent term specifically is at fault:

| case | \(\Delta q\) | analytic vs re-SCC FD |
| ---- | ------------ | --------------------- |
| O₂ triplet (`reg_dftb2-spin-analytic-forces`) | 0 | **3.0e-07** eV/Å |
| H₂ (unit test) | 0 | 5e-03 eV/Å |
| H₂O (this case) | ≈ ±1.9 e | **0.23 / 0.75** eV/Å |

Every \(\Delta q = 0\) system is essentially exact; only the polar one is off.
The culprit is therefore the \(-\Delta q_A \Delta q_B\,
\partial\gamma_{AB}/\partial R\) contribution (and/or the shift term riding on
\(\partial S\)), not a missing response term.

The cheap explanations were ruled out by A/B on `dev` (2026-07-28), rebuilding
`analytic_forces_restricted` with the term perturbed:

| variant of the \(\gamma\) gradient term | \(\|d\|\) on \(F_{O,x}\) |
| --------------------------------------- | ------------------------ |
| as shipped, \(-\Delta q_A \Delta q_B \, d\gamma/dR\) | **0.228** eV/Å |
| term dropped (coefficient 0) | 1.020 eV/Å |
| sign flipped (\(+\Delta q_A \Delta q_B\)) | 2.268 eV/Å |

So the shipped sign and magnitude are the best of the three — it is **not** a
sign slip or a missing factor of 2, and the fix needs a real re-derivation.
Route: `/mol:fix` on the DFTB2 restricted gradient; do **not** widen the
tolerance further, and tighten it to the O₂ class (≤1e-5 eV/Å) once fixed.
Until then the `MioSketch` analytic forces are quantitatively unreliable for
polar systems, and `reg_dftb2-spin-analytic-forces` inherits the same formula.

## Perf

`perf_exempt: true` — no DFTB+ binary / fair wall on pin host. Analytic path
replaces \(6N\) full re-SCC force evaluations with one post-SCC contraction.
