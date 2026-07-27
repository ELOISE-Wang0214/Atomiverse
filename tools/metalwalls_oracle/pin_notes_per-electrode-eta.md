# Pin notes — per-electrode-eta

**Date:** 2026-07-27
**Slug:** `per-electrode-eta`
**oracle:** literature (MetalWalls per-species Gaussian width \(\eta_e\); no MW
binary on the agent host, no LAMMPS ELECTRODE per-electrode-eta twin available
offline)
**oracle_version:** Reed, Lanning & Madden, J. Chem. Phys. **126**, 084704
(2007), DOI 10.1063/1.2464084; Scalfi, Dufils, Reeves, Rotenberg & Salanne,
J. Chem. Phys. **153**, 174704 (2020), DOI 10.1063/5.0028232; Marin-Laflèche
et al. (MetalWalls), J. Open Source Softw. **5**(53), 2373 (2020),
DOI 10.21105/joss.02373

## Method

MetalWalls assigns each electrode species its own Gaussian charge width
\(\eta_e\) (per-`species` `charge gaussian` block in `runtime.inpt`). The
constant-potential elastance matrix elements are the published closed forms
(Reed–Lanning–Madden 2007; Scalfi et al. 2020; MW docs):

- self: \(A_{ii} = k\,\sqrt{2/\pi}\,\eta_i\)
- cross: \(A_{ij} = k\,\mathrm{erf}(\eta_{ij} r_{ij})/r_{ij}\),
  \(\eta_{ij} = \eta_i\eta_j/\sqrt{\eta_i^2+\eta_j^2}\)
- electrode–point-charge RHS: \(b_i = k\,q_j\,\mathrm{erf}(\eta_i r)/r\)

with \(k =\) `units::coulomb` \(= 14.399645\) eV·Å/e².

The pin is the exact constrained BO solve for a **two-site vacuum capacitor**
with dissimilar widths (electroneutral, \(Q_2=-Q_1\)):

\[
Q = \frac{\Delta\Psi}{A_{11}+A_{22}-2A_{12}}
\]

Inputs: \(\eta_1 = 1.979\) Å⁻¹ (LAMMPS ELECTRODE default / Ahrens-Iwers 2022),
\(\eta_2 = 0.955\) Å⁻¹ (Merlet/Salanne graphite value used in MW examples),
\(d = 4\) Å, \(\Delta\Psi = 1\) V.

\[
A_{11} = 22.737234510126967,\quad
A_{22} = 10.972237977347778,\quad
\eta_{12} = 0.860091344088193\ \text{Å}^{-1},
\]
\[
A_{12} = 3.599907138217925,\qquad
Q = \mathbf{0.037722100829787}\ e.
\]

Shared-η twins of the same geometry (ordering check):
\(Q(\eta{=}1.979\ \text{shared}) = 0.026126955854989\) e,
\(Q(\eta{=}0.955\ \text{shared}) = 0.067823...\) e, so
\(Q_{\mathrm{shared,1.979}} < Q_{\mathrm{mixed}} < Q_{\mathrm{shared,0.955}}\).

Recorded in `regressions/md/per-electrode-eta/REFERENCE.md`; the regression
also asserts the mixed-η multi-atom two-wall system (2D Gaussian xy lattice
images) stays electroneutral, keeps the sign of \(\Delta\Psi\), and lands
strictly between its two shared-η twins.

**ctest never launches an oracle binary.** The would-be MetalWalls twin input
is preserved at `inputs/metalwalls/per-electrode-eta/runtime.inpt` (not run —
no binary; kept so a host with MW installed can reproduce the pin).

## perf

`perf_exempt=true` — no runnable oracle binary implements this physics on the
host (MetalWalls not installed; LAMMPS ELECTRODE fix here has a single `eta`
per fix and no offline per-electrode-eta twin was runnable), so no fair
wall-clock ratio exists.
