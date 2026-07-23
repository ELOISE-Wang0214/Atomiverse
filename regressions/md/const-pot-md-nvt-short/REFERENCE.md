# Reference — constant-potential MD NVT short trajectory

```text
oracle: metalwalls
oracle_version: algorithm Siepmann–Sprik / Reed / Scalfi BO (MetalWalls-class);
  no MetalWalls binary or LAMMPS ELECTRODE package available offline on worker
  (LAMMPS 30 Mar 2026 packages: EXTRA-*, KSPACE, MANYBODY — no ELECTRODE)
```

## Map slice

E1 product (capability map §1 P0): fluctuating electrode charges at fixed
potential, frozen walls, classical Coulomb MD. Composition on `cpu::md::MD`.

## Algorithm pin (offline, not a live binary)

Electrode charges \(Q\) minimise the Gaussian-charge energy (Siepmann & Sprik,
*J. Chem. Phys.* **102**, 511 (1995); Reed, Lanning, Madden, *J. Chem. Phys.*
**126**, 084704 (2007); Scalfi *et al.*, *Chem. Rev.* **122**, 10899 (2022);
MetalWalls JOSS 2020, 10.21105/joss.02373):

\[
A Q = \Psi - \Phi_{\mathrm{electrolyte}},\qquad
\sum_i Q_i = -Q_{\mathrm{electrolyte}}
\]

with Gaussian width \(w\):

\[
A_{ii} = k_e\sqrt{2/\pi}/w,\quad
A_{ij}=k_e\,\mathrm{erf}\!\big(r_{ij}/(w\sqrt2)\big)/r_{ij}\ (i\neq j).
\]

\(k_e = \texttt{units::coulomb} = 14.399645\,\mathrm{eV·Å}/e^2\).
\(\Phi\) is the min-image Gaussian-smeared potential of non-electrode charges.
Force path remains point-charge `CoulLong` + `KSpacePME` (Gaussian force
kernels are `gaussian-electrode-charge`).

## Geometry (regression)

| Species | type | role |
| --- | ---: | --- |
| Na⁺ | 0 | mobile electrolyte |
| Cl⁻ | 1 | mobile electrolyte |
| left wall | 2 | frozen electrode at \(\Psi_L\) |
| right wall | 3 | frozen electrode at \(\Psi_R\) |

- Box: orthorhombic \(L_x=L_y=12\) Å, \(L_z=20\) Å
- Walls: \(2\times2\) atoms at \(z=2\) Å (type 2) and \(z=18\) Å (type 3)
- Electrolyte: 4 Na⁺ + 4 Cl⁻ on a staggered mid-slab lattice
- \(\Psi_L=-0.5\,\mathrm{V}\), \(\Psi_R=+0.5\,\mathrm{V}\), \(w=0.5\) Å
- Soft LJ (`ε=0.05` eV, `σ=2.0` Å) + `CoulLong`/`PME` (`α=0.3` Å⁻¹, cutoff 6 Å)
- NVT: `VelocityVerlet` + `NoseHoover`, \(T=300\) K, `τ_T=100` fs, `dt=0.5` fs,
  400 steps

## Observables and tolerances

| Check | Bound | Rationale |
| --- | --- | --- |
| Vacuum-style polarity after first BO solve | \(Q_L<0<Q_R\) | Higher-\(\Psi\) plate is positive (capacitor) |
| Electroneutrality | \(\|Q_L+Q_R\|\le 10^{-8}\) e | Neutral electrolyte ⇒ electrode neutrality |
| Finite thermo | \(E,K\) finite every chunk | Smoke stability / no NaN |
| Finite charges | all \(q_i\) finite | Charge solve stability |
| Freeze | wall \(\Delta x\le 10^{-10}\) Å, \|v\|≤10^{-12} | FixFreeze |
| End polarity | \(Q_L<0<Q_R\) | Polarity survives short NVT |

Continuum order-of-magnitude (optional diagnostic, not a hard gate): vacuum
parallel-plate \(C\approx A/(4\pi k_e d)\) with \(A=L_x L_y\), \(d\approx 16\) Å,
\(\Delta\Psi=1\) V gives \(|Q|\sim 0.05\,e\) scale; atomistic Gaussian BO is
the same sign and similar magnitude on this toy grid.

## Perf

```text
perf_exempt: true
perf_exempt_reason: no MetalWalls binary and no LAMMPS ELECTRODE package on
  the offline worker; cannot fair-bench vs oracle twin for this MW-only path
```

## Units

eV, Å, fs, amu, eV/Å, K, charge e, potential V (= eV/e).
