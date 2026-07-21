# Chapter 19. Large molecules and the Coulomb bottleneck

The Fock matrix has \(O(N^2)\) elements, but the two-electron integrals that
build it scale formally as \(O(N^4)\). That mismatch is the central cost of
large-molecule HF. This chapter follows successive reductions of the Coulomb
work: direct integrals with screening, density fitting, periodic Ewald/PME,
open-boundary multipoles, and continuous fast multipoles (CFMM) for the
Coulomb matrix \(J\).

---

## 19.1 Four-centre integrals and direct SCF

Storing the full four-centre tensor is usually impossible, so **direct SCF**
recomputes integrals each iteration. The Cauchy–Schwarz bound

\[
|(\mu\nu|\lambda\sigma)|
  \le \sqrt{(\mu\nu|\mu\nu)}\,\sqrt{(\lambda\sigma|\lambda\sigma)}
\]

allows whole shell pairs with exponentially small contributions to be
discarded, bringing the practical cost near \(O(N^2)\) at fixed density and
precision. Shell-block Schwarz screening and Rys quadrature in this engine
implement that path.

Screening does not remove the prefactor of repeatedly evaluating surviving
shells. Further speed-ups need a compressed representation of Coulomb—low-rank
fitting or multipole/grid methods for the far field.

---

## 19.2 Density fitting

Product densities
\(\rho_{\mu\nu}(\mathbf r)=\chi_\mu(\mathbf r)\chi_\nu(\mathbf r)\) are highly
linearly dependent. Expanding them in an auxiliary set of size \(O(N)\) and
fitting in the **Coulomb metric** (Vahtras–Almlöf–Feyereisen) yields

\[
(\mu\nu|\lambda\sigma)
  \approx \sum_{PQ}
  (\mu\nu|P)\,[V^{-1}]_{PQ}\,(Q|\lambda\sigma),
\qquad
V_{PQ}=(P|Q).
\]

The four-centre tensor never appears; three-centre quantities
\((\mu\nu|P)\) are stored instead. The Coulomb metric makes the energy error
second order in the fitting residual. Auxiliary bases used only for \(J\)
(JFIT) differ from those for \(J\) and \(K\) (JKFIT); exchange is more
sensitive to auxiliary quality. This engine ships a single **def2 universal
JKFIT** table (including \(g\) shells) for both `RI::J` and `RI::JK`.

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2TZVP;
rhf.ri = cpu::scf::RI::JK; // None | J | JK
rhf.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;
rhf.run(frame);
```

`RI::None` is the exact four-centre path. Unsupported combinations (for
example analytic forces with RI) raise errors rather than degrading silently.

For water/def2-TZVP, RI-JK and exact ERI differ by about
\(5\,\mu\mathrm{Ha}\) (\(\sim 2\,\mu\mathrm{Ha}\) per atom) and match PySCF
`density_fit()` in the same fitting space to \(\sim10^{-12}\,\mathrm{Ha}\).
Fitting errors are systematic shifts; they often cancel in energy differences.
Auxiliary shells must be counted in the spherical convention to match the
standard def2-JKFIT definition.

---

## 19.3 Periodic systems: Ewald and PME

The periodic Coulomb sum is conditionally convergent: a naive
\(\sum q_i q_j/r_{ij}\) depends on the order of summation. Ewald splits
\(1/r\) into short- and long-range parts,

\[
\frac1r
  = \frac{\mathrm{erfc}(\alpha r)}{r}
  + \frac{\mathrm{erf}(\alpha r)}{r},
\]

with the erfc term evaluated in real space (neighbour lists) and the erf term
as a Poisson problem in reciprocal space. Both halves must share the same
\(\alpha\).

Direct Ewald accumulates a structure factor for each \(\mathbf k\) and, at
fixed precision, approaches \(O(N^2)\). Smooth particle-mesh Ewald (SPME;
Essmann et al., 1995) spreads charges onto a grid with B-splines, performs a
3D FFT convolution, and interpolates forces back to particles, giving
\(O(N\log N)\) overall.

| Method | Scaling | Use |
|---|---|---|
| Direct Ewald | \(\sim O(N^2)\) | Small boxes, reference, no FFT |
| SPME | \(O(N\log N)\) | Production periodic MD |
| FMM / CFMM | \(\sim O(N)\) | Open-boundary molecular Coulomb (§19.4–19.5) |

Implementation is split: pure mathematics under
`include/atomiverse/longrange/` (no Frame, MPI, or FFT); CPU wrappers
`KSpaceEwald` / `KSpacePME` own grids, communication, and FFTW. The build
flag remains `-DATV_BUILD_KSPACE`. Electronic-structure code that needs the
operators should include only `atomiverse/longrange/*`.

```cpp
const double alpha = 0.3;  // 1/Å
const double cutoff = 8.0; // Å
sim.add_potential(cpu::CoulLong{alpha, cutoff});
sim.add(cpu::KSpacePME{alpha, /*grid_spacing=*/1.0, /*order=*/6});
```

The rock-salt Madelung constant \(M\approx 1.7476\) per ion pair is a standard
check of the reciprocal implementation (`tests/cpu/longrange/` and related
regressions).

---

## 19.4 Open boundaries: the fast multipole method

Without periodic images the far field can be compressed with multipole series.
For a cluster of charges observed at a well-separated point,

\[
\Phi(\mathbf r)
  = \sum_{n=0}^{p}\sum_{m=-n}^{n}
    M_n^m\, I_n^m(\mathbf r-\mathbf c)
  + O\!\Big(\frac1{c-1}\Big(\frac1c\Big)^p\Big),
\]

truncation error falls geometrically with separation ratio \(c\)
(Greengard–Rokhlin). On an adaptive octree the standard pipeline is P2M on
leaves, M2M upward, M2L between well-separated cells, L2L downward, L2P on
leaves, and direct near-field sums for adjacent leaves. M2M and L2L are exact
within order \(p\); truncation enters only through M2L.

PME and FMM are complementary: PME handles **periodic** long-range Coulomb;
FMM handles **open-boundary** far fields. They are not interchangeable.

```cpp
#include "atomiverse/longrange/fmm_tree.h"

ATV::longrange::FmmTree tree(/*order=*/8, /*max_per_leaf=*/8);
double e = tree.energy(positions, charges); // dimensionless Coulomb sum
```

Relative to direct \(O(N^2)\) sums on random and lattice charge sets, the
energy error decreases with \(p\); at \(p=16\) relative errors below
\(10^{-10}\) are typical (`regressions/longrange/fmm_point_charges`).
Scaling \(N\) to \(4N\) should give a wall-time ratio well below the quadratic
expectation of 16.

---

## 19.5 Continuous FMM and the HF Coulomb matrix

The HF Coulomb matrix \(J\) is built from Gaussian product distributions, not
point charges. White, Johnson, Gill, and Head-Gordon (1994) assign each
product (total exponent \(\zeta=\zeta_a+\zeta_b\)) an extent

\[
r_{\mathrm{ext}}
  = \sqrt{2/\zeta}\,\mathrm{erfc}^{-1}(\varepsilon_{\mathrm{thr}}).
\]

Multipoles are allowed only when the centre–centre distance exceeds the sum of
extents; otherwise analytic ERIs are used. The predicate lives in
`longrange/gaussian_extent.h`.

In HF, CFMM is a **third route for \(J\)** (exact four-centre, density fitting,
CFMM). **Exchange \(K\) never uses multipoles.**

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::STO3G;
rhf.coulomb = cpu::scf::Coulomb::CFMM;
rhf.cfmm_order = 10;
rhf.cfmm_threshold = 1e-9;
rhf.run(frame);
```

Default `Coulomb::Direct` matches the historical exact path. Setting
`Coulomb::CFMM` together with `ri = RI::J` or `RI::JK` is rejected.

Correctness must precede any timing claim: CFMM \(J\) should reproduce
Schwarz direct \(J\) within a stated tolerance. Molecular PySCF has no
FMM/CFMM, so the external reference can only pin the direct path; CFMM’s
oracle is this engine’s Direct. On water/STO-3G a single \(J\) agrees
elementwise to \(\le 10^{-8}\,\mathrm{Ha}\); on \(n\)-decane C₁₀H₂₂/STO-3G the
converged total energies differ by \(\le 10^{-7}\,\mathrm{Ha}\).

`cfmm_order` is the multipole truncation; `cfmm_threshold` sets extents and
the near/far partition. Small molecules at default thresholds may be entirely
near-field, so CFMM and Direct coincide numerically and any wall-time gain
appears only for larger systems. Scaling methodology is documented under
`benchmarks/scf/validation/cfmm_scaling/`.

---

## Summary of numerical checks

| System | Note |
|---|---|
| H₂O def2-TZVP | RI-JK vs exact ERI: \(\Delta E\sim 5\,\mu\mathrm{Ha}\) |
| Rock salt | Madelung \(\approx 1.7476\) per ion pair |
| Point-charge FMM | vs direct sum: \(\lesssim 10^{-10}\) at \(p=16\) |
| C₁₀H₂₂ STO-3G | CFMM vs Direct total energy \(\le 10^{-7}\,\mathrm{Ha}\) |

Broader external energy comparisons are collected in the
[appendix](appendix-pyscf-oracle.md).

### Further reading

- Jensen, *Introduction to Computational Chemistry*, SCF bottleneck and linear scaling.
- Weigend, *Phys. Chem. Chem. Phys.* **8**, 1057 (2006) (auxiliary bases).
- Essmann et al., *J. Chem. Phys.* **103**, 8577 (1995) (SPME).
- Greengard & Rokhlin, *J. Comput. Phys.* **73**, 325 (1987) (FMM).
- White et al., *Chem. Phys. Lett.* **230**, 8 (1994) (CFMM).
