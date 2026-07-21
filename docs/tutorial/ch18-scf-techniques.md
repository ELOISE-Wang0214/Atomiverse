# Chapter 18. Convergence of the SCF

The Roothaan–Hall process

\[
F(P_k)\,C = SC\,\varepsilon
\quad\rightarrow\quad
P_{k+1}
\]

is a fixed-point iteration. Nothing guarantees that it converges. This chapter
treats three complementary controls, in the order they are usually applied:
a better initial guess, better density or Fock mixing, and—near
self-consistency—a second-order orbital update.

---

## 18.1 Why iterations fail

Two failure modes dominate.

**Charge oscillation.** For two nearly degenerate configurations \(A\) and
\(B\), the Fock operator built from \(P_A\) prefers \(B\), and that from
\(P_B\) prefers \(A\). The density alternates with period two. Small
HOMO–LUMO gaps, stretched bonds, and diradical character favour this
behaviour.

**Overshoot.** Each update is too large, so the new density is farther from
the fixed point than the old one. A bare-nuclear (Hcore) guess is a common
trigger: without electronic screening the first iteration pushes charge too
strongly toward the nuclei.

Two classical remedies appear earlier in the tutorial.

- **Density damping**,
  \(P\leftarrow(1-\alpha)P_{\mathrm{new}}+\alpha P_{\mathrm{old}}\)
  (member `damping_factor`; \(\alpha\) is the fraction of the *old* density
  retained). Oscillation is reduced at the cost of slower progress.
- **Level shift** (Saunders–Hillier): virtual levels are raised by
  \(\lambda\). Larger \(\lambda\) is more stable and slower. After
  convergence the engine re-diagonalises the unshifted Fock matrix so that
  the shift does not remain in the final eigenvalues.

A systematic recipe is staged: a good guess (§18.2), a far-from-SCF mixer
(EDIIS/ADIIS, §18.3), then CDIIS or a second-order method near convergence
(§18.3–18.4).

---

## 18.2 Initial guesses

SCF is local. The guess selects which basin is entered and how far the
iteration starts from a stationary point. In order of increasing physical
content:

| Guess | Role |
|---|---|
| Hcore | Diagonalise \(T+V_{ne}\). No screening; poorest, always available. |
| MINAO / SAD | Superposition of atomic densities; default. Sensible density, not an eigenstate of any Fock operator. |
| Extended Hückel | Empirical \(F_{\mu\nu}\) (VSIP and Wolfsberg–Helmholz \(K=1.75\)); one diagonalisation yields orbitals. |
| SAP | Superposition of atomic potentials; \(F=T+V_{\mathrm{SAP}}\); density-free. |

Hückel tables cover H–Ar; unsupported elements raise an error. SAP covers
\(Z\le 118\). For open shells a spin-symmetric guess cannot by itself produce
spin polarisation; symmetry breaking still relies on mixing and occupation.

```cpp
cpu::scf::RHF rhf;
rhf.basis = cpu::scf::BasisSet::Def2SVP;
rhf.init_guess = cpu::scf::InitGuess::Hueckel; // or SAP, MINAO, Hcore
rhf.run(frame);
```

For water/def2-SVP RHF at \(10^{-9}\) tolerance, Hückel and SAP need fewer
iterations than Hcore (about 14 versus 16). Different codes parameterise
Hückel differently; compare relative ordering, not absolute iteration counts.

---

## 18.3 Mixing: CDIIS, EDIIS, and ADIIS

**CDIIS** (Pulay) forms an unconstrained least-squares combination of past
Fock matrices from commutator error vectors
\(e_i=X^\top(F_i P_i S-S P_i F_i)X\). Near self-consistency the error is
approximately linear in the density residual and convergence is nearly
quadratic. Far from SCF the linearisation fails; large coefficients can
drive the density out of the physical region.

**EDIIS** (Kudin–Scuseria–Cancès) minimises a quadratic energy model on the
simplex \(\sum c_i=1\), \(c_i\ge 0\). The new density is a convex combination
of history densities and resists overshoot. **ADIIS** (Hu–Yang) uses a
related augmented Roothaan–Hall model that does not need an explicit energy
channel. Both are for the **far** region. Near convergence the optimum
collapses to the latest step and they stall around errors of \(10^{-4}\).
The intended use is a permanent switch: EDIIS/ADIIS far away, then CDIIS
once \(\max|e|\) drops below a threshold.

```cpp
cpu::scf::UHF uhf(5, 5);
uhf.mixer = cpu::scf::MixerKind::ADIIS_CDIIS;
uhf.mixer_switch_tol = 1e-2;
uhf.run(frame);
```

Example: water stretched by a factor 2.6, def2-SVP UHF, Hcore guess. The
spin-restricted solution is unstable; the true minimum is a broken-symmetry
diradical (\(\langle S^2\rangle\approx 1.97\)). Pure CDIIS oscillates for 60
steps; ADIIS→CDIIS converges in about 38 steps to
\(-75.718966\,\mathrm{Ha}\) (within \(10^{-6}\,\mathrm{Ha}\) of PySCF).

In this engine EDIIS energies lag the density by one iteration—an acceptable
far-region approximation. ADIIS avoids that lag. The near-region switch is
controlled by `mixer_switch_tol` (default \(10^{-2}\)).

---

## 18.4 Second-order SCF

The previous sections improve the fixed-point *path*. A qualitative change is
to minimise energy on the orbital manifold. Updates that preserve
\(C^\top S C=I\) can be written

\[
C(\kappa)=C\,e^{\kappa},
\qquad
\kappa^\top=-\kappa,
\]

with independent parameters only in the occupied–virtual block
\(\kappa_{ai}\). Expanding the energy to second order gives
\(E_0+g^\top\kappa+\frac12\kappa^\top H\kappa+\cdots\). For closed-shell real
orbitals \(g_{ai}=4F_{ai}\); the Hessian is the same \({}^1(A+B)\) operator
that appears in stability analysis. The Newton step \(\kappa=-H^{-1}g\)
converges quadratically near a minimum.

Building \(H\) explicitly is expensive, and far from the minimum \(H\) need
not be positive definite. **CIAH** (co-iterative augmented Hessian; Sun)
solves

\[
\begin{pmatrix}0&g^\top\\ g&H\end{pmatrix}
\begin{pmatrix}1\\ \kappa\end{pmatrix}
=\lambda
\begin{pmatrix}1\\ \kappa\end{pmatrix}
\]

for the lowest root. Davidson iteration needs only products \(H\cdot v\),
which reuses the stability \(\sigma\)-vector. A trust radius bounds
\(\|\kappa\|_2\); energy increases trigger a reduced trust radius and
rollback. The exponential map is evaluated with scaling-and-squaring Taylor
series so orthonormality is retained.

A Newton macro-iteration is more expensive than one CDIIS extrapolation
(occupied–virtual integral transformation). On well-behaved molecules CDIIS
suffices. On pathological cases the same staging as for mixers is
recommended: **DIIS first, SOSCF later**.

```cpp
cpu::scf::UHF uhf(5, 5);
uhf.soscf = true;
uhf.soscf_start = 0.25;   // switch after max|DIIS error| falls below this
uhf.soscf_max_micro = 12;
uhf.soscf_trust = 0.4;
uhf.run(frame);
```

Switching is not allowed before iteration 1: iteration 0 performs a Roothaan
diagonalisation so that density-only guesses are consistent with the
coefficients. After convergence an unshifted Fock matrix is diagonalised once
more to obtain canonical orbital energies (eigenvalues during SOSCF are only
semi-canonical).

For the stretched water example, Newton engages near step 4 and reaches the
broken-symmetry minimum in about nine macro-iterations, agreeing with PySCF
`mf.newton()` to better than \(10^{-11}\,\mathrm{Ha}\). For N₂ at
\(R=2.5\,\mathrm{Å}\)/6-31G RHF, CDIIS already converges; the Newton phase
of SOSCF is about four steps, comparable to PySCF’s macro-iteration count.

SOSCF is implemented for HF (RHF/UHF) only. Kohn–Sham needs the XC second
derivative kernel; enabling `soscf` on RKS/UKS is rejected.

---

## Numerical checks

| System | Observation |
|---|---|
| H₂O def2-SVP RHF | Hückel/SAP need fewer iterations than Hcore |
| 2.6× H₂O def2-SVP UHF | ADIIS→CDIIS converges; pure CDIIS stalls |
| Same + SOSCF | Energy matches PySCF newton to \(\lesssim10^{-8}\,\mathrm{Ha}\) |
| N₂ \(R=2.5\,\mathrm{Å}\) 6-31G RHF + SOSCF | Energy matches PySCF newton |

### Further reading

- Jensen, *Introduction to Computational Chemistry*, SCF techniques.
- Pulay, *Chem. Phys. Lett.* **73**, 393 (1980) (DIIS).
- Kudin, Scuseria, Cancès, *J. Chem. Phys.* **116**, 8255 (2002) (EDIIS).
- Hu & Yang, *J. Chem. Phys.* **132**, 054109 (2010) (ADIIS).
- Sun, arXiv:1610.08423 (CIAH).
- Helgaker, Jørgensen, Olsen, *Molecular Electronic-Structure Theory*, Ch. 10.
