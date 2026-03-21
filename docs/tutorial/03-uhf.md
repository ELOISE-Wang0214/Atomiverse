# Unrestricted Hartree-Fock

> **Prerequisite:** None. This tutorial assumes C++ and basic linear algebra
> (matrices, eigenvalues) but no quantum chemistry.
>
> **Related:** [`docs/tutorial/02-rhf.md`](02-rhf.md) covers the restricted (closed-shell)
> case. Reading it first is helpful but not required.

---

## Table of Contents

1. [What problem are we solving?](#1-what-problem-are-we-solving)
2. [The wavefunction and mean-field idea](#2-the-wavefunction-and-mean-field-idea)
3. [Restricted HF and its fatal flaw](#3-restricted-hf-and-its-fatal-flaw)
4. [Unrestricted HF: independent spin channels](#4-unrestricted-hf-independent-spin-channels)
5. [The SCF algorithm step by step](#5-the-scf-algorithm-step-by-step)
6. [Spin contamination](#6-spin-contamination)
7. [Atomiverse API walkthrough](#7-atomiverse-api-walkthrough)
8. [Key code paths annotated](#8-key-code-paths-annotated)
9. [End-to-end example and plots](#9-end-to-end-example-and-plots)
10. [Interpreting the output](#10-interpreting-the-output)
11. [Molecular visualization with MolRec](#11-molecular-visualization-with-molrec)
12. [Common pitfalls](#12-common-pitfalls)

---

## 1. What problem are we solving?

A molecule is made of nuclei (protons and neutrons) and electrons. The nuclei are heavy
and slow; in the Born-Oppenheimer approximation we treat them as fixed points in space.
The electrons are light and fast — they rearrange almost instantly when anything moves.

The fundamental question of quantum chemistry is: **given a fixed arrangement of nuclei,
what is the lowest-energy state of the electrons?**

Why do we care about the energy? Because forces on nuclei are minus the gradient of the
energy with respect to nuclear positions. If you want to know whether a molecule is
stable, how it vibrates, or how it reacts, you need the energy surface. The ground-state
energy — the lowest possible — is what a real molecule settles into at low temperature,
and it is what UHF computes.

In Atomiverse, the answer is reported in **Hartree** (1 Hartree = 27.211 eV = 627.5 kcal/mol).
For reference, the hydrogen atom has a ground-state energy of \f$-0.5\f$ Hartree exactly,
or \f$-0.4666\f$ Hartree in the small STO-3G basis set we use here. Negative means the
electron is bound; more negative means more stable.

---

## 2. The wavefunction and mean-field idea

### The wavefunction

Quantum mechanics describes \f$N\f$ electrons through a *wavefunction*
\f$\Psi(\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_N)\f$, a complex-valued function of
all \f$N\f$ electron positions. Its square \f$|\Psi|^2\f$ is the joint probability density of
finding the electrons at those positions simultaneously.

For more than two electrons, \f$\Psi\f$ lives in a \f$3N\f$-dimensional space and is completely
intractable to represent exactly. Hartree-Fock is the simplest systematic approximation
that makes the problem tractable.

### Mean-field approximation

The expensive part of the exact calculation is that every electron interacts with every
other electron simultaneously. Hartree-Fock replaces this with a **mean-field**: each
electron is treated as if it moves in an average electric field created by all the other
electrons. This decouples the electrons and turns a single intractable \f$3N\f$-dimensional
problem into \f$N\f$ independent three-dimensional problems.

The price of this approximation is that the average field depends on the orbitals, and
the orbitals depend on the average field. The calculation must therefore be done
**self-consistently** (SCF = Self-Consistent Field): guess the orbitals, compute the
field, recompute the orbitals, update the field, repeat until nothing changes.

### Basis functions and matrices

In practice, orbitals are expanded in a finite set of **basis functions**
\f$\{\chi_1, \chi_2, \ldots, \chi_K\}\f$ (Gaussian functions centered on nuclei).
An orbital \f$\varphi_i\f$ is then a linear combination:

\f[
\varphi_i(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu i}\, \chi_\mu(\mathbf{r})
\f]

where \f$C_{\mu i}\f$ are the coefficients we solve for. This turns the
differential-equation problem into a matrix problem of size \f$K \times K\f$.

The three key matrices are:

- **\f$S\f$** (overlap, \f$K \times K\f$): \f$S_{\mu\nu} = \int \chi_\mu(\mathbf{r})\, \chi_\nu(\mathbf{r})\, d\mathbf{r}\f$. Measures how much basis functions overlap. The identity matrix would mean orthogonal basis functions; \f$S\f$ is generally not the identity.

- **\f$H^\text{core}\f$** (core Hamiltonian, \f$K \times K\f$): kinetic energy of each electron plus its attraction to all nuclei, ignoring electron-electron repulsion entirely.

- **ERI** (electron repulsion integrals, \f$K^4\f$): \f$(\mu\nu|\lambda\sigma) = \iint \chi_\mu(1)\, \chi_\nu(1)\, r_{12}^{-1}\, \chi_\lambda(2)\, \chi_\sigma(2)\, dV_1\, dV_2\f$. The four-index tensor encoding all pairwise electron-electron Coulomb repulsion.

In Atomiverse, these three objects are pre-computed and packaged into an `IntegralPack`
by calling `cpu::scf::build_integrals_sto3g(frame)`. The HF solver takes this pack and
does everything else.

---

## 3. Restricted HF and its fatal flaw

### What RHF does

Restricted Hartree-Fock (RHF) makes one extra simplifying assumption: **every electron
is paired**. Each spatial orbital \f$\varphi_i\f$ is occupied by exactly two electrons — one
spin-up (\f$\alpha\f$) and one spin-down (\f$\beta\f$) — sharing the same spatial wavefunction.

With \f$N\f$ electrons, \f$N/2\f$ orbitals are occupied. The density matrix is:

\f[
P_{\mu\nu} = 2\sum_{i=1}^{N/2} C_{\mu i}\, C_{\nu i} \qquad \text{(RHF)}
\f]

The factor of 2 accounts for the two electrons per orbital. The Fock matrix — the
effective Hamiltonian for one electron in the mean field of all others — is:

\f[
F_{\mu\nu} = H^\text{core}_{\mu\nu} + G_{\mu\nu},
\qquad G_{\mu\nu} = \sum_{\lambda\sigma} P_{\lambda\sigma}
\left[(\mu\nu|\lambda\sigma) - \tfrac{1}{2}(\mu\lambda|\nu\sigma)\right]
\f]

where \f$(\mu\nu|\lambda\sigma)\f$ is the Coulomb repulsion integral and
\f$(\mu\lambda|\nu\sigma)\f$ is the exchange integral. Solving the Roothaan-Hall equations
\f$F\mathbf{C} = S\mathbf{C}\boldsymbol{\varepsilon}\f$ gives the orbitals and their energies.

### The flaw: bond dissociation

RHF works well for molecules near equilibrium geometry, where electrons naturally pair
up. It breaks down catastrophically for **bond dissociation** — pulling a molecule apart.

Consider H\f$_2\f$ at a large bond length \f$R\f$. The physically correct answer is two isolated
hydrogen atoms, each with one electron of well-defined spin. But RHF insists that the two
electrons share the same delocalized bonding orbital
\f$(\varphi_A + \varphi_B)/\sqrt{2}\f$. Even at \f$R = 10\f$ Å, RHF keeps both electrons in
this orbital, so the wavefunction contains equal weight for "both electrons on atom A"
and "both electrons on atom B". Those ionic contributions are physically wrong at large
\f$R\f$ but RHF cannot remove them.

The result: in STO-3G, the RHF energy at large \f$R\f$ plateaus at approximately \f$-0.71\f$
Hartree instead of the correct dissociation limit \f$2E_\text{H} = 2 \times (-0.4666) = -0.9332\f$ Hartree.
That is an error of \f$\sim\!140\f$ kcal/mol — not a minor correction, but a qualitatively
wrong description of chemistry.

RHF also cannot handle any system with an odd number of electrons. A hydrogen atom has
one electron; there is no partner for it. A hydroxyl radical (OH) has 9 electrons, one
unpaired. All such *open-shell* systems are unreachable by RHF.

---

## 4. Unrestricted HF: independent spin channels

### The key idea

UHF lifts the RHF constraint by giving \f$\alpha\f$ and \f$\beta\f$ electrons their own
independent spatial orbitals. Instead of one set of coefficients \f$\mathbf{C}\f$, UHF
maintains two:

- \f$\mathbf{C}^\alpha\f$ — MO coefficients for spin-up (\f$\alpha\f$) electrons
- \f$\mathbf{C}^\beta\f$  — MO coefficients for spin-down (\f$\beta\f$) electrons

These can be completely different. At equilibrium H\f$_2\f$ they happen to be the same; at
stretched H\f$_2\f$ they diverge, with \f$\alpha\f$ localizing on one atom and \f$\beta\f$ on the
other.

### Density matrices

With \f$N_\alpha\f$ alpha electrons and \f$N_\beta\f$ beta electrons (\f$N_\alpha \geq N_\beta\f$ by
convention):

\f[
P^\alpha_{\mu\nu} = \sum_{i=1}^{N_\alpha} C^\alpha_{\mu i}\, C^\alpha_{\nu i}, \qquad
  P^\beta_{\mu\nu}  = \sum_{j=1}^{N_\beta}  C^\beta_{\mu j}\,  C^\beta_{\nu j}
\f]

\f[
P_{\mu\nu} = P^\alpha_{\mu\nu} + P^\beta_{\mu\nu} \qquad \text{(total density)}
\f]

Note the absence of the factor of 2 that appears in the RHF density: each occupied
orbital contributes exactly one electron per spin channel.

The total density \f$P\f$ captures the average position of all electrons regardless of spin.
It enters the Coulomb part of the interaction, which is spin-blind — an electron is
repelled by all other electrons, not just those of opposite spin.

### Fock matrices

Each spin channel gets its own Fock matrix:

\f[
G^\alpha_{\mu\nu} = \sum_{\lambda\sigma} P_{\lambda\sigma}(\mu\nu|\lambda\sigma)
  - P^\alpha_{\lambda\sigma}(\mu\lambda|\nu\sigma)
\f]

\f[
G^\beta_{\mu\nu}  = \sum_{\lambda\sigma} P_{\lambda\sigma}(\mu\nu|\lambda\sigma)
  - P^\beta_{\lambda\sigma}(\mu\lambda|\nu\sigma)
\f]

\f[
F^\alpha = H^\text{core} + G^\alpha, \qquad F^\beta = H^\text{core} + G^\beta
\f]

The Coulomb term \f$\sum_{\lambda\sigma} P_{\lambda\sigma}(\mu\nu|\lambda\sigma)\f$ is the
same for both channels: it depends on the total density, so an \f$\alpha\f$ electron feels
repulsion from all electrons. The exchange term \f$P^\sigma_{\lambda\sigma}(\mu\lambda|\nu\sigma)\f$
is **spin-specific**: quantum mechanics (Pauli exclusion) only allows exchange between
same-spin electrons, so \f$\alpha\f$ electrons only exchange with \f$\alpha\f$ electrons.

### Total energy

The total electronic energy is (Szabo & Ostlund eq. 2.266):

\f[
E_\text{elec} = \frac{1}{2}\left[
  \text{Tr}\!\left((H^\text{core} + F^\alpha) P^\alpha\right)
  + \text{Tr}\!\left((H^\text{core} + F^\beta) P^\beta\right)
\right]
\f]

where \f$\text{Tr}(A \cdot B) = \sum_{\mu\nu} A_{\mu\nu} B_{\mu\nu}\f$ is the Frobenius
inner product. The factor of \f$\tfrac{1}{2}\f$ avoids double-counting. The total energy
adds the nuclear repulsion constant:

\f[
E_\text{total} = E_\text{elec} + V_\text{nuc}
\f]

---

## 5. The SCF algorithm step by step

The UHF equations are coupled: \f$F^\alpha\f$ depends on \f$P^\beta\f$ through the total density
\f$P\f$, and \f$F^\beta\f$ depends on \f$P^\alpha\f$. The self-consistent field loop resolves this
iteratively.

**Step 0 — Orthogonalization.** Compute \f$X = S^{-1/2}\f$, the inverse square root of the
overlap matrix. This transforms the non-orthogonal AO basis into an orthonormal working
basis where the eigenvalue problem becomes standard (not generalized). \f$X\f$ is computed
once and reused for both spin channels at every iteration.

**Step 1 — Initial guess.** Diagonalize \f$H^\text{core}\f$ (ignoring electron-electron
repulsion) to get starting orbitals for both channels. For open-shell systems, the
\f$\alpha\f$ HOMO and LUMO are rotated by \f$20°\f$ to break the \f$\alpha/\beta\f$ symmetry;
without this rotation, both channels converge to the same restricted solution, missing
the lower-energy unrestricted minimum.

**Step 2 — Build the Fock matrices.** Given current \f$P^\alpha\f$ and \f$P^\beta\f$, compute
\f$G^\alpha\f$ and \f$G^\beta\f$ from the ERI tensor, then assemble \f$F^\alpha\f$ and \f$F^\beta\f$.

**Step 3 — Diagonalize each channel.** Transform each Fock matrix:

\f[
F'^\sigma = X^T F^\sigma X
\f]

then solve the standard eigenvalue problem \f$F'^\sigma \mathbf{v} = \varepsilon \mathbf{v}\f$
to get orbital energies \f$\varepsilon^\sigma\f$ and transformed coefficients \f$C'^\sigma\f$.
Back-transform: \f$C^\sigma = X C'^\sigma\f$. The same \f$X\f$ is used for both channels because
\f$X\f$ depends only on the basis set, not the spin state.

**Step 4 — Build new density matrices.** Fill the \f$N_\alpha\f$ (or \f$N_\beta\f$) lowest orbitals:

\f[
P^\alpha_\text{new} = C^\alpha_\text{occ} \left(C^\alpha_\text{occ}\right)^T, \qquad
  P^\beta_\text{new}  = C^\beta_\text{occ}  \left(C^\beta_\text{occ}\right)^T
\f]

**Step 5 — Linear mixing (damping).** To prevent oscillation, the new density is mixed
with the old:

\f[
P^\alpha_\text{next} = 0.5\, P^\alpha_\text{new} + 0.5\, P^\alpha_\text{old}
\f]
\f[
P^\beta_\text{next}  = 0.5\, P^\beta_\text{new}  + 0.5\, P^\beta_\text{old}
\f]

This conservative 50% mixing doubles the number of iterations compared to a direct
update but avoids oscillation when starting from an asymmetric guess far from
convergence. For small STO-3G bases (1–5 functions) it reliably converges to \f$10^{-9}\f$
within 60–150 iterations.

**Step 6 — Check convergence.** Measure the density change:

\f[
\text{residual} = \max\!\left(
  \|P^\alpha_\text{next} - P^\alpha_\text{old}\|_F,\;
  \|P^\beta_\text{next}  - P^\beta_\text{old}\|_F
\right)
\f]

If residual \f$<\f$ tolerance (default \f$10^{-8}\f$), the calculation is converged. Otherwise
return to Step 2 with the new densities.

---

## 6. Spin contamination

### Pure spin states

In exact quantum mechanics, a spin state is characterized by total spin \f$S\f$ and its
\f$z\f$-projection \f$M_S = (N_\alpha - N_\beta)/2\f$. The expected value of the spin-squared
operator \f$\hat{S}^2\f$ is exactly \f$S(S+1)\f$:

| State | \f$N_\alpha\f$ | \f$N_\beta\f$ | \f$S\f$ | \f$\langle\hat{S}^2\rangle_\text{exact}\f$ |
|-------|-----------|----------|-----|----------------------------------------|
| Singlet | 1 | 1 | 0 | 0.00 |
| Doublet | 1 | 0 | 0.5 | 0.75 |
| Doublet | 2 | 1 | 0.5 | 0.75 |
| Triplet | 2 | 0 | 1.0 | 2.00 |

### UHF spin contamination

Because UHF uses a single Slater determinant (not a full linear combination of
determinants as in the exact theory), the UHF wavefunction is generally **not** a pure
spin eigenstate. It is a mixture of the intended spin state with higher spin states.

The UHF value of \f$\langle\hat{S}^2\rangle\f$ is (Szabo & Ostlund eq. 2.273):

\f[
\langle\hat{S}^2\rangle_\text{UHF} = S(S+1) + N_\beta - \sum_{i,j} |O_{ij}|^2
\f]

where \f$O_{ij}\f$ is the overlap between occupied \f$\alpha\f$ MO \f$i\f$ and occupied \f$\beta\f$ MO \f$j\f$:

\f[
O_{ij} = \sum_{\mu\nu} C^\alpha_{\mu i}\, S_{\mu\nu}\, C^\beta_{\nu j}
         \equiv \left(\mathbf{C}^\alpha_\text{occ}\right)^T \mathbf{S}\, \mathbf{C}^\beta_\text{occ}
\f]

The **spin contamination** is defined as the deviation from the pure-state value:

\f[
\Delta\langle\hat{S}^2\rangle
  = \langle\hat{S}^2\rangle_\text{UHF} - S(S+1)
  = N_\beta - \text{Tr}(O O^T)
\f]

A spin contamination of 0 means the wavefunction is a pure spin state. A contamination
of 1.0 means the singlet-aimed wavefunction has equal contributions from the singlet and
triplet. Values below 0.1 are generally considered acceptable; larger values signal that
the result should be interpreted with caution or that spin-projection techniques are
needed.

In the H\f$_2\f$ dissociation example: near equilibrium, contamination is \f$\approx 0\f$. Past
the Coulson-Fischer point (\f$\sim\!1.2\f$ Å in STO-3G), it rises toward 1.0 as the \f$\alpha\f$
and \f$\beta\f$ electrons localize on opposite atoms. This is not a bug — it is the correct
physics of a broken-symmetry solution. The energy is lower than RHF because UHF uses a
more flexible wavefunction, even though that wavefunction is not a pure spin eigenstate.

---

## 7. Atomiverse API walkthrough

The UHF solver lives in the CPU backend (`ATV::cpu::scf`). It takes pre-computed
molecular integrals and returns a `UHFResult` struct.

### Includes

```cpp
#include "atomiverse/cpu/frame.h"          // cpu::Frame, Mat3
#include "atomiverse/cpu/scf/integrals.h"  // IntegralPack, build_integrals_sto3g
#include "atomiverse/cpu/scf/uhf_solve.h"  // UHFResult, uhf_solve
```

### Setting up the geometry

Geometry is stored in `cpu::Frame`. Each atom needs a 3D position in **Angstrom** and
a `type_id` equal to its atomic number (H = 1, He = 2, Li = 3):

```cpp
using namespace ATV;

// Hydrogen atom at the origin
cpu::Frame h_atom(1);
h_atom.x[0] = 0.0;  h_atom.y[0] = 0.0;  h_atom.z[0] = 0.0;
h_atom.type_id[0] = 1;          // Z = 1
h_atom.box = Mat3::cubic(100.0); // large box = isolated molecule

// H2 at equilibrium bond length
cpu::Frame h2(2);
h2.x[0] = 0.0;   h2.y[0] = 0.0;  h2.z[0] = 0.0;
h2.x[1] = 0.74;  h2.y[1] = 0.0;  h2.z[1] = 0.0;  // 0.74 Angstrom
h2.type_id[0] = 1;  h2.type_id[1] = 1;
h2.box = Mat3::cubic(100.0);
```

### Computing the integrals

```cpp
ATV::cpu::scf::IntegralPack ints = cpu::scf::build_integrals_sto3g(h2);
// ints.S       -- overlap matrix       [n_basis x n_basis]
// ints.Hcore   -- core Hamiltonian     [n_basis x n_basis]
// ints.ERI     -- two-electron integrals [n_basis^4]
// ints.V_nuc   -- nuclear repulsion energy [Hartree]
// ints.n_basis -- number of basis functions
```

For STO-3G: H contributes 1 basis function, He 1, Li 5. H\f$_2\f$ gives `n_basis = 2`.

### Running UHF

```cpp
int n_alpha = 1;    // spin-up electrons
int n_beta  = 1;    // spin-down electrons
Real tol    = 1e-8; // density convergence threshold (Frobenius norm)
int max_iter = 200;

ATV::cpu::scf::UHFResult result =
    cpu::scf::uhf_solve(ints, n_alpha, n_beta, tol, max_iter);
```

**Choosing \f$N_\alpha\f$ and \f$N_\beta\f$:** Total electrons \f$= N_\alpha + N_\beta = \sum Z\f$.
Spin \f$S = (N_\alpha - N_\beta)/2\f$.

| System | \f$N_\text{elec}\f$ | Spin state | \f$N_\alpha\f$ | \f$N_\beta\f$ |
|--------|----------------|------------|------------|----------|
| H atom | 1 | doublet (\f$S = 1/2\f$) | 1 | 0 |
| He atom | 2 | singlet (\f$S = 0\f$) | 1 | 1 |
| H\f$_2\f$ @ eq | 2 | singlet (\f$S = 0\f$) | 1 | 1 |
| Li atom | 3 | doublet (\f$S = 1/2\f$) | 2 | 1 |

### Reading the result

```cpp
if (!result.converged)
    printf("SCF did not converge in %d iterations\n", result.iterations);

printf("E_total:            %.8f Hartree\n", (double)result.energy);
printf("E_electronic:       %.8f Hartree\n", (double)result.electronic_energy);
printf("Iterations:         %d\n",            result.iterations);
printf("Residual:           %.2e\n",          (double)result.residual);
printf("Spin contamination: %.4f\n",          (double)result.spin_contamination);

// Orbital energies: one per basis function, sorted ascending
for (int i = 0; i < (int)result.eigenvalues_alpha.size(); ++i)
    printf("  alpha orbital %d: %.6f Ha\n", i, (double)result.eigenvalues_alpha(i));

// Electron count check: Tr(P^sigma * S) = N_sigma
Real n_alpha_check = xt::sum(result.P_alpha * ints.S)();
Real n_beta_check  = xt::sum(result.P_beta  * ints.S)();
```

`result.energy` is the total energy including nuclear repulsion. `result.electronic_energy`
is the electronic part only; their difference equals `ints.V_nuc`. For a single atom,
\f$V_\text{nuc} = 0\f$ so both fields are identical.

---

## 8. Key code paths annotated

The implementation is in `src/cpu/scf/uhf_solve.cpp`.

### `build_uhf_G` — constructing the two-electron part

```cpp
// G^alpha_ij = sum_kl  P_total_kl*(ij|kl) - P_alpha_kl*(ik|jl)
// G^beta_ij  = sum_kl  P_total_kl*(ij|kl) - P_beta_kl*(ik|jl)
for (i, j):
    for (k, l):
        coulomb = P_total(k,l) * ERI(i,j,k,l);   // (ij|kl): same for both channels
        ga += coulomb - P_alpha(k,l) * ERI(i,k,j,l);  // (ik|jl): alpha exchange only
        gb += coulomb - P_beta(k,l)  * ERI(i,k,j,l);  // (ik|jl): beta exchange only
```

**Why two different ERI index orderings?** \f$(\mu\nu|\lambda\sigma)\f$ is the Coulomb
integral: electron 1 in AO pair \f$(\mu,\nu)\f$ repels electron 2 in pair \f$(\lambda,\sigma)\f$.
This is spin-blind. The exchange \f$(\mu\lambda|\nu\sigma)\f$ arises from the antisymmetry
of the wavefunction (Pauli principle) and only acts between same-spin electrons. Total
density \f$P\f$ enters Coulomb; spin density \f$P^\sigma\f$ enters exchange.

### `ortho_diag` — solving the eigenvalue problem

```cpp
F_prime = X^T * F * X;              // transform to orthogonal basis
[eig, C_prime] = eigh(F_prime);     // standard symmetric eigenproblem
C = X * C_prime;                    // back-transform to AO basis
```

Without \f$X\f$, the problem is \f$F\mathbf{C} = S\mathbf{C}\boldsymbol{\varepsilon}\f$ — a
generalized eigenproblem. The transformation \f$X = S^{-1/2}\f$ converts it to standard
form \f$F'\mathbf{C}' = \mathbf{C}'\boldsymbol{\varepsilon}\f$. The same \f$X\f$ is used for
both \f$F^\alpha\f$ and \f$F^\beta\f$.

### `compute_spin_contamination` — measuring spin purity

The contamination is \f$N_\beta - \text{Tr}(OO^T)\f$ where
\f$O = (C^\alpha_\text{occ})^T \mathbf{S}\, C^\beta_\text{occ}\f$:

```cpp
SCb = S_mat * Cb;          // S * C^beta_occ
O   = Ca^T * SCb;          // O_ij = <phi^alpha_i | phi^beta_j>
return n_beta - sum(O * O); // N_beta - Tr(O O^T)
```

When \f$\alpha\f$ and \f$\beta\f$ orbitals are identical (closed-shell), \f$O\f$ is the identity and
\f$\text{Tr}(OO^T) = N_\beta\f$, giving contamination = 0. When they fully localize on
opposite atoms (dissociation limit), \f$O \approx 0\f$ and contamination \f$\to N_\beta = 1\f$.

### Energy evaluation — which density enters the trace

At each iteration, the energy is evaluated using the **freshly diagonalized** densities
\f$P_\text{new}\f$, not the linearly mixed \f$P_\text{next}\f$:

```cpp
E_elec = 0.5 * (sum((Hcore + F_alpha) * P_alpha_new) +
                sum((Hcore + F_beta)  * P_beta_new));
```

This ensures \f$F^\sigma\f$ and \f$P^\sigma\f$ in the energy expression are from the same SCF
state (Szabo & Ostlund eq. 2.266). The mixed \f$P_\text{next}\f$ is only used as the
starting density for the next iteration.

---

## 9. End-to-end example and plots

The example is `examples/hf.cpp`. It runs three calculations and, when built with
`ATV_USE_MATPLOT=ON`, saves two figures.

### Building

```bash
# Without plots (stdout only):
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON \
      -DATV_BUILD_EXAMPLES=ON -DATV_REAL=double
cmake --build build --target hf
./build/examples/hf

# With plots (requires gnuplot installed):
cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON \
      -DATV_BUILD_EXAMPLES=ON -DATV_USE_MATPLOT=ON -DATV_REAL=double
cmake --build build --target hf
./build/examples/hf
# Produces h2_pes.png and h2_spin_contam.png in the working directory.
```

### What the example computes

**Case 1 — H atom.** One electron, no electron-electron repulsion: \f$G^\alpha = 0\f$, so
\f$F^\alpha = H^\text{core}\f$. The total energy equals the \f$1s\f$ orbital energy: \f$-0.4666\f$
Hartree (STO-3G). This case is unreachable by RHF.

**Case 2 — H\f$_2\f$ at equilibrium (0.74 Å).** Both methods give essentially the same
energy. This verifies that UHF degenerates to RHF for well-behaved closed-shell
molecules at equilibrium.

**Case 3 — H\f$_2\f$ PES scan.** Bond length scanned from 0.5 to 4.0 Å in 0.25 Å steps.
This is the key validation: UHF must dissociate correctly while RHF must not.

---

## 10. Interpreting the output

### Figure 1: `h2_pes.png` — Potential energy surface

The figure shows \f$E_\text{RHF}\f$ (gray), \f$E_\text{UHF}\f$ (blue), and the exact
dissociation limit \f$2E_\text{H} = -0.9332\f$ Ha (red dashed) as a function of bond length.

**What to look for:**

- **Near 0.74 Å** (equilibrium): the RHF and UHF curves overlap completely. This
  confirms the UHF implementation degenerates to RHF for a closed-shell singlet.

- **Near 1.2 Å** (the *Coulson-Fischer point*): the UHF curve dips below RHF. This is
  where \f$\alpha\f$ and \f$\beta\f$ electrons first become energetically favorable to localize
  on separate atoms.

- **At large \f$R\f$ (> 2.5 Å)**: UHF asymptotes to the dashed reference at \f$-0.9332\f$ Ha.
  RHF plateaus near \f$-0.71\f$ Ha. The gap (\f$\sim\!0.22\f$ Ha \f$=\f$ 140 kcal/mol) is the
  RHF dissociation error.

If the UHF and RHF curves are identical everywhere, UHF is not breaking symmetry —
check that \f$N_\alpha \neq N_\beta\f$ is set, or that the HOMO-LUMO mixing is active.

### Figure 2: `h2_spin_contam.png` — Spin contamination

The figure shows \f$\Delta\langle\hat{S}^2\rangle = \langle\hat{S}^2\rangle_\text{UHF} - S(S+1)\f$
as a function of bond length, with a dashed gray line at 0 (pure singlet) and a dashed
red line at 1 (maximum contamination for a singlet-aimed wavefunction).

**What to look for:**

- **0.5–1.1 Å**: contamination \f$\approx 0\f$. The UHF solution is essentially a pure
  singlet, identical to RHF.

- **1.2–1.5 Å**: contamination rises rapidly. This is the Coulson-Fischer region where
  the broken-symmetry solution first becomes lower in energy.

- **\f$> 2.5\f$ Å**: contamination approaches 1.0. The \f$\alpha\f$ and \f$\beta\f$ electrons have
  fully separated onto different atoms. A value of 1.0 exactly corresponds to a 50/50
  mixture of singlet and triplet character — the expected behavior of an unrestricted
  broken-symmetry solution in the dissociation limit.

**Caution:** Spin contamination does not mean the energy is wrong. The lower UHF energy
is real and is guaranteed by the variational principle. Spin contamination matters for
properties that are spin-sensitive (hyperfine coupling constants, g-tensors). For
energies and geometries, the contamination values seen here are generally acceptable.

---

## 11. Molecular visualization with MolRec

After a calculation converges it is useful to inspect the per-atom charge
distribution and export the geometry to a visualization tool such as OVITO
or Avogadro 2.  Atomiverse provides three cooperating classes for this:

| Component | Header | Role |
|-----------|--------|------|
| `cpu::MolRec` | `cpu/molrec.h` | Trajectory container — stores frames, fields, scalars |
| `cpu::mulliken_charge` / `cpu::mulliken_spin` | `cpu/observer/mulliken.h` | Compute per-atom Mulliken populations |
| `cpu::XYZ` | `cpu/recorder/xyz.h` | Write plain XYZ or extended-XYZ (extxyz) files |

### 11.1 What is Mulliken population analysis?

After the SCF converges, the density matrix \f$P_{\mu\nu}\f$ tells us how the electrons
are distributed across the basis functions.  Each basis function \f$\mu\f$ is centred on
a specific atom \f$A\f$.  The **Mulliken gross population** of atom \f$A\f$ is

\f[

N_A = \sum_{\mu \in A} (PS)_{\mu\mu}

\f]

where the sum runs over the basis functions centred on \f$A\f$ and \f$S\f$ is the overlap matrix.
The **Mulliken net charge** is then

\f[

q_A = Z_A - N_A

\f]

where \f$Z_A\f$ is the nuclear charge.  For a neutral isolated molecule \f$\sum_A q_A = 0\f$.

For UHF, we compute this separately for each spin channel.  The **spin density** at atom
\f$A\f$ is

\f[

\sigma_A = \sum_{\mu \in A} \bigl((P^\alpha - P^\beta) S\bigr)_{\mu\mu}
= N_A^\alpha - N_A^\beta

\f]

A positive \f$\sigma_A\f$ means the atom has more \f$\alpha\f$ (up) electron density.  At H\f$_2\f$
dissociation, one H carries \f$\sigma \approx +1\f$ and the other \f$\sigma \approx -1\f$, showing
that the \f$\alpha\f$ and \f$\beta\f$ electrons have localized on separate nuclei.

### 11.2 Building a MolRec record

`MolRec` uses a builder pattern: call `set_geometry`, attach fields and scalars, then
call `commit_frame` to push the snapshot into the internal buffer.

```cpp
#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cpu/recorder/xyz.h"

// After uhf_solve(ints, 1, 1, tol, 200) returns result:
xt::xtensor<Real,2> P_total = result.P_alpha + result.P_beta;

cpu::MolRec rec;
rec.set_geometry(frame);                         // positions + element symbols from cpu::Frame
rec.add_string("method", "UHF");
rec.add_scalar("energy",             result.energy);
rec.add_scalar("iterations",         Real(result.iterations));
rec.add_scalar("n_alpha",            Real(result.n_alpha));
rec.add_scalar("n_beta",             Real(result.n_beta));
rec.add_scalar("spin_contamination", result.spin_contamination);
rec.add_field("charge",
    cpu::mulliken_charge(frame, P_total, ints.S));          // q_A = Z_A - N_A
rec.add_field("spin_density",
    cpu::mulliken_spin(frame, result.P_alpha, result.P_beta, ints.S));  // N^alpha_A - N^beta_A
rec.commit_frame();

cpu::print_summary(rec);                         // stdout table
cpu::XYZ("h2_uhf.xyz").write_ext(rec);           // extxyz — opens in OVITO / Avogadro 2
```

For RHF, pass `rhf.P` as both the density and set spin density to zero:

```cpp
rec.add_field("charge",
    cpu::mulliken_charge(frame, rhf_result.P, ints.S));
rec.add_field("spin_density",
    std::vector<Real>(frame.n_atoms(), Real(0)));    // RHF: no spin polarization
```

### 11.3 Trajectory output for PES scans

`MolRec` accumulates frames.  For a PES scan, stream each frame directly to disk with
`append_ext` so memory does not grow with the scan:

```cpp
cpu::MolRec scan;
cpu::XYZ xyz("h2_scan_uhf.xyz");

for (double r = 0.5; r <= 4.0; r += 0.25) {
    auto frame = diatomic_h2(r);
    auto ints  = cpu::scf::build_integrals_sto3g(frame);
    auto res   = cpu::scf::uhf_solve(ints, 1, 1, tol, 400);

    xt::xtensor<Real,2> P_total = res.P_alpha + res.P_beta;
    scan.set_geometry(frame);
    scan.add_scalar("energy", res.energy);
    scan.add_scalar("spin_contamination", res.spin_contamination);
    scan.add_field("charge",       cpu::mulliken_charge(frame, P_total, ints.S));
    scan.add_field("spin_density", cpu::mulliken_spin(frame, res.P_alpha, res.P_beta, ints.S));
    scan.commit_frame();
    xyz.append_ext(scan);     // write this frame now; do not buffer all frames in memory
}
```

The resulting multi-frame `.xyz` file can be animated in OVITO.  The `charge` and
`spin_density` columns appear as per-atom properties that can be colour-mapped.

### 11.4 Plain XYZ vs extended XYZ

| Method | Comment line | Per-atom columns | Viewer support |
|--------|-------------|-----------------|----------------|
| `write(rec)` / `append(rec)` | `E=<energy>` | symbol, x, y, z | Universal |
| `write_ext(rec)` / `append_ext(rec)` | Full extxyz metadata | symbol, x, y, z, + all `add_field` columns | OVITO, Avogadro 2, ASE, VMD |

Use plain XYZ when you only need geometry (e.g. to feed into another program).
Use extxyz when you want to visualize per-atom properties or read the file back
with ASE (`ase.io.read`).

---

## 12. Common pitfalls

### Wrong electron count

**Symptom:** \f$\text{Tr}(P^\alpha S)\f$ reports 2.0 when you expected 1.0, or the energy
is much lower than the reference.

**Cause:** Passing the total electron count as `n_alpha` instead of the spin-up count.
For a 2-electron singlet, the correct call is `uhf_solve(ints, 1, 1, ...)` — one alpha
and one beta, not two alpha.

**Fix:** \f$N_\alpha + N_\beta = N_\text{elec}\f$, \f$N_\alpha - N_\beta = 2S\f$.

### No symmetry breaking for open-shell systems

**Symptom:** For an open-shell system, `spin_contamination = 0` and `energy` equals
the RHF value.

**Cause:** Both channels converged to the restricted (symmetric) solution. The 20°
HOMO-LUMO mixing is applied only when \f$N_\alpha > 0\f$, \f$N_\alpha < N_\text{basis}\f$, and
\f$N_\beta > 0\f$. If `n_beta = 0` (H atom), contamination is 0 by construction — that is
correct, not a bug.

**Fix:** Verify \f$N_\alpha > N_\beta\f$ is set correctly for open-shell systems. For Li
(doublet), use `n_alpha = 2, n_beta = 1`.

### Non-convergence near the Coulson-Fischer point

**Symptom:** `result.converged = false` for bond lengths around 1.0–1.5 Å.

**Cause:** Near the Coulson-Fischer point, the symmetric and broken-symmetry solutions
are close in energy. The initial 20° perturbation can oscillate between them.

**Fix:** Increase `max_iter` to 400–600 for PES scans. The example uses `max_iter_scan = 400`.

### Incorrect geometry units

**Symptom:** Energy far from the reference (H\f$_2\f$ equilibrium should be \f$\approx -1.117\f$ Ha).

**Cause:** Positions in Bohr instead of Angstrom, or vice versa.

**Fix:** All positions in `cpu::Frame` are in **Angstrom**. The H\f$_2\f$ equilibrium bond
length is 0.74 Å \f$= 1.397\f$ Bohr. Passing 1.397 as Angstrom places atoms too close,
giving an energy that is too negative.

### Build fails: `uhf_solve` not found

**Symptom:** Link error: `undefined symbol cpu::scf::uhf_solve`.

**Cause:** `uhf_solve.cpp` is compiled only when `ATV_BUILD_RHF=ON`.

**Fix:** `cmake -B build -DATV_BUILD_RHF=ON ...`

---

## References

1. A. Szabo & N. S. Ostlund, *Modern Quantum Chemistry*, Dover (1996).
   Eqs. 2.254–2.266 (UHF density, Fock, energy), eq. 2.273 (spin contamination).
2. H. B. Schlegel, *J. Chem. Phys.* **84**, 4530 (1986). Spin contamination formula.
3. W. J. Hehre, R. F. Stewart & J. A. Pople, *J. Chem. Phys.* **51**, 2657 (1969). STO-3G basis.
4. C. C. J. Roothaan, *Rev. Mod. Phys.* **23**, 69 (1951). Roothaan-Hall equations.
