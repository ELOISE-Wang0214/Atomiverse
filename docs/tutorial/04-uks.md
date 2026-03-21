# Tutorial 4: Unrestricted Kohn-Sham DFT (UKS)

This tutorial is for readers who know basic linear algebra but have **no background
in computational chemistry**.

> **Prerequisites:**
> Read [Tutorial 1: Toy SCF from First Principles](01-scf.md) first.
>
> **Recommended:**
> [Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md) and
> [Tutorial 3: Unrestricted Hartree-Fock (UHF)](03-uhf.md).

The goal of this page is to make three ideas fit into one mental picture:

- RHF: the simplest closed-shell Hartree-Fock.
- UHF: what changes when alpha and beta electrons are allowed to differ.
- UKS: what changes when we replace Hartree-Fock exchange with a density functional.

If you are new, focus on the SCF loop structure. The details of the integrals and
the exact form of the exchange-correlation functional can be learned later.

## What Problem Are We Solving?

Given fixed nuclei positions \f$\{\mathbf{R}_A\}\f$ and charges \f$\{Z_A\}\f$, we want the
lowest (ground-state) electronic energy. In practice we do not solve the exact
many-electron Schr\"odinger equation. Instead we solve a mean-field fixed point:

\f[

\text{guess a density} \;\to\; \text{build an effective one-electron matrix}
\;\to\; \text{diagonalize} \;\to\; \text{rebuild the density} \;\to\; \text{repeat}.

\f]

That feedback loop is SCF, and it is the shared computational skeleton of RHF, UHF,
and UKS.

## From Orbitals to Matrices: The One Translation Step

In all methods here, we expand orbitals in a finite basis \f$\{\chi_\mu\}\f$:

\f[

\varphi_i(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu i}\,\chi_\mu(\mathbf{r}),

\f]

so solving for orbitals becomes solving for the coefficient matrices \f$\mathbf{C}\f$.
The basic inputs are molecular integrals in that basis:

- Overlap \f$S_{\mu\nu} = \int \chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r})\,d\mathbf{r}\f$.
- Core Hamiltonian \f$H^{\text{core}}\f$ (kinetic + electron-nucleus attraction).
- Two-electron repulsion integrals (ERI) \f$(\mu\nu|\lambda\sigma)\f$.
- Nuclear repulsion energy \f$V_{\text{nuc}}\f$.

Atomiverse packages these into an `IntegralPack`.

## Quick Recap: RHF in One Page

RHF assumes electrons are paired: alpha and beta share the same spatial orbitals.
The Roothaan-Hall equation is a generalized eigenproblem:

\f[

\mathbf{F}\mathbf{C}=\mathbf{S}\mathbf{C}\boldsymbol{\varepsilon}.

\f]

For a closed-shell system with \f$N\f$ electrons, the RHF density matrix is

\f[

P_{\mu\nu} = 2\sum_{i=1}^{N/2} C_{\mu i} C_{\nu i}.

\f]

The Fock matrix is

\f[

F_{\mu\nu} = H^{\text{core}}_{\mu\nu} + G_{\mu\nu}(\mathbf{P}),

\f]

with the two-electron mean-field term

\f[

G_{\mu\nu}(\mathbf{P}) =
\sum_{\lambda\sigma}P_{\lambda\sigma}
\left[(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)\right].

\f]

At convergence, a common total-energy expression is

\f[

E_{\text{RHF}} =
\frac{1}{2}\mathrm{Tr}\left[\mathbf{P}\left(\mathbf{H}^{\text{core}}+\mathbf{F}\right)\right]
+ V_{\text{nuc}}.

\f]

The SCF loop is exactly the Toy SCF loop, except the "matrix update" is now a real
integral contraction, not \f$H_0 + \alpha P\f$.

## Quick Recap: UHF in One Page

UHF removes the pairing constraint by allowing separate orbitals for alpha and beta:

\f[

\varphi_i^\alpha(\mathbf{r})=\sum_\mu C_{\mu i}^\alpha\chi_\mu(\mathbf{r}),
\qquad
\varphi_i^\beta(\mathbf{r})=\sum_\mu C_{\mu i}^\beta\chi_\mu(\mathbf{r}).

\f]

The densities split:

\f[

\mathbf{P}^\alpha=\mathbf{C}^\alpha_{\text{occ}}\left(\mathbf{C}^\alpha_{\text{occ}}\right)^T,
\qquad
\mathbf{P}^\beta=\mathbf{C}^\beta_{\text{occ}}\left(\mathbf{C}^\beta_{\text{occ}}\right)^T,
\qquad
\mathbf{P}=\mathbf{P}^\alpha+\mathbf{P}^\beta.

\f]

The UHF Fock matrices are

\f[

\mathbf{F}^\alpha=\mathbf{H}^{\text{core}}+\mathbf{J}(\mathbf{P})-\mathbf{K}(\mathbf{P}^\alpha),
\qquad
\mathbf{F}^\beta=\mathbf{H}^{\text{core}}+\mathbf{J}(\mathbf{P})-\mathbf{K}(\mathbf{P}^\beta),

\f]

where the Coulomb term uses total density, while exchange is spin-specific:

\f[

J_{\mu\nu}=\sum_{\lambda\sigma}P_{\lambda\sigma}(\mu\nu|\lambda\sigma),
\qquad
K^\alpha_{\mu\nu}=\sum_{\lambda\sigma}P^\alpha_{\lambda\sigma}(\mu\lambda|\nu\sigma),
\qquad
K^\beta_{\mu\nu}=\sum_{\lambda\sigma}P^\beta_{\lambda\sigma}(\mu\lambda|\nu\sigma).

\f]

The algorithmic structure is the same as RHF, just duplicated per spin channel
and coupled through \f$\mathbf{J}(\mathbf{P})\f$.

## UKS: What Changes When We Switch to Kohn-Sham DFT?

Hartree-Fock has an explicit exchange matrix \f$\mathbf{K}\f$ that comes from the
antisymmetry of a single Slater determinant. It also misses most correlation energy.

Kohn-Sham DFT replaces "explicit exchange" with an exchange-correlation functional:
we compute a scalar energy \f$E_{\text{xc}}[\rho_\alpha,\rho_\beta]\f$ and its
corresponding potentials, and we never build \f$\mathbf{K}\f$ explicitly.

### Step 1: From AO density matrices to real-space densities

At a point \f$\mathbf{r}\f$, the spin densities are

\f[

\rho_\sigma(\mathbf{r})
=\sum_{\mu\nu}P^\sigma_{\mu\nu}\,\chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r}),
\qquad \sigma\in\{\alpha,\beta\}.

\f]

Total density and magnetization density are

\f[

\rho(\mathbf{r})=\rho_\alpha(\mathbf{r})+\rho_\beta(\mathbf{r}),
\qquad
m(\mathbf{r})=\rho_\alpha(\mathbf{r})-\rho_\beta(\mathbf{r}).

\f]

For GGA functionals you also need gradients. A practical AO expression is

\f[

\nabla\rho_\sigma(\mathbf{r})
=2\sum_{\mu\nu}P^\sigma_{\mu\nu}\,\left(\nabla\chi_\mu(\mathbf{r})\right)\chi_\nu(\mathbf{r}).

\f]

### Step 2: The UKS (Kohn-Sham) matrices

UKS solves two coupled eigenproblems, one per spin:

\f[

\mathbf{F}^\alpha=\mathbf{H}^{\text{core}}+\mathbf{J}(\mathbf{P})+\mathbf{V}_{\text{xc}}^\alpha,
\qquad
\mathbf{F}^\beta=\mathbf{H}^{\text{core}}+\mathbf{J}(\mathbf{P})+\mathbf{V}_{\text{xc}}^\beta.

\f]

Compared to UHF, the key difference is:

\f[

\text{UHF has } -\mathbf{K}(\mathbf{P}^\sigma), \qquad \text{UKS does not.}

\f]

All exchange and correlation are absorbed into \f$\mathbf{V}_{\text{xc}}^\sigma\f$.

### Step 3: The UKS total energy, and the easiest mistake

A common energy decomposition is

\f[

E_{\text{KS}}
=\sum_{\sigma\in\{\alpha,\beta\}}
\mathrm{Tr}\left[\mathbf{H}^{\text{core}}\mathbf{P}^\sigma\right]
+\frac{1}{2}\mathrm{Tr}\left[\mathbf{J}(\mathbf{P})\mathbf{P}\right]
+E_{\text{xc}}[\rho_\alpha,\rho_\beta]
+V_{\text{nuc}}.

\f]

The easiest implementation mistake is to assume
\f$E_{\text{xc}} = \mathrm{Tr}(\mathbf{V}_{\text{xc}}\mathbf{P})\f$.
In general they are not equal. In real-space language:

\f[

E_{\text{xc}}=\int \rho(\mathbf{r})\,\varepsilon_{\text{xc}}(\mathbf{r})\,d\mathbf{r},
\qquad
\mathrm{Tr}(\mathbf{V}_{\text{xc}}\mathbf{P})
=\int \rho(\mathbf{r})\,v_{\text{xc}}(\mathbf{r})\,d\mathbf{r},

\f]

and typically

\f[

v_{\text{xc}}(\mathbf{r})
=\varepsilon_{\text{xc}}(\mathbf{r})
+\rho(\mathbf{r})\frac{\partial\varepsilon_{\text{xc}}}{\partial\rho}(\mathbf{r})
\neq
\varepsilon_{\text{xc}}(\mathbf{r}).

\f]

So \f$E_{\text{xc}}\f$ must be accumulated separately from the grid energy density,
not reconstructed from matrix traces.

## Grid Integration: Turning \f$\int d\mathbf{r}\f$ into a Sum

DFT needs real-space integrals. Numerically, we replace

\f[

\int f(\mathbf{r})\,d\mathbf{r}
\approx
\sum_g w_g f(\mathbf{r}_g),

\f]

where \f$\mathbf{r}_g\f$ are grid points and \f$w_g\f$ are weights. Molecular DFT codes
typically use atom-centered grids and a partition-of-unity (e.g. Becke weights) so the
grid adapts naturally to molecular geometry.

In an LDA/LSDA setting, an intuitive view of the matrix assembly is:

\f[

\left(V_{\text{xc}}^\sigma\right)_{\mu\nu}
\approx
\sum_g w_g\,v_{\text{xc}}^\sigma(\mathbf{r}_g)\,\chi_\mu(\mathbf{r}_g)\chi_\nu(\mathbf{r}_g),

\f]

while the XC energy is

\f[

E_{\text{xc}}
\approx
\sum_g w_g\,\rho(\mathbf{r}_g)\,\varepsilon_{\text{xc}}(\mathbf{r}_g).

\f]

GGA adds gradient-dependent terms, but the computational pattern remains the same:
evaluate local quantities at points, then sum with weights.

## A Concrete XC Example: LSDA Exchange and the Idea of PBE

In code, the XC interface at a point looks like:

- inputs: \f$\rho_\alpha, \rho_\beta\f$ (and for GGA, gradient invariants),
- outputs: \f$\varepsilon_{\text{xc}}\f$ and \f$v_{\text{xc}}^\alpha, v_{\text{xc}}^\beta\f$
  (and for GGA, derivatives w.r.t. gradient invariants).

LSDA exchange (Dirac-Slater) is a classic example:

\f[

E_x^{\text{LSDA}}
=-C_x\int\left[\rho_\alpha(\mathbf{r})^{4/3}+\rho_\beta(\mathbf{r})^{4/3}\right]d\mathbf{r},
\qquad
C_x=\frac{3}{4}\left(\frac{3}{\pi}\right)^{1/3}.

\f]

The corresponding exchange potentials are simple:

\f[

v_x^\sigma(\mathbf{r})
=\frac{\delta E_x}{\delta\rho_\sigma(\mathbf{r})}
=-\left(\frac{3}{\pi}\right)^{1/3}\rho_\sigma(\mathbf{r})^{1/3}.

\f]

PBE is a GGA functional: it uses density gradients to correct LDA. One way to see this
is that PBE multiplies the LDA exchange energy density by an enhancement factor
\f$F_x(s)\f$ depending on a reduced gradient \f$s\f$:

\f[

F_x(s)=1+\kappa-\frac{\kappa}{1+\mu s^2/\kappa},
\qquad
\kappa=0.804,\ \mu=0.21951.

\f]

The important takeaway is not the constants. It is that GGA requires gradients, which
means your grid evaluation must provide \f$\nabla\rho_\sigma\f$ and your functional
must return extra derivative information.

## The UKS SCF Algorithm Step by Step

The iteration pattern is almost the same as UHF. The main new work is the grid loop.

1. Build the molecular grid once from the nuclear coordinates.
2. Build the orthogonalization matrix \f$\mathbf{X}=\mathbf{S}^{-1/2}\f$ once.
3. Build an initial guess for \f$\mathbf{P}^{\alpha(0)}\f$ and \f$\mathbf{P}^{\beta(0)}\f$.
4. For iteration \f$k=0,1,2,\ldots\f$:

\f[

\mathbf{P}^{(k)}=\mathbf{P}^{\alpha(k)}+\mathbf{P}^{\beta(k)},
\qquad
\mathbf{J}^{(k)}=\mathbf{J}(\mathbf{P}^{(k)}).

\f]

Compute \f$E_{\text{xc}}^{(k)}\f$, \f$\mathbf{V}_{\text{xc}}^{\alpha(k)}\f$,
\f$\mathbf{V}_{\text{xc}}^{\beta(k)}\f$ from a grid quadrature.

Build spin-channel matrices and diagonalize:

\f[

\mathbf{F}^{\sigma(k)}=\mathbf{H}^{\text{core}}+\mathbf{J}^{(k)}+\mathbf{V}_{\text{xc}}^{\sigma(k)},
\qquad
\mathbf{F}^{\prime\sigma(k)}=\mathbf{X}^T\mathbf{F}^{\sigma(k)}\mathbf{X}.

\f]

Rebuild densities from occupied orbitals, check the residual
\f$\|\mathbf{P}^{\sigma(k+1)}-\mathbf{P}^{\sigma(k)}\|_F\f$, and stop when it is below a
tolerance.

## Common Numerical Pitfalls

- **Density near zero.** Far from nuclei \f$\rho(\mathbf{r})\to 0\f$ and many XC formulas can
  produce NaN/Inf. In practice you skip XC evaluation when \f$\rho\f$ is below a small
  cutoff.
- **The \f$E_{\text{xc}}\f$ mistake.** Always accumulate \f$E_{\text{xc}}\f$ from the grid
  energy density. Do not replace it by \f$\mathrm{Tr}(\mathbf{V}_{\text{xc}}\mathbf{P})\f$.
- **Grid too coarse.** After convergence, a cheap diagnostic is the electron-count check
  \f$\int \rho_\alpha d\mathbf{r}\f$ and \f$\int \rho_\beta d\mathbf{r}\f$. Significant
  deviation means the grid is not accurate enough.
