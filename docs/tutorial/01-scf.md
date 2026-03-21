# Tutorial 1: Toy SCF from First Principles

This tutorial is for readers who know basic linear algebra but have **no
background in computational chemistry**.

The goal is not to learn chemistry jargon first. The goal is to understand one
simple idea:

> sometimes the matrix you want to solve depends on the solution you got from
> the previous solve.

That feedback loop is the core of a **self-consistent field** procedure, or
SCF.

This page uses the smallest possible model in Atomiverse:

- a dense symmetric matrix \f$H_0\f$,
- a density matrix \f$P\f$,
- a coupling constant \f$\alpha\f$,
- repeated diagonalization until the density stops changing.

Once this toy model makes sense, the real Hartree-Fock tutorial becomes much
easier to read.

## Why Start Here?

Real electronic-structure methods contain three layers of complexity at once:

1. the self-consistent iteration,
2. the quantum-chemistry objects such as overlap and Fock matrices,
3. the implementation details of CPU and GPU backends.

The toy model removes layer 2 on purpose. That lets you focus on the iteration
pattern that Atomiverse is built around.

## The One Idea Behind SCF

Suppose you want a matrix \f$H\f$, but \f$H\f$ depends on another matrix \f$P\f$, and
\f$P\f$ itself is built from the eigenvectors of \f$H\f$. Then you cannot solve the
problem in one shot. You have to iterate:

\f[
\text{guess } P^{(0)}
\f]

\f[
H^{(k)} = H_0 + \alpha P^{(k)}
\f]

\f[
H^{(k)} C^{(k)} = C^{(k)} \varepsilon^{(k)}
\f]

\f[
P^{(k+1)} = C_{\mathrm{occ}}^{(k)} \left(C_{\mathrm{occ}}^{(k)}\right)^T
\f]

and stop when

\f[
\left\|P^{(k+1)} - P^{(k)}\right\|_F < \mathrm{tol}.
\f]

The meaning of each symbol here is simple:

| Symbol | Meaning |
|--------|---------|
| \f$H_0\f$ | A fixed starting matrix |
| \f$\alpha\f$ | How strongly the current density feeds back into the next Hamiltonian |
| \f$H\f$ | The matrix we diagonalize at each iteration |
| \f$C\f$ | Eigenvectors of \f$H\f$ |
| \f$\varepsilon\f$ | Eigenvalues of \f$H\f$ |
| \f$C_{\mathrm{occ}}\f$ | The columns of \f$C\f$ corresponding to occupied orbitals |
| \f$P\f$ | Density matrix rebuilt from those occupied vectors |

The toy model is not intended to be physically realistic. It is intended to
teach the computational pattern cleanly.

## What the Density Matrix Means in This Tutorial

In real quantum chemistry, the density matrix tells you how electrons occupy
orbitals and how that occupancy influences the effective Hamiltonian.

Here, you can treat \f$P\f$ more abstractly:

- it is a symmetric matrix built from the currently occupied eigenvectors,
- it summarizes the current state of the solution,
- it feeds back into the next iteration through \f$H = H_0 + \alpha P\f$.

That feedback is the entire point.

## What the Code Actually Does

This tutorial corresponds to two concrete programs in the repository:

| File | What it is | Why it exists |
|------|------------|---------------|
| `examples/toy-scf.cpp` | A CPU-only, explicit implementation with xtensor-blas | Best place to see the math step by step |
| `examples/scf.cu` | The same toy SCF workflow expressed through Atomiverse's `scf::SCF` driver | Best place to see how the same iteration is scheduled in the framework |

There is also a library CPU reference used for parity checks:

| File | Role |
|------|------|
| `include/atomiverse/cpu/scf/solve.h` + `src/cpu/scf/solve.cpp` | CPU library implementation of the same toy SCF algorithm |

So the point is not "three unrelated levels." The point is:

- one mathematical loop,
- one explicit CPU example,
- one CPU library implementation,
- one Driver-based GPU implementation.

## Step 1: Build the Fixed Matrix \f$H_0\f$

In [`examples/toy-scf.cpp`](../../examples/toy-scf.cpp), the code builds a
deterministic symmetric matrix:

```cpp
xt::xtensor<double, 2> H0 = xt::zeros<double>({n, n});
for (int i = 0; i < n; ++i) {
    H0(i, i) = double(n + 1 - i);
    for (int j = i + 1; j < n; ++j) {
        double val = 1.0 / (1.0 + (j - i));
        H0(i, j) = val;
        H0(j, i) = val;
    }
}
```

This construction is artificial, but useful:

- the diagonal entries dominate, so the matrix is well behaved,
- the off-diagonal entries are nonzero, so eigenvectors are mixed rather than
  trivial,
- the matrix is reproducible, so CPU and GPU runs can be compared.

No chemistry enters yet. This is just a stable test problem.

## Step 2: Start with a Guess for \f$P\f$

The code begins with

```cpp
xt::xtensor<double, 2> P     = xt::zeros<double>({n, n});
xt::xtensor<double, 2> P_old = xt::zeros<double>({n, n});
```

Starting from zero is the simplest possible initial guess. In real SCF codes,
better guesses are often used to accelerate convergence, but zero is enough to
explain the algorithm.

## Step 3: Build the Current Hamiltonian

Each iteration begins with

```cpp
xt::xtensor<double, 2> H = H0 + alpha * P;
```

This is the self-consistent step. The Hamiltonian is no longer fixed; it is
updated using the density from the previous iteration.

Interpretation:

- if \f$\alpha = 0\f$, the problem is just a one-shot eigensolve of \f$H_0\f$,
- if \f$\alpha > 0\f$, the current occupied subspace changes the matrix you will
  diagonalize next.

That is the feedback loop.

## Step 4: Diagonalize the Matrix

The next line is

```cpp
auto [eig, C] = xt::linalg::eigh(H);
```

Mathematically, this solves

\f[
H C = C \varepsilon,
\f]

where:

- the columns of \f$C\f$ are eigenvectors,
- the entries of \f$\varepsilon\f$ are eigenvalues.

Because \f$H\f$ is symmetric, this is a standard dense symmetric eigensolve. On
 CPU, xtensor-blas calls LAPACK underneath. In the Atomiverse GPU path, the
same role is played by cuSOLVER.

## Step 5: Rebuild the Density Matrix

The occupied orbitals are taken from the first `n_occ` columns:

```cpp
P_old = P;
auto C_occ = xt::view(C, xt::all(), xt::range(0, n_occ));
P = xt::linalg::dot(C_occ, xt::transpose(C_occ));
```

This is the key update:

\f[
P = C_{\mathrm{occ}} C_{\mathrm{occ}}^T.
\f]

Why this form?

- \f$C_{\mathrm{occ}}\f$ contains the vectors we currently regard as occupied,
- multiplying by its transpose builds the projector onto that occupied
  subspace,
- that projector becomes the new density matrix used in the next Hamiltonian.

So the density is not an arbitrary matrix. It is built from the current
solution itself.

## Step 6: Check Whether the Iteration Has Settled

The convergence check in the example is

```cpp
double residual = xt::norm_l2(P - P_old)();
if (iter > 0 && residual < tol) break;
```

The mathematical condition is

\f[
\left\|P^{(k+1)} - P^{(k)}\right\|_F < \mathrm{tol}.
\f]

The idea is straightforward:

- if the density changes a lot, the solution is still moving,
- if the density barely changes, the loop is self-consistent.

The first iteration is skipped because the initial `P_old = 0` guess is not yet
meaningful as a convergence baseline.

## What the Reported Energy Means

The example also prints

```cpp
double energy = xt::sum(H * P)();
```

which corresponds to

\f[
E = \mathrm{Tr}(H P).
\f]

This is a toy energy expression. It is useful because it gives a scalar number
that changes as the SCF loop converges, but it is not meant to be a chemically
accurate Hartree-Fock total energy.

That distinction matters:

- in the toy model, the energy is an easy diagnostic,
- in the RHF tutorial, the energy expression becomes physically meaningful.

## How This Maps onto Atomiverse

The GPU example [`examples/scf.cu`](../../examples/scf.cu) uses the same
algorithm, but instead of writing the loop manually, it registers components
with the SCF driver:

```cpp
scf::SCF solver;
solver.add(scf::Hamiltonian{n_occ, alpha});
solver.add(scf::EigenSolver{});
solver.add(scf::DensityBuilder{});
solver.add(scf::ConvergenceCheck{tol});
solver.on(scf::Phase::end_of_step).add(Thermo{5});
solver.build();
solver.run(frame, max_iter);
```

The phase sequence is fixed by the driver:

\f[
\texttt{build\_hamiltonian}
\rightarrow
\texttt{diagonalize}
\rightarrow
\texttt{build\_density}
\rightarrow
\texttt{convergence\_check}
\rightarrow
\texttt{end\_of\_step}.
\f]

The mapping from math to components is:

| Mathematical step | Atomiverse phase | Component |
|-------------------|------------------|-----------|
| \f$H = H_0 + \alpha P\f$ | `build_hamiltonian` | `scf::Hamiltonian` |
| diagonalize \f$H\f$ | `diagonalize` | `scf::EigenSolver` |
| \f$P = C_{\mathrm{occ}} C_{\mathrm{occ}}^T\f$ | `build_density` | `scf::DensityBuilder` |
| test \f$\|P - P_{\mathrm{old}}\|_F\f$ | `convergence_check` | `scf::ConvergenceCheck` |

This is the real architectural lesson of the tutorial:

- the CPU example makes the mathematics transparent,
- the Driver example shows how Atomiverse schedules the same mathematics.

## What to Look for When You Compare CPU and GPU Versions

When reading [`examples/toy-scf.cpp`](../../examples/toy-scf.cpp) and
[`examples/scf.cu`](../../examples/scf.cu), focus on what stays the same:

- the same fixed matrix \f$H_0\f$,
- the same update rule \f$H = H_0 + \alpha P\f$,
- the same density rebuild rule,
- the same convergence test.

What changes is only the implementation strategy:

- explicit xtensor code on CPU in the tutorial example,
- reusable components scheduled by `scf::SCF` on GPU in Atomiverse.

That separation is intentional.

## Before Moving On

If the toy model is clear, you are ready for the next tutorial. The main new
ideas there are:

1. \f$H\f$ is replaced by a physically meaningful Fock matrix,
2. the density is built from molecular orbitals of a real basis set,
3. the energy becomes a real RHF total energy for \f$\mathrm{H}_2\f$.

## Next Tutorial

Continue with [Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md).
