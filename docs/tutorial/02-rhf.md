# Tutorial 2: Restricted Hartree-Fock (RHF)

This tutorial assumes you have read
[Tutorial 1: Toy SCF from First Principles](01-scf.md).

The first tutorial explained the **iteration pattern** of SCF. This tutorial
adds the missing chemistry.

The guiding question is:

> how do we turn the toy feedback loop into a physically meaningful
> Restricted Hartree-Fock (RHF) calculation?

To keep everything as small as possible, we will use a single concrete example
throughout: the hydrogen molecule in the minimal STO-3G basis. The method is the
point; H2/STO-3G is only the smallest test case that makes the method fully real.

In this repository, the answer is split into two parts:

1. compute molecular integrals on the CPU,
2. solve the RHF fixed-point problem on either CPU or GPU.

That division matches the actual code in Atomiverse. This page follows the code
as it really exists, rather than describing an imaginary standalone textbook
implementation.

## What Problem Are We Solving?

We want the ground-state electronic energy of a molecule in a finite basis. In
this tutorial, that molecule is \f$\mathrm{H}_2\f$ and the basis is STO-3G.

At the Hartree-Fock level, we do not solve for the full many-electron
wavefunction directly. Instead, we solve for one-electron orbitals that are
self-consistent with the average field created by all electrons.

That is why this is still an SCF problem:

- the orbitals determine the density,
- the density determines the effective Hamiltonian,
- the effective Hamiltonian determines the next orbitals.

## Why a Minimal Example (H2 / STO-3G)?

\f$\mathrm{H}_2\f$ is the smallest nontrivial molecule:

- two nuclei,
- two electrons,
- one occupied spatial orbital in a closed-shell RHF picture.

STO-3G is a minimal basis:

- each hydrogen \f$1s\f$ orbital is approximated by a fixed combination of three
  Gaussian primitives,
- the basis is small enough to inspect by hand,
- the resulting matrices are only \f$2 \times 2\f$ for \f$\mathrm{H}_2\f$.

That makes it an ideal first physically meaningful example.

## The Key Objects You Need to Know

Before looking at code, it helps to know what the matrices mean.

| Symbol | Name | Role |
|--------|------|------|
| \f$S\f$ | Overlap matrix | Basis functions are not orthogonal, so this measures how much they overlap |
| \f$T\f$ | Kinetic matrix | One-electron kinetic energy contribution |
| \f$V\f$ | Nuclear attraction matrix | Attraction between electrons and nuclei |
| \f$H_{\mathrm{core}}\f$ | Core Hamiltonian | Defined as \f$T + V\f$ |
| \f$(\mu\nu|\lambda\sigma)\f$ | Two-electron repulsion integrals | Electron-electron interaction tensor |
| \f$F\f$ | Fock matrix | Effective one-electron matrix used in RHF |
| \f$P\f$ | Density matrix | Built from occupied molecular orbitals |
| \f$X\f$ | Orthogonalization matrix | Defined from \f$S^{-1/2}\f$ to convert the generalized eigenproblem into a standard one |

The idea is:

- \f$H_{\mathrm{core}}\f$ describes the one-electron part,
- the two-electron integrals describe electron-electron repulsion,
- together they define the Fock matrix \f$F\f$,
- the occupied eigenvectors of \f$F\f$ define the next density matrix.

## The RHF Equations in Plain Language

In a non-orthogonal atomic-orbital basis, RHF solves

\f[

F C = S C \varepsilon.

\f]

This differs from the toy tutorial in one important way: the basis functions
overlap, so the right-hand side is not just \f$C \varepsilon\f$.

To make the problem look like a standard symmetric eigensolve, we introduce
\f$X = S^{-1/2}\f$ and transform:

\f[

F' = X^T F X

\f]

\f[

F' C' = C' \varepsilon

\f]

\f[

C = X C'.

\f]

The density matrix for closed-shell RHF is then

\f[

P = 2 C_{\mathrm{occ}} C_{\mathrm{occ}}^T.

\f]

The factor of \f$2\f$ is because each occupied spatial orbital holds two
electrons, one with \f$\alpha\f$ spin and one with \f$\beta\f$ spin.

The Fock matrix itself is

\f[

F_{\mu\nu} = H_{\mathrm{core},\mu\nu} + G_{\mu\nu}(P),

\f]

with

\f[

G_{\mu\nu}(P) =
\sum_{\lambda\sigma}
P_{\lambda\sigma}
\left[
(\mu\nu|\lambda\sigma) - \frac{1}{2}(\mu\lambda|\nu\sigma)
\right].

\f]

This is the RHF analogue of the toy update \f$H = H_0 + \alpha P\f$. The
difference is that the feedback now comes from real molecular integrals rather
than an artificial scalar coupling.

## What the Code in This Repository Actually Does

For \f$\mathrm{H}_2\f$, the repository has two relevant examples:

| File | What it actually does |
|------|------------------------|
| `examples/h2-rhf.cpp` | Builds \f$\mathrm{H}_2\f$ in `cpu::Frame`, computes STO-3G integrals with Atomiverse, and solves RHF on the CPU |
| `examples/h2_rhf.cu` | Builds the same `cpu::Frame`, computes the same integrals on the CPU, runs a CPU reference solve, then runs the Driver-based GPU RHF solve |

This matters because the old wording was misleading. The CPU example is **not**
a handwritten standalone chemistry program. It already uses the Atomiverse
library:

- `cpu::scf::build_integrals_sto3g(frame)`
- `cpu::scf::rhf_solve(integrals, tol, max_iter)`

So the real comparison in this repository is:

- same geometry,
- same integral data,
- CPU solve versus GPU Driver solve.

That is a much better comparison anyway, because it isolates the backend and
execution model instead of mixing in unrelated example code.

## Step 1: Build the Molecular Geometry

Both CPU and GPU examples start from a user-facing `cpu::Frame`:

```cpp
cpu::Frame frame(2);
frame.x[0] = Real(0);     frame.y[0] = Real(0); frame.z[0] = Real(0);
frame.x[1] = Real(0.74);  frame.y[1] = Real(0); frame.z[1] = Real(0);
frame.mass[0] = Real(1.008); frame.mass[1] = Real(1.008);
frame.type_id[0] = 1; frame.type_id[1] = 1;
frame.box = Mat3::cubic(100.0);
```

Several details are worth making explicit:

- positions are given in Angstrom,
- `type_id = 1` means hydrogen, so the nuclear charge is \f$Z=1\f$,
- a very large cubic box is used only as a harmless container here,
- the user-facing data structure is still `cpu::Frame`, even for the GPU path.

That last point is central to Atomiverse: the user describes the problem once,
then the selected backend handles execution.

## Step 2: Compute Molecular Integrals on the CPU

The next line is

```cpp
auto integrals = cpu::scf::build_integrals_sto3g(frame);
```

This returns an `IntegralPack` containing:

- \f$S\f$,
- \f$T\f$,
- \f$V\f$,
- \f$H_{\mathrm{core}} = T + V\f$,
- the four-index ERI tensor,
- the nuclear repulsion energy \f$V_{\mathrm{nuc}}\f$,
- basis and electron counts.

In the codebase this API is declared in
[`include/atomiverse/cpu/scf/integrals.h`](../../include/atomiverse/cpu/scf/integrals.h).

For a beginner, the important conceptual point is:

- before you can start the RHF iteration, you need fixed molecular data derived
  from the geometry and basis set,
- that fixed data is what the integral builder provides.

For \f$\mathrm{H}_2\f$ in STO-3G, the problem is tiny:

- \f$n_{\mathrm{basis}} = 2\f$,
- \f$n_{\mathrm{electrons}} = 2\f$,
- \f$n_{\mathrm{occ}} = 1\f$.

So the matrices are small enough to print and inspect directly, which the
example does.

## Step 3: Solve RHF on the CPU

The CPU example continues with

```cpp
auto result = cpu::scf::rhf_solve(integrals, tol, 100);
```

This is the CPU reference implementation of the RHF algorithm. Its job is to:

1. compute \f$X = S^{-1/2}\f$,
2. build the Fock matrix from the current density,
3. transform to the orthogonal basis,
4. diagonalize,
5. back-transform the coefficients,
6. rebuild the density,
7. compute the energy,
8. test convergence.

That is the whole RHF loop in one library call.

This is a much more accurate description of the current code than pretending the
example manually re-implements all primitive Gaussian formulas itself.

## The RHF Iteration as a Sequence of Steps

The mathematical loop performed by `rhf_solve()` and by the GPU driver is:

\f[

F = H_{\mathrm{core}} + G(P)

\f]

\f[

F' = X^T F X

\f]

\f[

F' C' = C' \varepsilon

\f]

\f[

C = X C'

\f]

\f[

P = 2 C_{\mathrm{occ}} C_{\mathrm{occ}}^T

\f]

\f[

E_{\mathrm{elec}} =
\frac{1}{2}\mathrm{Tr}\left[P(H_{\mathrm{core}} + F)\right]

\f]

\f[

E_{\mathrm{total}} = E_{\mathrm{elec}} + V_{\mathrm{nuc}}.

\f]

Convergence is checked from the density change:

\f[

\left\|P^{(k+1)} - P^{(k)}\right\|_F < \mathrm{tol}.

\f]

This should look structurally familiar from the toy tutorial:

- build an effective matrix from the current density,
- diagonalize,
- rebuild the density,
- repeat until stable.

The difference is that every piece now has a physical interpretation.

## Step 4: Run the Same RHF Solve Through the GPU Driver

The CUDA example
[`examples/h2_rhf.cu`](../../examples/h2_rhf.cu)
uses the same `IntegralPack`, but schedules the solve through Atomiverse's SCF
driver:

```cpp
scf::SCF solver;
solver.add(scf::RHFFock{integrals});
solver.add(scf::EigenSolver{});
solver.add(scf::RHFDensity{integrals.n_occ()});
solver.add(scf::ConvergenceCheck{tol, false});
solver.on(scf::Phase::end_of_step).add(Thermo{1});
solver.build();
solver.run(frame, 100);
```

This is the cleanest way to understand the framework:

- `RHFFock` replaces the toy Hamiltonian builder,
- `EigenSolver` is reused,
- `RHFDensity` replaces the toy density builder,
- `ConvergenceCheck` is reused,
- the `scf::SCF` phase order stays the same.

So the major architectural claim is not "a completely different algorithm."
It is:

- the same SCF driver structure survives,
- only the problem-specific components become richer.

## CPU Path and GPU Path: What Is Shared, What Is Different?

Here is the practical split:

| Part of the workflow | CPU example | GPU example |
|----------------------|-------------|-------------|
| Build `cpu::Frame` | yes | yes |
| Compute `IntegralPack` on CPU | yes | yes |
| Solve RHF with `cpu::scf::rhf_solve` | yes | used as reference check |
| Solve RHF with `scf::SCF` driver | no | yes |

That means the GPU example is not "chemistry on GPU from scratch." It is:

- chemistry data prepared on CPU,
- iterative linear algebra and density updates executed through GPU components.

This is a sensible design for an early implementation because it keeps the
physics and backend concerns separated.

## How the GPU Components Map onto the RHF Mathematics

The relevant headers are:

- [`include/atomiverse/cuda/scf/rhf_fock.cuh`](../../include/atomiverse/cuda/scf/rhf_fock.cuh)
- [`include/atomiverse/cuda/scf/rhf_density.cuh`](../../include/atomiverse/cuda/scf/rhf_density.cuh)
- [`include/atomiverse/cuda/scf/scf.cuh`](../../include/atomiverse/cuda/scf/scf.cuh)

The mapping is:

| RHF step | Atomiverse phase | Component |
|----------|------------------|-----------|
| build \f$F = H_{\mathrm{core}} + G(P)\f$ and transform to \f$F'\f$ | `build_hamiltonian` | `scf::RHFFock` |
| diagonalize \f$F'\f$ | `diagonalize` | `scf::EigenSolver` |
| back-transform \f$C = XC'\f$, rebuild \f$P\f$, compute energy | `build_density` | `scf::RHFDensity` |
| test \f$\|P - P_{\mathrm{old}}\|_F\f$ | `convergence_check` | `scf::ConvergenceCheck` |

This is a stronger and clearer statement than the old "three levels" section,
because it matches the code directly.

## What `RHFFock` Really Adds Beyond the Toy Model

In the toy example, the Hamiltonian builder only had to evaluate

\f[

H = H_0 + \alpha P.

\f]

In RHF, `RHFFock` must do much more:

1. own or upload the fixed integral data,
2. build the two-electron contribution \f$G(P)\f$,
3. add it to \f$H_{\mathrm{core}}\f$,
4. transform the resulting Fock matrix into the orthogonal basis.

That is why RHF feels more complicated. The SCF loop is still the same shape,
but the "build effective matrix" step now contains real chemistry.

## What `RHFDensity` Really Adds Beyond the Toy Model

In the toy example, density rebuild was just

\f[

P = C_{\mathrm{occ}} C_{\mathrm{occ}}^T.

\f]

In RHF, `RHFDensity` must:

1. back-transform coefficients from the orthogonal basis,
2. include the closed-shell factor of \f$2\f$,
3. compute the physically meaningful RHF energy expression.

So again, the loop shape is the same, but the physical content is richer.

## What Result Should You Expect?

For \f$\mathrm{H}_2\f$ near \f$R = 0.74\f$ Angstrom, the examples are intended to
recover the familiar STO-3G RHF energy near

\f[

E_{\mathrm{total}} \approx -1.1168 \text{ Hartree}.

\f]

The exact printed value may vary slightly with precision or implementation
details, but the CPU and GPU paths should agree closely.

That agreement is the practical check that matters in this repository.

## File Guide

If you want to keep reading code after this tutorial, use this order:

1. [`examples/h2-rhf.cpp`](../../examples/h2-rhf.cpp)
2. [`include/atomiverse/cpu/scf/integrals.h`](../../include/atomiverse/cpu/scf/integrals.h)
3. [`include/atomiverse/cpu/scf/rhf_solve.h`](../../include/atomiverse/cpu/scf/rhf_solve.h)
4. [`examples/h2_rhf.cu`](../../examples/h2_rhf.cu)
5. [`include/atomiverse/cuda/scf/rhf_fock.cuh`](../../include/atomiverse/cuda/scf/rhf_fock.cuh)
6. [`include/atomiverse/cuda/scf/rhf_density.cuh`](../../include/atomiverse/cuda/scf/rhf_density.cuh)

That order mirrors how a new reader usually understands the system:

- first the user-facing example,
- then the CPU chemistry APIs,
- then the Driver-based GPU realization.

## References

1. A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*, Dover, 1996.
2. W. J. Hehre, R. F. Stewart, and J. A. Pople, *J. Chem. Phys.* **51**, 2657
   (1969).
