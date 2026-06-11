# Tutorial 1: SCF From First Principles

This chapter is about one idea only: sometimes the matrix you want to
diagonalize depends on the occupied eigenvectors that come out of that same
diagonalization. When that happens, one solve is never enough. You must guess a
solution, build a matrix from that guess, solve again, and repeat until the
input and output agree. That feedback loop is what "self-consistent field" or
SCF means.

Open [`examples/01_scf.cpp`](../../examples/01_scf.cpp) while
you read. The calculation in that file is deliberately stripped down, not
because the code is the main point, but because a realistic molecular example
would hide the SCF pattern under too many new ideas at once. Later chapters
will add chemistry. This first chapter keeps only the feedback loop itself.

## Why SCF Has To Iterate

Suppose there is a fixed symmetric matrix $H_0$, a current density-like
matrix $P$, and a scalar coupling strength $\alpha$. The stripped-down model says
that the matrix we solve at iteration $k$ is

$$
H^{(k)} = H_0 + \alpha P^{(k)}.
$$

Once $H^{(k)}$ is available, we diagonalize it,

$$
H^{(k)} C^{(k)} = C^{(k)} \varepsilon^{(k)},
$$

take the occupied columns of $C^{(k)}$, and rebuild the next density as

$$
P^{(k+1)} = C_{\mathrm{occ}}^{(k)} \left(C_{\mathrm{occ}}^{(k)}\right)^T.
$$

Nothing in this loop is chemically realistic. That is intentional. The point is
to make the fixed-point structure obvious. The current density changes the next
matrix, the next matrix changes the occupied subspace, and the occupied
subspace changes the next density. The calculation stops only when the density
hardly moves anymore:

$$
\left\|P^{(k+1)} - P^{(k)}\right\|_F < \mathrm{tol}.
$$

If you keep only one picture from this chapter, keep this one: SCF is an
attempt to make the matrix used for the solve consistent with the orbitals that
the solve itself produces.

## Reading The Example And The Equations Together

The example begins by building a fixed matrix $H_0$ and setting the first
density guess to zero. That corresponds exactly to the mathematics above: at the
start of the calculation, the only thing you know is the fixed part of the
problem and a guess for the occupied subspace.

Once the loop starts, the file follows the same story as the equations. It
forms the current matrix $H^{(k)}$, diagonalizes it, keeps the occupied
columns, and rebuilds the next density from them. That is why the most
important object in this chapter is not the matrix $H_0$ by itself, but the
projector
$P = C_{\mathrm{occ}} C_{\mathrm{occ}}^T$. It is the mathematical record of
which states are occupied, and it is the object that gets fed back into the
next iteration.

The example also prints an energy-like scalar
$E = \mathrm{Tr}(H P)$ and the density residual at every step. In this first
chapter those numbers are only diagnostics. Their job is to let you watch the
iteration settle. In later chapters, almost the same pattern will turn into a
real molecular total energy.

If you build and run the file with

```bash
cmake -B build -DATV_BUILD_EXAMPLES=ON
cmake --build build --target 01_scf
./build/examples/01_scf
```

you should see the residual shrink until the loop stops changing. That is the
whole lesson of this chapter: SCF is not a single diagonalization. It is a
search for a density that reproduces itself.

## What To Take Into Tutorial 2

By the end of this introductory model you should be comfortable with two facts. First,
SCF is a fixed-point problem rather than a one-shot eigensolve. Second, the
problem is easiest to learn when the example stays as close as possible to the
matrix operations themselves. In
[Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md) the loop will stop being
abstract. The model matrix $H_0 + \alpha P$ will be replaced by a Fock
matrix built from real molecular integrals, and the trace-like energy will turn
into the first physically meaningful number in the sequence.
