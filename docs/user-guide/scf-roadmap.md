# SCF Learning Path

This page explains how to study the SCF part of Atomiverse without getting lost
in either chemistry jargon or internal implementation details.

The key idea is that all four SCF tutorials describe the same computational
pattern:

1. start from a density guess
2. build an effective matrix from that guess
3. diagonalize the matrix
4. rebuild the density
5. stop when the density stops changing

What changes from chapter to chapter is not the loop itself, but what the
matrix means physically and which Atomiverse APIs are responsible for building
it.

## Recommended Reading Order

Read these chapters in order:

1. [Tutorial 1: Toy SCF from First Principles](../tutorial/01-scf.md)
2. [Tutorial 2: Restricted Hartree-Fock (RHF)](../tutorial/02-rhf.md)
3. [Tutorial 3: Unrestricted Hartree-Fock (UHF)](../tutorial/03-uhf.md)
4. [Tutorial 4: Unrestricted Kohn-Sham DFT (UKS)](../tutorial/04-uks.md)

Each chapter adds exactly one new idea.

Tutorial 1 teaches feedback and convergence with a matrix model that has almost
no chemistry in it.

Tutorial 2 replaces the toy matrix with real molecular integrals and introduces
the RHF closed-shell approximation.

Tutorial 3 removes the “same spatial orbital for both spins” restriction and
shows how open-shell systems require separate alpha and beta channels.

Tutorial 4 keeps the spin-resolved structure of UHF but replaces explicit
Hartree-Fock exchange with exchange-correlation functionals evaluated on a grid.

## Which Example Matches Which Tutorial?

| Tutorial | Primary example or API | Why it matters |
|---|---|---|
| Toy SCF | `examples/toy-scf.cpp`, `examples/scf.cu`, `cpu::scf::solve` | Learn the loop before learning chemistry |
| RHF | `examples/h2-rhf.cpp`, `examples/h2_rhf.cu`, `cpu::scf::rhf_solve` | First physically meaningful calculation |
| UHF | `cpu::scf::uhf_solve` | Open-shell systems and spin contamination |
| UKS | `cpu::scf::uks_solve` | DFT, grids, and exchange-correlation potentials |

The examples and solver APIs are intentionally paired. The example shows how to
assemble a complete calculation. The solver header shows the exact input and
output structures you need to understand.

## When To Use CPU Solvers Versus CUDA SCF Components

If you are learning the method or validating a result, start with the CPU
solvers. The CPU code is easier to read because the linear algebra and the SCF
logic sit close together.

If you are learning how Atomiverse expresses workflows in its staged execution
model, read the CUDA SCF component path next. That path is more about software
architecture: the SCF loop is split into `build_hamiltonian`, `diagonalize`,
`build_density`, and `convergence_check` phases, and each phase is implemented
by one or more components.

## What To Read Alongside The Tutorials

Pair each tutorial with one code-reading target:

- after Tutorial 1, read [`cpu::scf::solve`](../../include/atomiverse/cpu/scf/solve.h)
- after Tutorial 2, read [`IntegralPack`](../../include/atomiverse/cpu/scf/integrals.h) and [`rhf_solve`](../../include/atomiverse/cpu/scf/rhf_solve.h)
- after Tutorial 3, read [`uhf_solve`](../../include/atomiverse/cpu/scf/uhf_solve.h)
- after Tutorial 4, read [`uks_solve`](../../include/atomiverse/cpu/scf/uks_solve.h) and [`xc_functional.h`](../../include/atomiverse/cpu/scf/xc_functional.h)

If you want to understand how the same ideas map onto the runtime system, then
read the [Developer Manual](../developer/index.md) after finishing Tutorial 2.
That is usually the point where `Driver`, `Plan`, and `Component` start to make
sense instead of feeling abstract.
