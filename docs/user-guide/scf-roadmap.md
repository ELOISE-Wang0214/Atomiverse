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
matrix means physically.

## Recommended Reading Order

Read these chapters in order:

1. [Tutorial 1: SCF From First Principles](../tutorial/01-scf.md)
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

| Tutorial | Primary example | Why it matters |
|---|---|---|
| SCF introduction | `examples/tutorial/01_scf.cpp` | Learn the loop before learning chemistry |
| RHF | `examples/tutorial/02_rhf.cpp` | First physically meaningful calculation |
| UHF | `examples/tutorial/03_uhf.cpp` | Open-shell systems and spin contamination |
| UKS | `examples/tutorial/04_uks.cpp` | DFT, grids, and exchange-correlation potentials |

The examples are intentionally numbered to match the chapters. Read the theory
and the matching code side by side.

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

Keep the matching tutorial example open while you read each chapter. If you
want to understand the runtime system after that, then read the
[Developer Manual](../developer/index.md). The tutorials themselves should be
enough to learn the chemistry without dropping straight into architecture.
