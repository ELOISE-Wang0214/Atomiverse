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
5. [Tutorial 5: Plane-Wave DFT for Periodic Systems](../tutorial/05-pw.md)

Each chapter adds exactly one new idea.

Tutorial 1 teaches feedback and convergence with a matrix model that has almost
no chemistry in it.

Tutorial 2 replaces the toy matrix with real molecular integrals and introduces
the RHF closed-shell approximation.

Tutorial 3 removes the "same spatial orbital for both spins" restriction and
shows how open-shell systems require separate alpha and beta channels.

Tutorial 4 keeps the spin-resolved structure of UHF but replaces explicit
Hartree-Fock exchange with exchange-correlation functionals evaluated on a grid.

Tutorial 5 moves from molecules to periodic solids, replacing atom-centered
Gaussian basis sets with plane waves and introducing pseudopotentials, FFT
acceleration, k-point sampling, and Ewald summation.

## Which Example Matches Which Tutorial?

| Tutorial | Primary example | Why it matters |
|---|---|---|
| SCF introduction | `examples/tutorial/01_scf.cpp` | Learn the loop before learning chemistry |
| RHF | `examples/tutorial/02_rhf.cpp` | First physically meaningful calculation |
| UHF | `examples/tutorial/03_uhf.cpp` | Open-shell systems and spin contamination |
| UKS | `examples/tutorial/04_uks.cpp` | DFT, grids, and exchange-correlation potentials |
| Plane-Wave DFT | `cpu::pw::solve_si` | Periodic systems, pseudopotentials, FFT acceleration |

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

Pair each tutorial with one code-reading target:

- after Tutorial 1, read [`cpu::scf::solve`](../../include/atomiverse/cpu/scf/solve.h)
- after Tutorial 2, read [`IntegralPack`](../../include/atomiverse/cpu/scf/integrals.h) and [`rhf_solve`](../../include/atomiverse/cpu/scf/rhf_solve.h)
- after Tutorial 3, read [`uhf_solve`](../../include/atomiverse/cpu/scf/uhf_solve.h)
- after Tutorial 4, read [`uks_solve`](../../include/atomiverse/cpu/scf/uks_solve.h) and [`xc_functional.h`](../../include/atomiverse/cpu/scf/xc_functional.h)
- after Tutorial 5, read [`pw_solver.h`](../../include/atomiverse/cpu/pw/pw_solver.h) and [`internal/pw_workspace.h`](../../include/atomiverse/cpu/pw/internal/pw_workspace.h)

If you want to understand how the same ideas map onto the runtime system, then
read the [Developer Manual](../developer/index.md) after finishing Tutorial 2.
That is usually the point where `Driver`, `Plan`, and `Component` start to make
sense instead of feeling abstract.
