# Tutorials

These chapters are written as a short book. They are organized in three
arcs that share the same `cpu::Frame` data model, the same
Driver/Component vocabulary, and the same units (eV, Å, fs, amu).

If you have not built the project yet, start with
[Getting Started](../user-guide/getting-started.md). If you are
completely new to the codebase, [Tutorial 0: A Tour of the CPU
Backend](00-cpu-overview.md) is the orientation chapter — it explains
the three driver categories and the role/capability split before
*any* physics is introduced.

## Arc 1: SCF and Quantum Forces (Chapters 1–4, 6, 8)

The same self-consistent loop appears in every chapter, but the meaning
of the matrices becomes richer as you move forward: first a
stripped-down feedback problem, then Hartree-Fock for a real molecule,
then open-shell spin physics, and finally Kohn-Sham DFT on a grid.

The core reading order is
[Tutorial 1: SCF From First Principles](01-scf.md),
[Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md),
[Tutorial 3: Unrestricted Hartree-Fock (UHF)](03-uhf.md), and
[Tutorial 4: Unrestricted Kohn-Sham DFT (UKS)](04-uks.md).
Those four chapters form one continuous path from "why does SCF need
iteration?" to "how does Atomiverse evaluate a real exchange-correlation
functional?"

[Tutorial 6: Ab Initio Molecular Dynamics](06-aimd.md) and [Tutorial 8:
Geometry Minimization with L-BFGS or CG](08-geometry-minimization.md)
both apply the SCF machinery to nuclear motion. They reuse the
same `SCFEvaluator` component as a force provider — Tutorial 6 plugs it
into `cpu::md::MD`, Tutorial 8 plugs it into `cpu::minimizer::Minimize`.
Tutorial 8 is written around `examples/05_geom_opt.cpp`;
Tutorial 6 is written around `examples/06_aimd.cpp`. The
recipe in both cases is "driver + SCFEvaluator at the force phase",
which is the architecture paying off.

If you want a one-page map of the SCF arc before opening the first
chapter, read [SCF Learning Path](../user-guide/scf-roadmap.md).

## Arc 2: Periodic Systems with Plane Waves (Chapter 5)

The molecular SCF arc uses atom-centered Gaussian basis functions.
[Tutorial 5: Plane-Wave DFT for Periodic Systems](05-pw.md) switches to the
other major basis family — plane waves — which is the natural choice for
crystals, surfaces, and bulk materials. It keeps the same SCF *shape*
(build effective Hamiltonian → diagonalize → rebuild density → check) but
introduces pseudopotentials, G-vectors, k-point sampling, FFT-accelerated
Hartree/XC potentials, the Davidson diagonalizer, and Ewald ion-ion energy.
It is written around `examples/08_pw.cpp` (bulk silicon) and needs FFTW3
(`-DATV_BUILD_PW=ON`).

## Arc 3: Classical Molecular Dynamics (Chapter 7)

[Tutorial 7: Classical Molecular Dynamics on the CPU](07-classical-md.md)
walks through a complete classical-MD setup — pair potential, bonded
terms, neighbor list, integrator, thermostat, barostat, and observers
— without any quantum chemistry. It uses the same `cpu::md::MD` driver
as Tutorial 6, with classical force-field components at the force
phase instead of `SCFEvaluator`. After reading it you will have a
working mental model for both halves of the AIMD recipe.

## Reading Order

Read each chapter with the matching numbered example open beside it.
Tutorial 1 is written around `examples/01_scf.cpp`. Tutorial 2
is written around `examples/02_rhf.cpp`, with
`examples/02_rhf.cu` as the optional CUDA translation of the
same calculation. Tutorial 3 follows `examples/03_uhf.cpp`,
and Tutorial 4 follows `examples/04_uks.cpp`. Tutorial 5 follows
`examples/08_pw.cpp`. Tutorial 6 follows `examples/06_aimd.cpp`. Tutorial 7
follows `examples/07_md_lj.cpp` (CPU; the CUDA twin is `examples/07_md_lj.cu`).
Tutorial 8 follows `examples/05_geom_opt.cpp`.

(Note: example numbers do not line up with tutorial numbers — the example files
have their own order. Follow the file named in each chapter rather than assuming
the numbers match.)

If you are deciding which chapter to start with: if quantum chemistry
is what brought you here, start at Tutorial 1. If classical force
fields and ensemble simulation are what brought you here, start at
Tutorial 7.
