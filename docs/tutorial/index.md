# SCF Tutorials

These chapters are meant to be read as a short book. The same self-consistent
loop appears in every chapter, but the meaning of the matrices becomes richer
as you move forward: first a stripped-down feedback problem, then Hartree-Fock for a real
molecule, then open-shell spin physics, and finally Kohn-Sham DFT on a grid.

If you have not built the project yet, start with
[Getting Started](../user-guide/getting-started.md). If you want a one-page map
of the SCF material before opening the first chapter, read
[SCF Learning Path](../user-guide/scf-roadmap.md).

The core reading order is
[Tutorial 1: SCF From First Principles](01-scf.md),
[Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md),
[Tutorial 3: Unrestricted Hartree-Fock (UHF)](03-uhf.md), and
[Tutorial 4: Unrestricted Kohn-Sham DFT (UKS)](04-uks.md).
Those four chapters form one continuous path from "why does SCF need
iteration?" to "how does Atomiverse evaluate a real exchange-correlation
functional?"

A sixth chapter applies the SCF machinery to molecular dynamics:
[Tutorial 6: Ab Initio Molecular Dynamics](06-aimd.md). It reuses the
same `SCFEvaluator` component as a force provider, this time plugged
into the `cpu::md::MD` driver. The pattern is the same as the geometry
optimization shown in `examples/tutorial/05_geom_opt.cpp` (driver +
SCFEvaluator at the force phase), only with a different driver and a
time-stepping outer loop.

Read each chapter with the matching numbered example open beside it. Tutorial 1
is written around `examples/tutorial/01_scf.cpp`. Tutorial 2 is written around
`examples/tutorial/02_rhf.cpp`, with `examples/tutorial/02_rhf.cu` as the
optional CUDA translation of the same calculation. Tutorial 3 follows
`examples/tutorial/03_uhf.cpp`, and Tutorial 4 follows
`examples/tutorial/04_uks.cpp`. Tutorial 6 follows
`examples/tutorial/06_aimd.cpp`.
