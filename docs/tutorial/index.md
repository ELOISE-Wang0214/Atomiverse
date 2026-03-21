# SCF Deep-Dive Tutorials

This sequence is now the deep-dive tutorial track inside the
[User Guide](../user-guide/index.md). It is written for readers who are new to
computational chemistry and new to Atomiverse.

Read [Getting Started](../user-guide/getting-started.md) first if you have not
built the project or run any examples yet.

If you want the high-level map first, read
[SCF Learning Path](../user-guide/scf-roadmap.md) before starting Tutorial 1.

## Reading Order

1. [Tutorial 1: Toy SCF from First Principles](01-scf.md)
2. [Tutorial 2: Restricted Hartree-Fock (RHF)](02-rhf.md)
3. [Tutorial 3: Unrestricted Hartree-Fock (UHF)](03-uhf.md)
4. [Tutorial 4: Unrestricted Kohn-Sham DFT (UKS)](04-uks.md)

## Why This Order?

Read Tutorial 1 first to understand the SCF feedback loop without chemistry
details. Then read Tutorial 2 to see how the same iteration becomes a real RHF
calculation for \f$\mathrm{H}_2\f$. Tutorial 3 extends RHF to open-shell systems
using unrestricted orbitals. Tutorial 4 replaces Hartree-Fock exchange with an
exchange-correlation functional and introduces the grid-integration ideas behind
Kohn-Sham DFT.

## Code Map

Use this map if you want to switch between the prose and the code while reading.

| Topic | Primary code target |
|---|---|
| Toy SCF loop | `examples/toy-scf.cpp`, `examples/scf.cu`, `include/atomiverse/cpu/scf/solve.h` |
| RHF | `examples/h2-rhf.cpp`, `examples/h2_rhf.cu`, `include/atomiverse/cpu/scf/rhf_solve.h` |
| UHF | `include/atomiverse/cpu/scf/uhf_solve.h` |
| UKS | `include/atomiverse/cpu/scf/uks_solve.h`, `include/atomiverse/cpu/scf/xc_functional.h` |
