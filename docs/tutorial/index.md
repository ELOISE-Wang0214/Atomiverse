# Tutorials

Each chapter teaches one scientific method. Install and project layout live in
the [User Guide](../user-guide/). Every chapter uses `cpu::Frame`, the
Driver/Component pattern, and units (eV, Å, fs, amu).

## How every chapter is organized

| Block | What it is |
|-------|------------|
| **Theory** | Why the method exists, equations, assumptions, limitations |
| **Practice** | Minimal C++ (members, `run`, outputs) |
| **Results** | Engine numbers vs literature / PySCF / other reference |

Read Theory with a pen; type or adapt the Practice snippet; check Results
against the table. Figures appear only for computed fields or quantitative
trends—not bare toy geometries.

Suggested paths:

- Electronic structure: [Ch. 2](ch02-rhf.md) → [3](ch03-open-shell.md) →
  [4](ch04-basis-sets.md) → [18](ch18-scf-techniques.md)–[20](ch20-hf-properties-environment.md)
- Classical MD: [Ch. 10](ch10-classical-md.md)
- Architecture: [Ch. 0](ch00-architecture.md)

## I. Foundations

| Ch | Title |
|----|-------|
| [0](ch00-architecture.md) | CPU Backend |
| [1](ch01-foundations.md) | Systems and Units |

## II. Hartree–Fock core

| Ch | Title |
|----|-------|
| [2](ch02-rhf.md) | Restricted Hartree–Fock |
| [3](ch03-open-shell.md) | Unrestricted Hartree–Fock |
| [4](ch04-basis-sets.md) | Basis Sets |

## III. DFT

| Ch | Title |
|----|-------|
| [5](ch05-kohn-sham-dft.md) | Kohn–Sham DFT |
| [6](ch06-hybrid-dft.md) | Hybrid Functionals |

## IV. Correlation

| Ch | Title |
|----|-------|
| [7](ch07-correlation.md) | Correlation Methods |

## V. Periodic electronic structure

| Ch | Title |
|----|-------|
| [8](ch08-plane-wave.md) | Plane-Wave DFT |
| [9](ch09-gpw.md) | Gaussian and Plane Waves |

## VI. Classical simulation

| Ch | Title |
|----|-------|
| [10](ch10-classical-md.md) | Classical Molecular Dynamics |
| [11](ch11-ensembles.md) | Ensembles |
| [12](ch12-geometry-optimization.md) | Geometry Optimization |

## VII. Analysis and response

| Ch | Title |
|----|-------|
| [13](ch13-observables.md) | Observables |
| [14](ch14-vibrations.md) | Vibrations |
| [15](ch15-cross-method.md) | Cross-Method Workflows |

## VIII–IX. Beyond pairs and semiempirical

| Ch | Title |
|----|-------|
| [16](ch16-many-body-potentials.md) | Many-Body Potentials |
| [17](ch17-semiempirical-am1.md) | Semiempirical AM1 |

## X. Hartree–Fock advanced

| Ch | Title |
|----|-------|
| [18](ch18-scf-techniques.md) | SCF Convergence |
| [19](ch19-hf-large-systems.md) | Large Molecules |
| [20](ch20-hf-properties-environment.md) | Properties and Environment |
| [Appendix](appendix-pyscf-oracle.md) | External References |

---

SCF reading map: [User Guide → SCF Learning Path](../user-guide/scf-roadmap.md).
