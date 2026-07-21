# Tutorials

Each chapter teaches one scientific method. Install and build first if needed:
[User Guide](../user-guide/).

All examples share `cpu::Frame`, the Driver/Component pattern, and units
(eV, Å, fs, amu).

## How every chapter is organized

Every chapter (and every major section inside multi-topic chapters) follows
the same three blocks:

| Block | What it is |
|-------|------------|
| **1. Theory** | Why the method exists, equations, assumptions, limitations |
| **2. Practice** | Minimal C++ that runs in Atomiverse (members, `run`, outputs) |
| **3. Results** | Numbers from the engine **vs** literature / PySCF / other reference |

How to read:

1. Skim **Theory** with a pen, not a compiler.
2. Build and run the **Practice** snippet (or the named example binary).
3. Check **Results** against the table; do not invent numbers.

Interactive figures (MolPlot curves, MolVis **density/orbital cubes**) appear
only when they show a computed field or a quantitative trend—not bare toy
geometries.

Suggested first paths:

- Electronic structure: [Ch. 2](ch02-rhf.md) → [3](ch03-open-shell.md) → [4](ch04-basis-sets.md) → [18](ch18-scf-techniques.md)–[20](ch20-hf-properties-environment.md)
- Classical MD: [Ch. 10](ch10-classical-md.md)
- Architecture: [Ch. 0](ch00-architecture.md)

## I. Foundations

| Ch | Title |
|----|-------|
| [0](ch00-architecture.md) | Tour of the CPU backend |
| [1](ch01-foundations.md) | Systems, units, precision, validation |

## II. Hartree–Fock core

| Ch | Title | Example |
|----|-------|---------|
| [2](ch02-rhf.md) | Restricted HF single point | `02_rhf`, `16_density_cube` |
| [3](ch03-open-shell.md) | Open shells (UHF, ROHF, GHF) | `03_uhf` |
| [4](ch04-basis-sets.md) | Gaussian basis sets | basis regressions |

## III. DFT

| Ch | Title | Example |
|----|-------|---------|
| [5](ch05-kohn-sham-dft.md) | Unrestricted Kohn–Sham | `04_uks` |
| [6](ch06-hybrid-dft.md) | Hybrids and forces | `05_geom_opt` |

## IV. Correlation

| Ch | Title | Example |
|----|-------|---------|
| [7](ch07-correlation.md) | MP2 and full CI | `13_mp2`, `14_fci` |

## V. Periodic electronic structure

| Ch | Title | Example |
|----|-------|---------|
| [8](ch08-plane-wave.md) | Plane-wave DFT | `08_pw` |
| [9](ch09-gpw.md) | Gaussian and plane waves | `09_gpw` |

## VI. Classical simulation

| Ch | Title | Example |
|----|-------|---------|
| [10](ch10-classical-md.md) | Classical MD | `07_md_lj` |
| [11](ch11-ensembles.md) | Thermostats and barostats | — |
| [12](ch12-geometry-optimization.md) | Geometry optimization | `05_geom_opt` |

## VII. Analysis and response

| Ch | Title |
|----|-------|
| [13](ch13-observables.md) | Observables and trajectories |
| [14](ch14-vibrations.md) | Vibrations and thermochemistry |
| [15](ch15-cross-method.md) | Swapping evaluators |

## VIII–IX. Beyond pairs and semiempirical

| Ch | Title | Example |
|----|-------|---------|
| [16](ch16-many-body-potentials.md) | EAM and Stillinger–Weber | — |
| [17](ch17-semiempirical-am1.md) | AM1 | `15_am1` |

## X. Hartree–Fock advanced

Continuation of §II: convergence tools, large-system Coulomb, properties and
environment (Jensen-style arc).

| Ch | Title |
|----|-------|
| [18](ch18-scf-techniques.md) | Convergence of the SCF |
| [19](ch19-hf-large-systems.md) | Large molecules and the Coulomb bottleneck |
| [20](ch20-hf-properties-environment.md) | Molecular properties and environment |
| [Appendix](appendix-pyscf-oracle.md) | External reference numbers |

---

External-engine parity policy:
[`benchmarks/README.md`](https://github.com/MolCrafts/Atomiverse/blob/dev/benchmarks/README.md).
One-page SCF map: [User Guide → SCF Learning Path](../user-guide/scf-roadmap.md).
