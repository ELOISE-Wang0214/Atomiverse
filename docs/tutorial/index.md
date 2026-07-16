# Tutorials

These chapters are written as a single book — one continuous progression from
the engine's architecture, through localized-basis and plane-wave electronic
structure, to classical dynamics, analysis, and response properties. Every
chapter shares the same `cpu::Frame` data model, the same Driver/Component
vocabulary, and the same units (eV, Å, fs, amu), and every chapter ends by
pointing at the one architectural pattern the whole book is built on: **one loop,
many physics** — a `Driver` owns the loop, a `Component` supplies the physics,
and you change method by swapping the evaluator.

If you have not built the project yet, start with
[Getting Started](../user-guide/getting-started.md). If quantum chemistry is what
brought you here, the natural entry point is **Chapter 2**; if classical force
fields and ensemble simulation are what you want, start at **Chapter 10**. Read
each chapter with the example named in its text open beside it — the example
files under [`examples/`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/) have their own numbering (follow the
filename each chapter names, not the number).

Completely new to the codebase? [Chapter 0](ch00-architecture.md) is the
orientation chapter — it explains the drivers, their termination handlers, and
the role/capability split before any physics is introduced.

## I. Foundations

- [Chapter 0: A Tour of the CPU Backend](ch00-architecture.md) — the driver
  presets and their loops, the Driver/Component/Frame pillars, and the
  role/capability split. The map for everything that follows.
- [Chapter 1: Systems, Units, Precision, and Validation](ch01-foundations.md) —
  the `cpu::Box`/`Atom`/units model, LAMMPS-style tiered precision
  (`StorageReal`/`WorkReal`/`AccReal`, `LocalInt`/`TagInt`/`ImageInt`/`BigInt`),
  and finite-difference gradient validation. Read this before the physics arcs.

## II. Localized-Basis Electronic Structure

- [Chapter 2: Atomic Orbitals to Single-Point RHF](ch02-rhf.md) — the SCF
  fixed-point idea taken from a bare feedback loop (`examples/01_scf.cpp`) to a
  real closed-shell Hartree-Fock energy for H₂ (`examples/02_rhf.cpp`).
- [Chapter 3: Unrestricted Hartree-Fock](ch03-open-shell.md) — dropping the
  closed-shell restriction: open-shell spin, spin contamination, Mulliken
  populations (`examples/03_uhf.cpp`).
- [Chapter 4: Gaussian Basis Sets and Convergence](ch04-basis-sets.md) — STO-3G
  → def2-SVP and how the energy converges with the basis.

## III. Density-Functional Theory

- [Chapter 5: Unrestricted Kohn-Sham DFT (UKS)](ch05-kohn-sham-dft.md) — the
  density representation, the KS flow, exchange-correlation functionals, and
  integration grids (`examples/04_uks.cpp`).
- [Chapter 6: Hybrid Functionals and Analytic Gradients](ch06-hybrid-dft.md) —
  B3LYP/PBE0 exact-exchange admixture and `scf_gradient`.

## IV. Correlation Methods

- [Chapter 7: Beyond Mean Field](ch07-correlation.md) — a **status chapter**:
  the correlation energy, the MP2/CI/coupled-cluster hierarchy, why DFT is
  Atomiverse's practical answer, and an honest statement that wavefunction
  correlation methods are planned, not implemented.

## V. Periodic / Plane-Wave Electronic Structure

- [Chapter 8: Plane-Wave DFT for Periodic Systems](ch08-plane-wave.md) — the
  other major basis family: pseudopotentials, G-vectors, k-points, FFT
  Hartree/XC, the Davidson diagonalizer, Ewald ion-ion energy. Bulk silicon
  (`examples/08_pw.cpp`; needs `-DATV_BUILD_PW=ON`).
- [Chapter 9: Gaussian and Plane Waves (GPW)](ch09-gpw.md) — CP2K Quickstep's
  method: a Gaussian wavefunction with FFT-grid Hartree/XC and GTH
  pseudopotentials, validated on H₂O and a small periodic wurtzite-ZnO cell
  (`examples/09_gpw.cpp`; needs `-DATV_BUILD_GPW=ON`).
  Experimental (closed-shell, Γ-point, PADE/PBE).

## VI. Classical Simulation

- [Chapter 10: Classical Molecular Dynamics on the CPU](ch10-classical-md.md) —
  a complete classical-MD setup: pair + bonded potentials, neighbor list, PBC,
  integrator, and NVE energy conservation (`examples/07_md_lj.cpp`; CUDA twin
  `examples/07_md_lj.cu`).
- [Chapter 11: Ensembles — Thermostats and Barostats](ch11-ensembles.md) —
  NVT/NPT with Nosé-Hoover, Langevin, CSVR, and Berendsen.
- [Chapter 12: Geometry Minimization with L-BFGS or CG](ch12-geometry-optimization.md)
  — optimizing on a potential-energy surface, classical or *ab initio*
  (`examples/05_geom_opt.cpp`).

## VII. Analysis & Observables

- [Chapter 13: Observables and Trajectory Analysis](ch13-observables.md) — energy
  conservation, the radial distribution function (`ATV::RDF`), and transport via
  `ATV::MSD` (Einstein self-diffusion) and `ATV::VACF` (Green-Kubo / VDOS). The
  engine measures through an async snapshot ring; molrs owns the analysis math.

## VIII. Response & Cross-Method Workflows

- [Chapter 14: Vibrational Analysis and Thermochemistry](ch14-vibrations.md) —
  the finite-difference Hessian (`cpu::vib::Hessian`), harmonic frequencies, and
  thermochemistry.
- [Chapter 15: Cross-Method Workflows — Swap the Evaluator](ch15-cross-method.md)
  — the payoff chapter: one system driven through single-point, optimization, and
  MD with classical then quantum forces, one line different each time. Includes
  the *ab initio* MD (AIMD) deep-dive — ASPC, ensembles, and pitfalls
  (`examples/06_aimd.cpp`).

## IX. Beyond Pair Potentials

- [Chapter 16: Many-Body Potentials — Metals and Covalent Solids](ch16-many-body-potentials.md)
  — why a pair potential is *qualitatively* wrong for a metal or a covalent solid
  (close packing, the Cauchy relation, the vacancy paradox), and the two potentials
  that fix it: the embedded-atom method (`cpu::EAM`, Sutton-Chen copper) and
  Stillinger-Weber silicon (`cpu::StillingerWeber`). Covers the full neighbor list,
  the two-pass force, and why every force needs a finite-difference oracle.

## X. Semiempirical Quantum Chemistry

- [Chapter 17: Semiempirical AM1 — Cheap Quantum Forces](ch17-semiempirical-am1.md)
  — NDDO/AM1 as `Semiempirical` kind + `NDDOSolver` (not a Driver): single-point,
  ΔHf, analytic forces, Minimize, and short BOMD (`examples/15_am1.cpp`).
  DFTB/xTB stay out of family.

---

The validation policy — what each acceptance check proves, and against which reference engine — lives in
[`benchmarks/README.md`](https://github.com/MolCrafts/Atomiverse/blob/dev/benchmarks/README.md). For a
one-page map of the SCF arc, see the
[SCF Learning Path](../user-guide/scf-roadmap.md).
