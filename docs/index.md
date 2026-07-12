---
title: Atomiverse
hide: [navigation, toc]
hero:
  kicker: Manual
  title: Atomiverse
  description: >-
    One Driver/Component API across classical MD, Monte Carlo, geometry
    optimization, vibrational analysis, and electronic structure — on CPU and
    CUDA. Change the physics by swapping a component, not the program.
  actions:
    - { label: Get started, href: user-guide/getting-started/, style: primary }
    - { label: Tutorial, href: tutorial/ }
  install:
    label: Build
    methods:
      - { label: cmake, command: "cmake -B build -DATV_BUILD_RHF=ON && cmake --build build -j" }
      - { label: python, command: "pip install ./python" }
      - { label: source, command: "git clone https://github.com/MolCrafts/Atomiverse" }
---

## The whole engine in one idea

A **Driver** owns the loop. A **Component** supplies one step of physics. A
**Frame** holds the state. Swap the evaluator and the same loop becomes a
different calculation.

```cpp
cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{0.5});          // fs

sim.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});            // classical MD
// or:
sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{cfg}); // ab-initio MD

sim.build();
sim.run(frame, 1000);
```

## What Atomiverse computes

| Method | What it gives you |
| ------ | ----------------- |
| [Classical MD](tutorial/ch10-classical-md.md) | LJ and bonded styles, NVE/NVT/NPT, neighbor lists, PME |
| [Optimization](tutorial/ch12-geometry-optimization.md) | L-BFGS and conjugate gradient over any force provider |
| [Monte Carlo](user-guide/methods.md#monte-carlo) | NVT Metropolis sampling on the same force pipeline |
| [Gaussian SCF](tutorial/ch02-rhf.md) | RHF, UHF, and Kohn–Sham DFT for molecules — with gradients |
| [Plane-wave DFT](tutorial/ch08-plane-wave.md) | Kohn–Sham LDA for periodic crystals, k-points, Davidson |
| [GPW DFT](tutorial/ch09-gpw.md) | Gaussian and plane waves: GTH pseudopotentials, PBE, D3 |
| [Vibrations](tutorial/ch14-vibrations.md) | Finite-difference Hessian and thermochemistry |
| [Semiempirics](user-guide/methods.md#semiempirical) | NDDO/AM1, plus AM1-BCC charge assignment |
| [Observables](tutorial/ch13-observables.md) | RDF, virial stress, Mulliken charges, densities |

Which methods are implemented, on which backend, and how far — that is the
[Scientific Methods Roadmap](roadmap/index.md).

## Find the right page

1. **[Getting started](user-guide/getting-started.md)** — from a fresh checkout
   to a converged calculation.
2. **[Methods](user-guide/methods.md)** — every method the engine ships, and
   which driver runs it.
3. **[Tour of the backend](tutorial/ch00-architecture.md)** — driver categories,
   roles, capabilities, phases.
4. **[The tutorial book](tutorial/index.md)** — sixteen chapters, from atomic
   orbitals to ab-initio MD.
5. **[Architecture](developer/architecture.md)** — how the facade drives two
   backends from one API.

## Three ways to read this

**[User Guide](user-guide/index.md)** — build it, run it, understand the output.
Task-oriented, no internals.

**[Tutorial](tutorial/index.md)** — one continuous book. Each chapter pairs with
a program under `examples/`.

**[Developer](developer/index.md)** — Driver/Component/Frame internals, the MD
step model, and how to extend them.
