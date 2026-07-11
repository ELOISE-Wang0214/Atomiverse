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

<h1 class="molcrafts-sr-only">Atomiverse</h1>

<div class="molcrafts-manual-home" markdown>

<section class="molcrafts-manual-section molcrafts-manual-section--compact" markdown>
  <div class="molcrafts-manual-section__header" markdown>
    <span class="molcrafts-manual-eyebrow">One loop, many physics</span>
    ## The whole engine in one idea

    A **Driver** owns the loop. A **Component** supplies one step of physics. A
    **Frame** holds the state. Swap the evaluator and the same loop becomes a
    different calculation.
  </div>

```cpp
cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{0.5});          // fs

sim.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});            // classical MD
// or:
sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{cfg}); // ab-initio MD

sim.build();
sim.run(frame, 1000);
```
</section>

<section class="molcrafts-manual-section molcrafts-manual-section--stack" markdown>
  <div class="molcrafts-manual-section__header" markdown>
    <span class="molcrafts-manual-eyebrow">Capabilities</span>
    ## What Atomiverse computes
  </div>
  <div class="molcrafts-manual-grid molcrafts-manual-grid--cols-3">
    <a href="tutorial/ch10-classical-md/"><strong>Classical MD</strong><em>LJ + bonded styles, NVE/NVT/NPT, neighbor lists, PME.</em></a>
    <a href="tutorial/ch12-geometry-optimization/"><strong>Optimization</strong><em>L-BFGS and conjugate gradient over any force provider.</em></a>
    <a href="user-guide/methods/#monte-carlo"><strong>Monte Carlo</strong><em>NVT Metropolis sampling on the same force pipeline.</em></a>
    <a href="tutorial/ch02-rhf/"><strong>Gaussian SCF</strong><em>RHF, UHF, and Kohn-Sham DFT for molecules — with gradients.</em></a>
    <a href="tutorial/ch08-plane-wave/"><strong>Plane-wave DFT</strong><em>Kohn-Sham LDA for periodic crystals, k-points, Davidson.</em></a>
    <a href="tutorial/ch09-gpw/"><strong>GPW DFT</strong><em>Gaussian-and-plane-waves: GTH pseudopotentials, PBE, D3.</em></a>
    <a href="tutorial/ch14-vibrations/"><strong>Vibrations</strong><em>Finite-difference Hessian and thermochemistry.</em></a>
    <a href="user-guide/methods/#semiempirical"><strong>Semiempirics</strong><em>NDDO/AM1, plus AM1-BCC charge assignment.</em></a>
    <a href="tutorial/ch13-observables/"><strong>Observables</strong><em>RDF, virial stress, Mulliken charges, densities.</em></a>
  </div>
</section>

<section class="molcrafts-manual-section" markdown>
  <div class="molcrafts-manual-section__header" markdown>
    <span class="molcrafts-manual-eyebrow">Start here</span>
    ## Find the right page
  </div>
  <div class="molcrafts-manual-index">
    <a href="user-guide/getting-started/"><span>01</span><strong>Getting started</strong><em>From a fresh checkout to a converged calculation.</em></a>
    <a href="user-guide/methods/"><span>02</span><strong>Methods</strong><em>Every method the engine ships, and which driver runs it.</em></a>
    <a href="tutorial/ch00-architecture/"><span>03</span><strong>Tour of the backend</strong><em>Driver categories, roles, capabilities, phases.</em></a>
    <a href="tutorial/"><span>04</span><strong>The tutorial book</strong><em>Sixteen chapters, from atomic orbitals to ab-initio MD.</em></a>
    <a href="developer/architecture/"><span>05</span><strong>Architecture</strong><em>How the facade drives two backends from one API.</em></a>
  </div>
</section>

<section class="molcrafts-manual-section molcrafts-manual-section--stack" markdown>
  <div class="molcrafts-manual-section__header" markdown>
    <span class="molcrafts-manual-eyebrow">Manuals</span>
    ## Three ways to read this
  </div>
  <div class="molcrafts-doc-map">
    <section markdown>
### [User Guide](user-guide/index.md)
Build it, run it, understand the output. Task-oriented, no internals.
    </section>
    <section markdown>
### [Tutorial](tutorial/index.md)
One continuous book. Each chapter pairs with a program under `examples/`.
    </section>
    <section markdown>
### [Developer](developer/index.md)
Driver/Component/Frame internals, the MD step model, and how to extend them.
    </section>
  </div>
</section>

</div>
