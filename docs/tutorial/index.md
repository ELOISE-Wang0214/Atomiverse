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

## Method tutorials

Campaign handbook pages (technical-report style, one per regression slug).
Full index and template: [methods/README](methods/README.md).

| Method | Page | Regression |
|--------|------|------------|
| UMP2 | [ump2-o2-sto3g](methods/ump2-o2-sto3g.md) | `regressions/scf/ump2-o2-sto3g` |
| Hybrid TDDFT-TDA (B3LYP) | [h2-b3lyp-hybrid-tda](methods/h2-b3lyp-hybrid-tda.md) | `…/h2-b3lyp-hybrid-tda` |
| Frozen-core RMP2 | [frozen-core-mp2-h2o](methods/frozen-core-mp2-h2o.md) | `…/frozen-core-mp2-h2o` |
| CIS(D) | [cisd-h2co-sto3g](methods/cisd-h2co-sto3g.md) | `…/cisd-h2co-sto3g` |
| BLYP GGA | [blyp-water-sto3g](methods/blyp-water-sto3g.md) | `…/blyp-water-sto3g` |
| GHF analytic gradient | [ghf-grad-h3-sto3g](methods/ghf-grad-h3-sto3g.md) | `…/ghf-grad-h3-sto3g` |
| C-PCM free-energy forces | [pcm-grad-h2o](methods/pcm-grad-h2o.md) | `…/pcm-grad-h2o` |
| SCS-MP2 | [scs-mp2-h2o-sto3g](methods/scs-mp2-h2o-sto3g.md) | `…/scs-mp2-h2o-sto3g` |
| Full LR-TDDFT (RKS) | [full-lr-tddft-h2-pbe](methods/full-lr-tddft-h2-pbe.md) | `…/full-lr-tddft-h2-pbe` |
| ROKS | [roks-o2-sto3g](methods/roks-o2-sto3g.md) | `…/roks-o2-sto3g` |
| RHF static polarizability | [static-polariz-h2o](methods/static-polariz-h2o.md) | `…/static-polariz-h2o` |
| IEF-PCM | [ief-pcm-h2o](methods/ief-pcm-h2o.md) | `…/ief-pcm-h2o` |
| ROMP2 | [romp2-o2-sto3g](methods/romp2-o2-sto3g.md) | `…/romp2-o2-sto3g` |
| ROKS analytic gradient | [roks-grad-o2-sto3g](methods/roks-grad-o2-sto3g.md) | `…/roks-grad-o2-sto3g` |
| Hybrid full-LR TDDFT (B3LYP RPA) | [h2-b3lyp-rpa-tddft](methods/h2-b3lyp-rpa-tddft.md) | `…/h2-b3lyp-rpa-tddft` |
| Spin-flip TDA (SF-CIS) | [spin-flip-tda-o2](methods/spin-flip-tda-o2.md) | `…/spin-flip-tda-o2` |
| Spin-flip full TDDFT (SF-RPA) | [sf-tddft-o2](methods/sf-tddft-o2.md) | `…/sf-tddft-o2` |
| ROHF nuclear Hessian | [rohf-hessian-o2](methods/rohf-hessian-o2.md) | `…/rohf-hessian-o2` |
| ADC(2) | [adc2-h2co-sto3g](methods/adc2-h2co-sto3g.md) | `…/adc2-h2co-sto3g` |
| ADC(2)-x | [adc2x-h2co-sto3g](methods/adc2x-h2co-sto3g.md) | `…/adc2x-h2co-sto3g` |
| IP-ADC(2) | [ip-adc2-h2o-sto3g](methods/ip-adc2-h2o-sto3g.md) | `…/ip-adc2-h2o-sto3g` |
| EA-ADC(2) | [ea-adc2-h2-sto3g](methods/ea-adc2-h2-sto3g.md) | `…/ea-adc2-h2-sto3g` |
| RMP3 | [mp3-h2-sto3g](methods/mp3-h2-sto3g.md) | `…/mp3-h2-sto3g` |
| RMP4(SDTQ) | [mp4-h2-sto3g](methods/mp4-h2-sto3g.md) | `…/mp4-h2-sto3g` |
| CCSD | [ccsd-h2-sto3g](methods/ccsd-h2-sto3g.md) | `…/ccsd-h2-sto3g` |
| CC2 | [cc2-h2-sto3g](methods/cc2-h2-sto3g.md) | `…/cc2-h2-sto3g` |
| CCSD(T) | [ccsd-t-h2o-sto3g](methods/ccsd-t-h2o-sto3g.md) | `…/ccsd-t-h2o-sto3g` |
| G0W0 | [g0w0-h2-sto3g](methods/g0w0-h2-sto3g.md) | `…/g0w0-h2-sto3g` |
| BSE@G0W0 | [bse-h2-sto3g](methods/bse-h2-sto3g.md) | `…/bse-h2-sto3g` |
| Direct RPA | [direct-rpa-h2-sto3g](methods/direct-rpa-h2-sto3g.md) | `…/direct-rpa-h2-sto3g` |
| One-shot GF2 (G0F2) | [gf2-h2-sto3g](methods/gf2-h2-sto3g.md) | `…/gf2-h2-sto3g` |
| B2PLYP double hybrid | [b2plyp-h2o](methods/b2plyp-h2o.md) | `…/b2plyp-h2o` |
| DFT-D4 (PBE params) | [dft-d4-h2o-pbe](methods/dft-d4-h2o-pbe.md) | `…/dft-d4-h2o-pbe` |
| CASCI | [casci-h2-sto3g](methods/casci-h2-sto3g.md) | `…/casci-h2-sto3g` |
| ΔSCF (IMOM) | [delta-scf-n2](methods/delta-scf-n2.md) | `…/delta-scf-n2` |
| EOM-EE-CCSD | [eom-ee-ccsd](methods/eom-ee-ccsd.md) | `…/eom-ee-ccsd` |
| EOM-EA-CCSD | [eom-ea-ccsd](methods/eom-ea-ccsd.md) | `…/eom-ea-ccsd` |
| EOM-SF-CCSD | [eom-sf-ccsd-o2](methods/eom-sf-ccsd-o2.md) | `…/eom-sf-ccsd-o2` |

SCF reading map: [User Guide → SCF Learning Path](../user-guide/scf-roadmap.md).
