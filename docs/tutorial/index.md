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

### MetalWalls capability wave 1 (classical electrochem + LAMMPS-overlap MD)

Offline LAMMPS / MetalWalls / literature pins only — see
[metalwalls-method-handoff](../developer/metalwalls-method-handoff.md) and
`tools/metalwalls_oracle/`. Narrative: [Ch. 10 Classical MD](ch10-classical-md.md).
Oracle priority: **LAMMPS** when the feature overlaps; MetalWalls for MW-only.

| Method | Page | Regression | Oracle |
|--------|------|------------|--------|
| Const-pot BO parallel plate (madelung, E1) | [const-pot-bo-parallel-plate](methods/const-pot-bo-parallel-plate.md) | `regressions/md/const-pot-bo-parallel-plate` | lammps (`perf_exempt`) |
| Freeze electrode atoms (`setforce 0`) | [freeze-electrode-atoms](methods/freeze-electrode-atoms.md) | `regressions/md/freeze-electrode-atoms` | lammps (ratio **3.98**) |
| Electrode total-charge observer | [const-pot-charge-observer](methods/const-pot-charge-observer.md) | `regressions/md/const-pot-charge-observer` | lammps (`perf_exempt`) |
| Charging work \(-\sum_i \Psi_i Q_i\) | [charging-work-energy](methods/charging-work-energy.md) | `regressions/md/charging-work-energy` | literature (`perf_exempt`) |
| Const-pot MD NVT short trajectory | [const-pot-md-nvt-short](methods/const-pot-md-nvt-short.md) | `regressions/md/const-pot-md-nvt-short` | metalwalls (`perf_exempt`) |
| Gaussian electrode charge kernels (E2) | [gaussian-electrode-charge](methods/gaussian-electrode-charge.md) | `regressions/md/gaussian-electrode-charge` | literature (`perf_exempt`) |
| Slab 2D Ewald foothold (Yeh–Berkowitz) | [slab-2d-ewald-foothold](methods/slab-2d-ewald-foothold.md) | `regressions/md/slab-2d-ewald-foothold` | lammps (ratio **6.49**) |
| External uniform E-field MD | [external-efield-md](methods/external-efield-md.md) | `regressions/md/external-efield-md` | lammps (ratio **0.60**) |
| Buckingham pair (LAMMPS `buck`) | [buckingham-pair-lammps](methods/buckingham-pair-lammps.md) | `regressions/forcefield/buckingham-pair-lammps` | lammps (`perf_exempt`) |
| Four-site TIP4P foothold | [four-site-tip4p-foothold](methods/four-site-tip4p-foothold.md) | `regressions/md/four-site-tip4p-foothold` | lammps (ratio **0.338**) |

### MetalWalls capability wave 2 (classical electrochem + LAMMPS-overlap MD)

Offline LAMMPS / MetalWalls / literature pins only — see
[metalwalls-method-handoff](../developer/metalwalls-method-handoff.md) and
`tools/metalwalls_oracle/`. Narrative: [Ch. 10 Classical MD](ch10-classical-md.md).
Oracle priority: **LAMMPS** when the feature overlaps; MetalWalls for MW-only.
Full index: [methods/README](methods/README.md).

| Method | Page | Regression | Oracle |
|--------|------|------------|--------|
| Const-pot CG large multi-site electrode (E2) | [const-pot-cg-large-electrode](methods/const-pot-cg-large-electrode.md) | `regressions/md/const-pot-cg-large-electrode` | metalwalls (`perf_exempt`) |
| Electrode charge predictor (Kolafa ASPC, E2) | [electrode-charge-predictor](methods/electrode-charge-predictor.md) | `regressions/md/electrode-charge-predictor` | metalwalls (`perf_exempt`) |
| Finite-field + single electrode (Dufils 3D, E3) | [finite-field-single-electrode](methods/finite-field-single-electrode.md) | `regressions/md/finite-field-single-electrode` | metalwalls (`perf_exempt`) |
| Bulk PIM self-consistent induced dipoles (E4) | [pim-induced-dipole-bulk](methods/pim-induced-dipole-bulk.md) | `regressions/md/pim-induced-dipole-bulk` | metalwalls (`perf_exempt`) |
| Rigid piston electrode constant \(P_z\) (E6) | [piston-electrode-pz](methods/piston-electrode-pz.md) | `regressions/md/piston-electrode-pz` | metalwalls (`perf_exempt`) |
| Thomas–Fermi finite metallicity electrode (E7) | [thomas-fermi-electrode](methods/thomas-fermi-electrode.md) | `regressions/md/thomas-fermi-electrode` | metalwalls (`perf_exempt`) |
| Born–Mayer–Huggins pair foothold (E8) | [born-repulsion-foothold](methods/born-repulsion-foothold.md) | `regressions/forcefield/born-repulsion-foothold` | lammps (`perf_exempt`) |
| Fumi–Tosi multi-type Born + Coulomb NaCl (E8) | [fumi-tosi-nacl](methods/fumi-tosi-nacl.md) | `regressions/forcefield/fumi-tosi-nacl` | lammps (`perf_exempt`) |
| Steele 10-4-3 wall (Ar single-point, E8) | [steele-wall-argon](methods/steele-wall-argon.md) | `regressions/forcefield/steele-wall-argon` | lammps (`perf_exempt`) |
| Capacitance from electrode charge fluctuations (E10) | [capacitance-charge-fluctuation](methods/capacitance-charge-fluctuation.md) | `regressions/md/capacitance-charge-fluctuation` | literature (`perf_exempt`) |

### MetalWalls capability wave 3 (classical electrochem + LAMMPS-overlap MD)

Offline LAMMPS / MetalWalls / literature pins only — see
[metalwalls-method-handoff](../developer/metalwalls-method-handoff.md) and
`tools/metalwalls_oracle/`. Narrative: [Ch. 10 Classical MD](ch10-classical-md.md).
Oracle priority: **LAMMPS** when the feature overlaps; MetalWalls for MW-only.
Full index: [methods/README](methods/README.md).

| Method | Page | Regression | Oracle |
|--------|------|------------|--------|
| Constant-force ion pull (`fix addforce`) | [constant-force-pull-ion](methods/constant-force-pull-ion.md) | `regressions/md/constant-force-pull-ion` | lammps (ratio **0.12**) |
| Oscillating external E-field MD (E3) | [oscillating-efield-md](methods/oscillating-efield-md.md) | `regressions/md/oscillating-efield-md` | lammps (ratio **2.12**) |
| Slab \(\Psi(z)\) planar Poisson observer (E10) | [potential-profile-observer](methods/potential-profile-observer.md) | `regressions/md/potential-profile-observer` | literature (`perf_exempt`) |
| Surface-tension / NPγT coupling | [surface-tension-coupling](methods/surface-tension-coupling.md) | `regressions/md/surface-tension-coupling` | lammps (ratio **2.88**) |
| Multi-topology virtual sites TIP4P+TIP5P (E9) | [virtual-site-multi-topology](methods/virtual-site-multi-topology.md) | `regressions/md/virtual-site-multi-topology` | lammps (ratio **0.084**) |
| Coupled PIM dipoles ⊕ electrode \(Q\) (E5) | [dipole-electrode-joint](methods/dipole-electrode-joint.md) | `regressions/md/dipole-electrode-joint` | metalwalls (`perf_exempt`) |

### CP2K capability wave 1 (GPW / PW / periodic)

Offline CP2K pins only — see [cp2k-method-handoff](../developer/cp2k-method-handoff.md)
and `tools/cp2k_oracle/`. Narrative chapters: [Ch. 8 PW](ch08-plane-wave.md),
[Ch. 9 GPW](ch09-gpw.md).

| Method | Page | Regression |
|--------|------|------------|
| GPW Γ-RKS H₂O PBE (MT) | [gpw-water-pbe-mt](methods/gpw-water-pbe-mt.md) | `regressions/gpw/gpw-water-pbe-mt` |
| GPW analytic forces H₂O | [gpw-forces-water](methods/gpw-forces-water.md) | `regressions/gpw/gpw-forces-water` |
| Γ-UKS GPW O₂ | [gpw-uks-o2-gamma](methods/gpw-uks-o2-gamma.md) | `regressions/gpw/gpw-uks-o2-gamma` |
| GPW PBE-D3(BJ) ice dimer | [gpw-pbe-d3bj-ice](methods/gpw-pbe-d3bj-ice.md) | `regressions/gpw/gpw-pbe-d3bj-ice` |
| PW LDA Si (k-mesh) | [pw-lda-si](methods/pw-lda-si.md) | `regressions/pw/pw-lda-si` |
| PW true PBE GGA Si | [pw-pbe-gga-si](methods/pw-pbe-gga-si.md) | `regressions/pw/pw-pbe-gga-si` |
| Metallic Fermi–Dirac Al | [pw-smearing-al](methods/pw-smearing-al.md) | `regressions/pw/pw-smearing-al` |
| PW collinear spin O₂ | [pw-spin-polarized-o2](methods/pw-spin-polarized-o2.md) | `regressions/pw/pw-spin-polarized-o2` |
| PW BOMD H (NVE foothold) | [pw-bomd-h](methods/pw-bomd-h.md) | `regressions/pw/pw-bomd-h` |
| OT SCF water (GPW PBE) | [ot-minimizer-water](methods/ot-minimizer-water.md) | `regressions/gpw/ot-minimizer-water` |

### CP2K capability wave 2 (GPW / PW / periodic)

Offline CP2K pins only — see [cp2k-method-handoff](../developer/cp2k-method-handoff.md)
and `tools/cp2k_oracle/`. Narrative chapters: [Ch. 8 PW](ch08-plane-wave.md),
[Ch. 9 GPW](ch09-gpw.md). Full index: [methods/README](methods/README.md).

| Method | Page | Regression |
|--------|------|------------|
| GPW cell stress strain-FD H₂O (NUMERICAL pin) | [gpw-stress-water](methods/gpw-stress-water.md) | `regressions/gpw/gpw-stress-water` |
| GPW BOMD H₂O NVE energy-drift foothold | [gpw-bomd-water](methods/gpw-bomd-water.md) | `regressions/gpw/gpw-bomd-water` |
| PW Nielsen–Martin analytic stress Si | [pw-analytic-stress-si](methods/pw-analytic-stress-si.md) | `regressions/pw/pw-analytic-stress-si` |
| FD supercell phonons diamond Si (Γ) | [fd-phonon-supercell](methods/fd-phonon-supercell.md) | `regressions/pw/fd-phonon-supercell` |
| Dudarev DFT+U NiO mol (Γ-UKS GPW) | [gpw-dftu-nio](methods/gpw-dftu-nio.md) | `regressions/gpw/gpw-dftu-nio` |
| PW bands + DOS bulk Si | [bands-dos-si](methods/bands-dos-si.md) | `regressions/pw/bands-dos-si` |
| PW k-point collinear UKS O₂ | [pw-uks-kmesh-fe](methods/pw-uks-kmesh-fe.md) | `regressions/pw/pw-uks-kmesh-fe` |
| GPW Γ-RKS PBEsol Si₂ (MT) | [gpw-pbesol-si](methods/gpw-pbesol-si.md) | `regressions/gpw/gpw-pbesol-si` |
| GPW PBE-D4 periodic water dimer | [gpw-d4-water](methods/gpw-d4-water.md) | `regressions/gpw/gpw-d4-water` |

### PySCF campaign (molecular post-HF / response)

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
| CCSD nuclear gradient | [ccsd-grad-h2o-sto3g](methods/ccsd-grad-h2o-sto3g.md) | `…/ccsd-grad-h2o-sto3g` |
| CC2 | [cc2-h2-sto3g](methods/cc2-h2-sto3g.md) | `…/cc2-h2-sto3g` |
| CCSD(T) | [ccsd-t-h2o-sto3g](methods/ccsd-t-h2o-sto3g.md) | `…/ccsd-t-h2o-sto3g` |
| G0W0 | [g0w0-h2-sto3g](methods/g0w0-h2-sto3g.md) | `…/g0w0-h2-sto3g` |
| BSE@G0W0 | [bse-h2-sto3g](methods/bse-h2-sto3g.md) | `…/bse-h2-sto3g` |
| Direct RPA | [direct-rpa-h2-sto3g](methods/direct-rpa-h2-sto3g.md) | `…/direct-rpa-h2-sto3g` |
| One-shot GF2 (G0F2) | [gf2-h2-sto3g](methods/gf2-h2-sto3g.md) | `…/gf2-h2-sto3g` |
| B2PLYP double hybrid | [b2plyp-h2o](methods/b2plyp-h2o.md) | `…/b2plyp-h2o` |
| DFT-D4 (PBE params) | [dft-d4-h2o-pbe](methods/dft-d4-h2o-pbe.md) | `…/dft-d4-h2o-pbe` |
| CASCI | [casci-h2-sto3g](methods/casci-h2-sto3g.md) | `…/casci-h2-sto3g` |
| CASSCF | [casscf-h2-sto3g](methods/casscf-h2-sto3g.md) | `…/casscf-h2-sto3g` |
| SC-NEVPT2 | [nevpt2-h2-sto3g](methods/nevpt2-h2-sto3g.md) | `…/nevpt2-h2-sto3g` |
| ΔSCF (IMOM) | [delta-scf-n2](methods/delta-scf-n2.md) | `…/delta-scf-n2` |
| EOM-EE-CCSD | [eom-ee-ccsd](methods/eom-ee-ccsd.md) | `…/eom-ee-ccsd` |
| EOM-EA-CCSD | [eom-ea-ccsd](methods/eom-ea-ccsd.md) | `…/eom-ea-ccsd` |
| EOM-IP-CCSD | [eom-ip-ccsd](methods/eom-ip-ccsd.md) | `…/eom-ip-ccsd` |
| EOM-SF-CCSD | [eom-sf-ccsd-o2](methods/eom-sf-ccsd-o2.md) | `…/eom-sf-ccsd-o2` |
| EOM-CC2 | [eom-cc2-h2](methods/eom-cc2-h2.md) | `…/eom-cc2-h2` |
| RKS nuclear Hessian (FD of forces) | [rks-hessian-h2o-sto3g](methods/rks-hessian-h2o-sto3g.md) | `…/rks-hessian-h2o-sto3g` |
| UKS nuclear Hessian (FD of forces) | [uks-hessian-o2-sto3g](methods/uks-hessian-o2-sto3g.md) | `…/uks-hessian-o2-sto3g` |
| UKS full-LR TDDFT (LSDA O₂) | [uks-full-lr-o2-lsda](methods/uks-full-lr-o2-lsda.md) | `…/uks-full-lr-o2-lsda` |
| UMP2 analytic gradient | [ump2-grad-o2-sto3g](methods/ump2-grad-o2-sto3g.md) | `…/ump2-grad-o2-sto3g` |

SCF reading map: [User Guide → SCF Learning Path](../user-guide/scf-roadmap.md).
