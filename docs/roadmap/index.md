# Scientific Methods Roadmap

This roadmap is the honest inventory of what Atomiverse can actually compute
today, method by method, on each execution backend. It is written against the
code on `dev`, not against intent: a row is only marked implemented if the
method is reachable from the public API and validated.

## Status

| Symbol | Meaning                                 |
| ------ | --------------------------------------- |
| ○      | Not implemented                         |
| ✓      | Implemented and numerically correct     |
| ⚡      | Numerically correct and fully optimized |

A CUDA implementation must execute its principal computation on CUDA. A CPU
fallback is not considered CUDA support.

### What this roadmap does not list

A method that a **strictly better algorithm dominates** — same purpose, same or lower
cost, no regime where the older one wins — is simply **not in these tables**. It is not
marked, not tracked, and not planned. Listing it, even as "will not implement", would
suggest the roadmap has an opinion worth reading about it; it does not.

So there is no BFGS row (L-BFGS matches its convergence at O(mN) memory instead of
O(N²)), no Fletcher–Reeves row (Polak–Ribière⁺ restarts strictly better at the same
cost), no Berendsen-thermostat row (CSVR samples the canonical ensemble; Berendsen does
not sample it at all).

This is only about *dominance*, never about our schedule. A method that survives because
it occupies a genuinely different cost/accuracy niche stays ○ and stays listed — steepest
descent is dominated by CG for minimization proper, but it is the robust first resort on a
badly broken start geometry, so it keeps its row. Brownian dynamics keeps its row too: it
is different physics from Langevin, not a worse integrator.

### How the symbols are assigned here

**✓** means the method runs from the public API and its numbers are pinned
against an external reference or an analytic result — a validation case under
`benchmarks/<domain>/<case>/validation/` (RHF, MP2, FCI, GPW, vibrations,
plane-wave Si, LJ MD, bonded MD vs LAMMPS, LJ MC vs NIST), or a unit test that
checks the physics rather than the plumbing.

**⚡** is deliberately scarce. It is reserved for kernels that are both correct
*and* carry a committed performance suite measured against a production engine.
Today that is exactly one method family: the Lennard–Jones 12-6 pair style,
benchmarked against LAMMPS/Kokkos on CPU (xsimd, `-march=native`) and on GPU
(`/* benchmarks removed; re-add later */ lj/performance/`). Everything else that is correct but has no
standing performance gate is a ✓, including code that is fast in practice.

**○** means not implemented. It does not distinguish "planned" from "not
planned"; the sections below are a method census, not a commitment schedule.

## Chapters

Each chapter is one **method family**. Periodicity, gradients, and excited-state
response are *axes*, not families — a periodic calculation is still DFT — so
those rows live inside the method that owns them: periodic DFT (plane-wave, GPW)
is a section of the DFT chapter, periodic HF a section of the Hartree–Fock one.

| Chapter | Covers |
| ------- | ------ |
| [Hartree–Fock](hartree-fock.md) | RHF/UHF/ROHF/GHF, HF gradients and Hessians, periodic HF, SCF infrastructure |
| [Density-Functional Theory](dft.md) | Kohn–Sham methods, functionals (LDA → double hybrid), dispersion, **periodic DFT: plane-wave, GPW, forces, stress, phonons** |
| [Post-Hartree–Fock](post-hf.md) | Møller–Plesset, configuration interaction, coupled cluster, multireference |
| [Excited States and Many-Body](excited-states.md) | CIS/TDHF/TDDFT, ADC, EOM-CC, GW/BSE/RPA |
| [Semiempirical and Tight Binding](semiempirical.md) | NDDO/AM1, DFTB, extended tight binding (xTB) |
| [Force Fields](force-fields.md) | Bonded, nonbonded, electrostatics, polarizable, many-body/reactive, coarse-grained, machine-learning potentials |
| [Molecular Dynamics](molecular-dynamics.md) | Integrators, thermostats, barostats, ensembles, constraints, nonequilibrium MD, AIMD, nuclear quantum, nonadiabatic, vibrations |
| [Geometry Optimization](optimization.md) | Minimizers, transition states, reaction paths |
| [Sampling and Free Energy](sampling.md) | Monte Carlo, enhanced sampling, free-energy methods |
| [Multiscale and Properties](multiscale.md) | Atomic charges, implicit solvent, QM/MM, embedding, relativistic, intermolecular interactions |

## Campaign principle (CPU method-first)

As of **2026-07-23** on **`dev`**. Closed architecture: **backend-mirror**,
**cuda-resurrect**, **cuda-pme-accuracy**. Active track is **CPU-only method
breadth** — textbook cases that prove Driver / Component / Frame can host
**electronic structure through mesoscale**, not product depth in one family
and not CUDA ports.

Each shipped method is a public API + `regressions/` gate (≤ 60 s, literature
or analytic oracle). Prefer a **new family or scale** over another variant of
a proven shape.

### Already framework-proven

| Scale | Methods (do not re-prove) |
| ----- | ------------------------- |
| Electronic mean-field | RHF/UHF/ROHF/GHF, RKS/UKS, AM1/RM1, PW, GPW; FD Hessian |
| Electronic TB / response | RHF TDA/TDHF; `RKS::tda` (+ LSDA ALDA diagonal); FD TDA forces; NAC strength; H-only DFTB2; GFN0 toy; public `RHF::mp2` |
| Classical FF / MD | LJ+bonded, FENE, WCA, Coul/PME, EAM/SW; NVE/NVT/NPT; FIRE/LBFGS/CG; RATTLE/SETTLE; Brownian |
| Soft / meso | soft LJ/Coul, Yukawa, DPD (+ multi-type \(a_{ij}\)); λ TI; Widom; H-REMD λ |
| Paths / replicas | NEB, dimer, string, IRC, T-REMD, H-REMD, PIMD+Rg² |
| Sampling | MC NVT/NPT/μVT; Distance/Angle/Dihedral CV + umbrella; 1D well-tempered metaD |
| Cell | isotropic + joint + anisotropic + triclinic stress relax; PW electronic stress |
| Multiscale | mechanical + electrostatic QM/MM (`QmmmMechanical`, `QmmmElectrostatic`) |

### Closed queues (2026-07-23)

| Queue | Items |
| ----- | ----- |
| breadth-first 1–6 | FIRE, soft LJ, Widom, PIMD, isotropic cell, dimer |
| breadth-next 1–9 | soft Coulomb, AngleCV, FEP/TI, joint cell, string, IRC, RPMD Rg², SETTLE, Brownian |
| method ladder 1–8 | RKS-TDA, SCC-DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD λ |
| open scales 1–8 | ALDA flag, GFN0, aniso cell, mech. QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| depth-next 1–4 | ALDA diagonal, ES QM/MM, TDA NAC strength, triclinic cell |

### Recommended next (still method-first, CPU)

| Theme | Candidates |
| ----- | ---------- |
| Response | Full grid MO-product ALDA; hybrid TDDFT; analytic TDA forces; true NAC vectors |
| TB product | Real GFN / mio SK tables (not toys); DFTB3 |
| Cell / MD | Parrinello–Rahman NPT; NPH |
| Multiscale | Polarizable QM/MM; link atoms |
| Sampling | Multi-CV metaD / OPES; multi-window WHAM |
| Post-HF | RI-MP2 production path; gradients |

**Parked until a consumer asks:** CUDA ports, HF/KS product depth, ReaxFF/ML
potentials, full GAFF polymer protocols.

## Work in progress — known gaps and caveats

The tables tell you *whether* a method exists. This section tells you what to
watch out for in the ones that do. Snapshot: **2026-07-23** on `dev`.

### Gaussian SCF

- **Cartesian is the default engine space**; a **spherical (pure) path exists**
  for sets defined that way (def2-SVP spherical regression
  `regressions/scf/h2o_rhf_def2svp_sph`). Native Rys ERI covers `l ≤ 3`
  (s/p/d/f).
- **No analytic Hessian anywhere.** Frequencies / thermochemistry are
  central FD of analytic forces (`cpu::vib::Hessian`). CPHF is not
  implemented; static polarizability is absent.
- **Open shell:** UHF is the default workhorse. **ROHF** is energy + **analytic
  gradient** (O₂ vs PySCF). **GHF** is energy-only (real two-component;
  `regressions/scf/h3_ghf_sto3g`). ROKS does not exist. No ROHF-MP2.
- **Shipped infra (CPU):** DIIS, level shift, SOSCF (hard cases), stability
  analysis, RI-J/RI-JK, CFMM-J, C-PCM, point-charge embedding + external
  field forces, ECP (semi-numerical forces), x2c-1e energy (forces throw),
  counterpoise, RHF TDA/TDHF (`RHF::tda` / `RHF::tdhf`), APT/IR dipole
  derivatives. Hard gates / docs honesty for some of these are still in the
  `hf-parity` parking lot.
- **CUDA SCF is still a restricted-HF energy engine.** Device Fock/density/
  diagonalization; no DIIS, no UHF/KS, **no gradients** → no GPU AIMD.
  Device-native integrals are s-type only.
- **MINAO guess covers H, C, N, O, F only**; other elements fall back to Hcore.

### Post-Hartree–Fock

- **MP2 is closed-shell canonical only** — RHF reference, full O(N⁵) AO→MO
  transform, no frozen core, no RI-MP2 product path yet, no gradients.
  Density fitting exists for SCF J/K (`RI::{J,JK}`), not for MP2.
- **FCI is a reference tool**, not production: dense CI matrix, hard cap
  (`n_mo > 14` or `N_det > 2000`).

### Excited states

- **RHF TDA/TDHF** (`reg_ch2o_tdhf`); **`RKS::tda`** with LSDA dense-A +
  Dirac ALDA *diagonal* shift (`reg_h2_alda_vs_cis`); **FD TDA forces** and
  **NAC strength** (`RHF::tda_forces` / `tda_nac_strength`). Still absent:
  full grid MO-product ALDA, hybrid TDDFT, analytic ES forces, true NAC
  vectors — surface hopping stays blocked.

### Exchange-correlation

- Molecular selectors: `{LSDA, PBE, B3LYP, PBE0}`. B88/LYP exist only as
  B3LYP pieces — no standalone BLYP/BP86.
- libxc is an optional *test oracle* only (`ATV_BUILD_XC_ORACLE`).
- GPU XC outside plane waves: none (PW: Slater + PW92 LDA on device).

### Plane-wave DFT

- Analytic local + Ewald + KB forces; FD oracle remains.
  **Cell stress via strain-FD** of SCF energy (`pw_stress`); Nielsen–Martin
  analytic stress is a follow-on.
- Effectively LDA on the PW path (`sigma = 0` into the PBE routine). No
  spin polarization.
- Norm-conserving UPF only (no USPP/PAW). Empty UPF = toy V_loc.
- CUDA: H·ψ / FFT / Poisson / V_eff / LDA on device; Ritz + SCF control on host.

### GPW

- Γ-point, closed-shell only. No k-points, spin, stress, or hybrids.
- Analytic forces available → can drive AIMD; keep
  `GpwConfig.local_pp_on_grid` consistent between energy and gradient.
- DFT-D3 zero-damping + PBE params only.

### Semiempirical

- **AM1 and RM1** on the shared NDDO path with analytic forces → Evaluate /
  MD / Minimize. Hessian remains FD. AM1-BCC charges via molrs typing.
- SCC-DFTB2 H-only (`DFTB2`) and GFN0-xTB toy (`XtbGfn0`) ship as separate
  TB families (not under `NDDOSolver`). Full mio/GFN tables / DFTB3 remain ○.

### Classical force fields

- **GAFF-capable on CPU**: multi-type LJ, CoulCut, CoulLong+PME, Amber
  specials (incl. 1–4 = 0.5), harmonic bond/angle, CHARMM + multi/harmonic
  dihedrals, harmonic improper. Gaps: multi-term Fourier dihedral as a
  first-class style, `.ff` text parser, CUDA Amber parity, paper-length
  polymer protocols.
- **Many-body:** EAM + Stillinger–Weber (CPU, MPI, ghost-of-ghost halo).
  No Tersoff/MEAM/ReaxFF/CUDA many-body.
- **Open-boundary FMM tree math** is validated (`reg_fmm_point_charges`) —
  not yet a MD `KSpace` / pair-style drop-in under PBC.
- **Constraints:** RATTLE + **SETTLE** (CPU; `reg_settle_tip3p`). No LINCS /
  rigid bodies / virtual sites.
- **Torch backend** is a differentiable classical FF (LJ + harmonic bond),
  not an ML potential and not a dynamics engine.

### Molecular dynamics

- NVE / NVT / NPT on CPU and CUDA. Standalone `NoseHoover` is a **single**
  chain; M-chain NHC lives inside MTK NPT only. Bussi/CSVR is CPU-only.
- MTK couples the **cell diagonal** only → no PR NPT / constant-stress MD.
  **0 K cell optimizers:** isotropic, anisotropic, **triclinic**
  (`TriclinicCellRelax`). **NPH unreachable** (thermostat fused into MTK).
- `Deform` uniaxial only. CUDA reciprocal space is single-GPU; CPU PME
  replicates the grid rather than distributing the FFT.
- **CUDA PME Madelung residual** closed 2026-07-23 (`cuda-pme-accuracy`).

### Sampling and paths

- MC: NVT / NPT / μVT GCMC on CPU (single rank; O(N) trial energy by design).
  No CBMC, Gibbs. **Widom** (`cpu::Widom`). Soft-core **λ TI** (`dU/dλ`).
- **Umbrella:** Distance / Angle / Dihedral CV + `HarmonicBias`. **1D
  well-tempered metadynamics** (`Metadynamics1D`). Multi-window WHAM /
  multi-CV metaD / OPES still ○.
- **REMD** (T ladder) + **H-REMD** (λ ladder). NEB / dimer / string / IRC /
  PIMD skeletons ship.
- **REMD exists** (`cpu::path::Remd`, per-replica owned MD programs,
  `reg_remd_harmonic`). Hamiltonian RE / REST2 / multidimensional RE absent.
- **CI-NEB** (`PathDriver` + `NebCoupling` + `BandLBFGS`) and **dimer TS**
  (`cpu::path::Dimer`, `reg_dimer_saddle`) ship. No string / IRC.
- **PIMD skeleton** ships (`cpu::path::Pimd` + ring springs +
  `IdentityReplication`, `reg_pimd_free_particle`). No RPMD estimators /
  ring-polymer contraction / PILE.

### Optimization

- CPU: L-BFGS + Polak–Ribière⁺ CG + **FIRE** (`cpu::minimizer::FIRE`).
  CUDA: L-BFGS / CG / steepest descent (no FIRE).
- **Isotropic variable-cell** ships (`IsotropicCellRelax` +
  `PressureConverged`). No anisotropic cell DOFs, no joint atom+cell.
- Cartesian only; no internal coordinates.
