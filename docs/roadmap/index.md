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

## Campaign principle (breadth-first)

As of **2026-07-23**, active development is on **`dev`**. Architecture work that
recently closed: **backend-mirror** (layout rewrite) + **cuda-resurrect**
(CUDA builds and mostly passes on GH200). The CPU breadth campaign
(frame-block + cpu-breadth-01…07) already proved the Driver / Component /
Frame spine across Evaluate, MD, Minimize, Hessian, MC, Path/NEB, REMD, and
umbrella bias.

**Principle: CPU-only method breadth for framework generality.**

The goal is **not** product depth in one family and **not** CUDA ports.
It is to show that one Driver / Component / Frame spine can host
**methods from electronic structure through mesoscale**, each as a small
textbook case (`regressions/`, literature/analytic oracle, ≤ 60 s).

Prefer a **new method family or scale** over another variant of an already-proven
shape (second soft-core style, second isotropic cell tweak, …). Depth
(Jensen/PySCF parity, multi-window WHAM product, ReaxFF, GAFF polymers,
metadynamics suites) waits until the scale ladder has more rungs.

| Scale | Already framework-proven (do not re-prove) |
| ----- | ------------------------------------------ |
| Electronic (mean-field) | RHF/UHF/ROHF/GHF, RKS/UKS, AM1/RM1, PW, GPW; FD Hessian; RHF TDA/TDHF |
| Electronic (new families) | `RKS::tda` CIS@RKS skeleton; H-only SCC-DFTB2; PW stress → isotropic cell |
| Classical atomistic | LJ+bonded, FENE, Coul/PME, EAM/SW; MD NVE/NVT/NPT; FIRE/LBFGS/CG; RATTLE/SETTLE |
| Soft / alchemical / meso | soft LJ/Coul, WCA, Yukawa, DPD; λ TI, Widom; H-REMD λ |
| Paths / replicas | NEB, dimer, string, IRC, T-REMD, H-REMD, PIMD+Rg², IdentityReplication |
| Sampling | MC NVT/NPT/μVT; Distance/Angle/**Dihedral** CV + umbrella; Brownian |
| Cell | isotropic + joint atom+cell (classical virial **and** PW electronic stress) |

| Scale | **Not** yet framework-proven (good next methods) |
| ----- | ------------------------------------------------ |
| Electronic | real ALDA/hybrid **TDDFT**; xTB/GFN; **GPW** cell/forces product; ROHF/UKS TDA |
| Response / vib | excited-state gradients; nonadiabatic coupling; analytic Hessian / CPHF polarizability |
| Classical | anisotropic/triclinic cell DOFs; Morse production pair; multi-type FENE melt |
| Mesoscale | multi-component DPD phase; full DLVO surface; CG polymer melt entanglement |
| Cross-scale | QM/MM embedding loop; AIMD with DFTB or RKS force family; metadynamics CV bias |

### Closed CPU breadth queues (2026-07-23)

| Queue | Items |
| ----- | ----- |
| breadth-first 1–6 | FIRE, soft LJ, Widom, PIMD, isotropic cell, dimer |
| breadth-next 1–9 | soft Coulomb, AngleCV, FEP/TI, joint cell, string, IRC, RPMD Rg², SETTLE, Brownian |
| **method ladder 1–8** | RKS-TDA, SCC-DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD λ |

### Recommended next queue — **CPU open scales** (still method-first)

**CPU only. No CUDA. Prefer a scale not already proven.** Each row opens a
**new framework question**, not product depth.

| # | Method | Scale | Framework question | Textbook gate (sketch) | Out of scope |
| - | ------ | ----- | ------------------ | ---------------------- | ------------ |
| 1 | **ALDA RKS-TDDFT TDA** | electronic / response | Can \(f^{\mathrm{xc}}\) sit in the same A-matrix path as CIS@RKS? | H₂ / H₂O first singlet vs PySCF TDDFT-TDA | Forces, hybrids |
| 2 | **xTB / GFN0-1 skeleton** | semi-empirical TB | Second TB family beyond DFTB2 toy tables | Minimal organics energy vs offline | Full GFN2 product |
| 3 | **Anisotropic / triclinic cell DOFs** | cell | Cell matrix beyond isotropic scale | Sheared LJ box pressure residual | Full PR barostat |
| 4 | **QM/MM mechanical embedding** | multiscale | ES force + classical region on one MD/Minimize | Small solute + LJ shell energy | Polarizable embedding |
| 5 | **Metadynamics (1D well-tempered)** | enhanced sampling | History-dependent Bias on existing CV protocol | 1D double-well reweight sketch | Multi-CV OPES product |
| 6 | **DPD multi-component / soft polymer** | mesoscale product | Composition or chain on DPD spine | A/B mix density or R_g of FENE-DPD chain | Phase diagram campaign |
| 7 | **Excited-state forces (TDA)** | response | Post-SCF gradient on Frame for surface hopping prep | Analytic or FD force on S1 for H₂ | Full NAC / TSH suite |
| 8 | **RI-MP2 energy skeleton** | post-HF | Correlated method as Evaluate force/energy kind | H₂O correlation energy pin | Gradients, triples |

Order is suggested. Skip rows that need large data tables; prefer any open
**scale** over deepening an already-proven one.

**Never on this track (until a consumer asks)**

- Any CUDA / dual-backend port
- HF/KS product depth (full PCM family, production analytic Hessians…)
- Multi-window WHAM service; ReaxFF / ML potentials / full GAFF polymer protocols

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

- **RHF TDA (CIS) and TDHF/RPA** ship as post-SCF members (`RHF::tda` /
  `RHF::tdhf`) with CH₂O regression. **`RKS::tda`** is CIS@RKS (H₂ gate);
  true ALDA/hybrid TDDFT, excited-state gradients, and nonadiabatic
  couplings are still absent — surface hopping stays blocked.

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
- SCC-DFTB2 H-only skeleton ships (`DFTB2`); xTB / full SK tables / DFTB3
  remain ○ (still not under `NDDOSolver`).

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
- **Constraints:** RATTLE (CPU, MPI) holds TIP3P geometry at `dt = 2 fs`.
  No SETTLE/LINCS/rigid bodies/virtual sites.
- **Torch backend** is a differentiable classical FF (LJ + harmonic bond),
  not an ML potential and not a dynamics engine.

### Molecular dynamics

- NVE / NVT / NPT on CPU and CUDA. Standalone `NoseHoover` is a **single**
  chain; M-chain NHC lives inside MTK NPT only. Bussi/CSVR is CPU-only.
- MTK couples the **cell diagonal** only → no PR / shear / constant-stress.
- **NPH unreachable** (thermostat fused into MTK barostat).
- `Deform` uniaxial only. CUDA reciprocal space is single-GPU; CPU PME
  replicates the grid rather than distributing the FFT.
- **CUDA PME Madelung residual** closed 2026-07-23 (`cuda-pme-accuracy`).

### Sampling and paths

- MC: NVT / NPT / μVT GCMC on CPU (single rank; O(N) trial energy by design).
  No CBMC, Gibbs. **Widom** ships as snapshot Observer (`cpu::Widom`,
  `reg_widom_hs`).
- **Umbrella sampling exists** (`bias::DistanceCV` + `HarmonicBias`,
  `reg_umbrella_harmonic`). Multi-window WHAM production path still thin
  (molrs-owned analysis incomplete). FEP/TI/BAR absent.
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
