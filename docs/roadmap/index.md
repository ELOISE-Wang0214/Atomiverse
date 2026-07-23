# Scientific Methods Roadmap

This roadmap is the honest inventory of what Atomiverse can actually compute
today, method by method, on each execution backend. It is written against the
code on `dev`, not against intent: a row is only marked implemented if the
method is reachable from the public API and validated.

## Backends

| Backend | Role |
| ------- | ---- |
| **CPU** | Primary portable path (OpenMP / MPI; xsimd where used). |
| **CUDA** | Device path for selected classical MD and SCF/PW kernels (GH200). |
| **Kokkos** | Planned portable performance backend (host + device via Kokkos); **○ everywhere today**. |

Status tables use three columns: `CPU | CUDA | Kokkos`. A CUDA or Kokkos mark
requires the **principal** work of that method to run on that backend — a CPU
fallback is not counted as CUDA/Kokkos support.

## Status

| Symbol | Meaning                                 |
| ------ | --------------------------------------- |
| ○      | Not implemented                         |
| ✓      | Implemented and numerically correct     |
| ⚡      | Numerically correct and fully optimized |

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
Today that is exactly one method family: the Lennard–Jones 12-6 pair style on
CPU (xsimd, `-march=native`) and CUDA, historically benchmarked against
LAMMPS/Kokkos (`/* benchmarks removed; re-add later */ lj/performance/`).
Everything else that is correct but has no standing performance gate is a ✓,
including code that is fast in practice. No Kokkos ⚡ yet (backend ○).

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
**cuda-resurrect**, **cuda-pme-accuracy**. Closed method campaign: **CPU next
A/B/C** (parallel worktrees `feat/cpu-next-{md,tb-arch,msamp,response}`,
merged and deleted). Active track remains **CPU-only method breadth** —
textbook cases that prove Driver / Component / Frame can host **electronic
structure through mesoscale**, not product depth in one family and not CUDA
ports.

Each shipped method is a public API + `regressions/` gate (≤ 60 s, literature
or analytic oracle). Prefer a **new family or scale** over another variant of
a proven shape.

### Already framework-proven

| Scale | Methods (do not re-prove) |
| ----- | ------------------------- |
| Electronic mean-field | RHF/UHF/ROHF/GHF, RKS/UKS, concrete AM1/RM1 / PW / GPW; FD Hessian |
| Electronic TB / response | RHF TDA/TDHF; `RKS::tda` (+ full MO-product LSDA ALDA + hybrid B3LYP/PBE0 exact-exchange kernel); analytic CIS/TDA forces; true NAC vectors; DFTB2 H/C/N/O SK; GFN0 toy; canonical + RI `RHF::mp2` |
| Classical FF / MD | LJ+bonded, FENE, WCA, Coul/PME, EAM/SW; NVE/NVT/NPT (MTK diag + **PR full-cell**); **NPH**; FIRE/LBFGS/CG; RATTLE/SETTLE; Brownian |
| Soft / meso | soft LJ/Coul, Yukawa, DPD (+ multi-type \(a_{ij}\)); λ TI; Widom; H-REMD λ |
| Paths / replicas | NEB, dimer, string, IRC, T-REMD, H-REMD, PIMD+Rg² |
| Sampling | MC NVT/NPT/μVT; Distance/Angle/Dihedral CV + umbrella; 1D + **multi-CV** well-tempered metaD; multi-window WHAM (molrs) |
| Cell | isotropic + joint + anisotropic + triclinic stress relax; PW electronic stress |
| Multiscale | mechanical + electrostatic + **polarizable** QM/MM; **H link atoms** |

### Closed queues (2026-07-23)

| Queue | Items |
| ----- | ----- |
| breadth-first 1–6 | FIRE, soft LJ, Widom, PIMD, isotropic cell, dimer |
| breadth-next 1–9 | soft Coulomb, AngleCV, FEP/TI, joint cell, string, IRC, RPMD Rg², SETTLE, Brownian |
| method ladder 1–8 | RKS-TDA, SCC-DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD λ |
| open scales 1–8 | ALDA flag, GFN0, aniso cell, mech. QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| depth-next 1–4 | ALDA diagonal, ES QM/MM, TDA NAC strength, triclinic cell |
| **cpu-next A (md)** | PR full-cell NPT, isotropic NPH |
| **cpu-next B (tb-arch)** | multi-element DFTB2 SK; concrete PW/GPW/AM1–RM1 (no public Config bags) |
| **cpu-next C (msamp)** | H link atoms, polarizable QM/MM, multi-CV metaD, multi-window WHAM |
| **cpu-next D (response)** | full MO-product ALDA, analytic TDA forces, true NAC vectors, RI-MP2 |

### Recommended next (still method-first, CPU)

| Theme | Candidates |
| ----- | ---------- |
| Response | Surface-hopping consumer of NAC; CIS(D) / ADC foothold; full RPA-TDDFT |
| TB product | Production GFN2 / 3ob tables (not toys); DFTB3; analytic DFTB forces |
| Cell / MD | Semi-isotropic / surface-tension PR; NPH anisotropic; CUDA PR/NPH parity |
| Multiscale | ONIOM / subtractive; adaptive QM/MM; long-range PME QM/MM |
| Sampling | OPES / multi-walker metaD; MBAR/PMF production curve; REST2 |
| Post-HF | MP2 gradients; ROMP2 (frozen-core ✓) |

**Parked until a consumer asks:** CUDA method ports, HF/KS product depth
(hf-parity parking lot), ReaxFF/ML potentials, full GAFF polymer protocols,
Kokkos backend (WIP archived on `archive/impl-kokkos-wip`, not on `dev`).

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
  gradient** (O₂ vs PySCF). **GHF** energy + analytic gradient (real two-component;
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

- **MP2** — closed-shell RMP2 on RHF; **UMP2** on UHF (`UHF::mp2`,
  `reg_ump2-o2-sto3g`). Canonical path: full O(N⁵) AO→MO transform.
  **RI-MP2** when `ri = RI::{J,JK}` + aux basis on RHF only
  (`reg_h2o_ri_mp2`, `MP2Result::ri`). Frozen-core via `rhf.frozen=N`
  (`reg_frozen-core-mp2-h2o`). No gradients, no ROMP2.
- **FCI is a reference tool**, not production: dense CI matrix, hard cap
  (`n_mo > 14` or `N_det > 2000`).

### Excited states

- **RHF TDA/TDHF** (`reg_ch2o_tdhf`); **`RKS::tda`** with **full MO-product
  LSDA ALDA** (`reg_h2_full_grid_alda`) and **hybrid TDDFT-TDA** (B3LYP exact-
  exchange kernel; `reg_h2-b3lyp-hybrid-tda`). Full RPA-TDDFT still ○.
- **Analytic CIS/TDA forces** (default) and **true NAC vectors**
  (`RHF::tda_forces` / `tda_nac_vectors` / `tda_nac_strength`;
  `reg_h2_analytic_tda_force`, `reg_h2_nac_vectors`). FD fallback remains.
  Surface hopping still ○ (no dynamics consumer of NAC).

### Exchange-correlation

- Molecular selectors: `{LSDA, PBE, BLYP, B3LYP, PBE0}`. BP86 still ○.
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
- Concrete method `GpwDft` (parameters as members; no public Config bag).
  Analytic forces available → can drive AIMD; keep member
  `local_pp_on_grid` consistent between energy and gradient.
- DFT-D3 zero-damping + PBE params only.

### Semiempirical

- **Concrete AM1 / RM1** (`cpu::semiempirical::AM1` / `RM1`) on the shared
  NDDO path with analytic forces → Evaluate / MD / Minimize. Hessian remains
  FD. AM1-BCC charges via molrs typing. `NDDOSolver` is internal.
- **SCC-DFTB2** multi-element H/C/N/O static SK (`DFTB2`,
  `reg_h2o_dftb2_mio`); forces are **central FD of E**, not analytic.
  **GFN0-xTB toy** (`XtbGfn0`) remains a separate TB family. Production
  GFN2 / 3ob / DFTB3 remain ○.

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
- **No machine-learning potentials** in-tree (see force-fields ML table).

### Molecular dynamics

- NVE / NVT / NPT on CPU and CUDA. Standalone `NoseHoover` is a **single**
  chain; M-chain NHC lives inside MTK NPT only. Bussi/CSVR is CPU-only.
- **Two NPT integrators on CPU:** diagonal MTK (`NoseHooverNPT`) and
  **full-cell Parrinello–Rahman** (`ParrinelloRahmanNPT`, `reg_pr_npt_lj`).
  **Isotropic NPH** ships as thermostat-free `NphIntegrator` (`reg_nph_lj`).
  **0 K cell optimizers:** isotropic, joint, anisotropic, **triclinic**.
- `Deform` uniaxial only. CUDA reciprocal space is single-GPU; CPU PME
  replicates the grid rather than distributing the FFT. CUDA PR/NPH not
  ported.
- **CUDA PME Madelung residual** closed 2026-07-23 (`cuda-pme-accuracy`).

### Sampling and paths

- MC: NVT / NPT / μVT GCMC on CPU (single rank; O(N) trial energy by design).
  No CBMC, Gibbs. **Widom** (`cpu::Widom`). Soft-core **λ TI** (`dU/dλ`).
- **Umbrella:** Distance / Angle / Dihedral CV + `HarmonicBias`. **1D
  well-tempered metadynamics** (`Metadynamics1D`) and **multi-CV**
  (`Metadynamics`, `reg_metad_2d`). **Multi-window WHAM** via molrs
  (`cpu::Wham`, `reg_wham_double_well`). OPES / multi-walker still ○.
- **REMD** (T ladder) + **H-REMD** (λ ladder). REST2 / multidimensional RE
  absent.
- **CI-NEB**, **dimer**, **string**, **IRC** ship
  (`reg_{neb,dimer,string,irc}_*`).
- **PIMD skeleton** ships (`cpu::path::Pimd` + ring springs +
  `IdentityReplication`, `reg_pimd_free_particle`). No RPMD / PILE /
  ring-polymer contraction product path.

### Optimization

- CPU: L-BFGS + Polak–Ribière⁺ CG + **FIRE** (`cpu::minimizer::FIRE`).
  CUDA: L-BFGS / CG / steepest descent (no FIRE).
- **Variable-cell:** isotropic, joint atom+cell, anisotropic, and
  triclinic (`*CellRelax` + `PressureConverged`).
- Cartesian only; no internal coordinates.
