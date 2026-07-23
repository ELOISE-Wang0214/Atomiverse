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

**Principle going forward: breadth before depth.** Prefer simple, diverse
methods that exercise a *new* combination of (workflow × backend × physics
axis) over productizing one method family. Depth (GAFF polymers, full
Jensen/PySCF HF parity, multi-window WHAM campaigns, CUDA SCF gradients) is
parked until the framework has been poked from more angles.

What already counts as **framework-proven** (do not re-prove with another
variant of the same shape):

| Axis | Proven on |
| ---- | --------- |
| Loop shapes | Evaluate, MD, Minimize, Hessian, MC (NVT/NPT/μVT), Path/NEB, REMD, single-window umbrella |
| Serial + MPI DD | classical MD / NPT-PME / RATTLE |
| Energy families | LJ+bonded, CoulCut/Long+PME, EAM/SW, RHF/UHF/ROHF/GHF, RKS/UKS, AM1/RM1, PW, GPW |
| Dual backend (MD/Minimize/Evaluate classical) | CUDA resurrected; residual **cuda-pme-accuracy** debt |
| Replicas | `ReplicaSet` + expansions; NEB coupled band; REMD N owned MD programs; IdentityReplication seed |

What is **not** yet framework-proven (good breadth targets — see next section):

- Cell DOFs in Minimize (variable-cell / joint atom+cell)
- Nuclear quantum beads (PIMD / RPMD)
- Free-energy *estimation* beyond Widom snapshot + single biased window (FEP/TI one-leg, multi-window WHAM)
- A second transition-state algorithm (dimer / string) — NEB alone is not enough
- Soft-core / λ-dependent potentials (alchemical parameterization)
- Sampling / path / bias / constraints on CUDA
- A non-MD Minimizer algorithm beyond L-BFGS/CG (FIRE)

### Recommended next queue (CPU breadth-first)

**This track is CPU-only.** Dual-backend / CUDA residuals
(`cuda-pme-accuracy`, CUDA PairCoulCut, CUDA SCF gradients, CUDA MC/path/bias
mirrors) are a **separate ownership lane** — they do not appear in the
numbered queue below and are not assigned to the CPU breadth campaign.

Ordered for **cheap diversity** that stress-tests the API, not for production
completeness. Each item should ship one textbook regression under
`regressions/` with a literature/analytic oracle. **CPU only.**

| # | Task | Why it is breadth | Stay out of scope |
| - | ---- | ----------------- | ----------------- |
| 1 | ~~**FIRE minimizer**~~ | **SHIPPED 2026-07-23** (`cpu::minimizer::FIRE`, `reg_fire_lj_dimer`) | — |
| 2 | ~~**Soft-core LJ (λ)**~~ | **SHIPPED 2026-07-23** (`cpu::PairLJSoft`, `reg_lj_soft_dimer`) | — |
| 3 | ~~**Widom insertion**~~ | **SHIPPED 2026-07-23** (`cpu::Widom`, `reg_widom_hs`) | — |
| 4 | **PIMD skeleton** (harmonic or free particle first) | Second real consumer of `IdentityReplication` + bead coupling; nuclear-quantum axis | RPMD estimators, ring-polymer contraction, multi-GPU |
| 5 | **Variable-cell Minimize** (classical virial or PW strain-FD stress) | Cell DOFs in the optimizer — never exercised | Fully anisotropic PR cell dynamics, joint cell+atom fancy coords |
| 6 | **Dimer TS** (or growing string) | Second path method; proves PathDriver is not NEB-only | IRC, conical intersections, CUDA path |

**Explicitly deprioritized right now** (depth / productization / gated backlog):

- Entire **CUDA ownership lane** (PME accuracy residual, device coulomb cut,
  CUDA SCF gradients, sampling/path/bias ports) — not this track
- `hf-parity` chain past hard-gates of *already-shipped* code (CFMM-J / C-PCM /
  x2c / embed docs honesty) — only promote when a concrete consumer needs it
- Polarizability, full UHF/KS analytic Hessian, TDDFT gradients, RI-MP2
  (`hf-parity-10+`)
- Multi-window umbrella + production molrs WHAM campaign
- ReaxFF / AIREBO / ML potentials / full GAFF polymer protocols
- Metadynamics / OPES / ABF (wait until CV layer has ≥2 CVs + free-energy path)

Active parking lots (do not `/mol:impl` the umbrella): `hf-deferred-backlog`,
`hf-parity-*` drafts. CUDA residuals live under their own specs
(`cuda-pme-accuracy`, …) and are out of this queue.

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
  `RHF::tdhf`) with CH₂O regression. No TDDFT, no excited-state gradients,
  no nonadiabatic couplings — so surface hopping stays blocked.

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
- No DFTB / xTB (separate method families when they land — not `NDDOSolver`).

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
- **CUDA residual:** PME reciprocal energy accuracy debt
  (`cuda-pme-accuracy` spec) — classical MD path otherwise green after
  cuda-resurrect.

### Sampling and paths

- MC: NVT / NPT / μVT GCMC on CPU (single rank; O(N) trial energy by design).
  No CBMC, Widom, Gibbs.
- **Umbrella sampling exists** (`bias::DistanceCV` + `HarmonicBias`,
  `reg_umbrella_harmonic`). Multi-window WHAM production path still thin
  (molrs-owned analysis incomplete).
- **REMD exists** (`cpu::path::Remd`, per-replica owned MD programs,
  `reg_remd_harmonic`). Hamiltonian RE / REST2 / multidimensional RE absent.
- **CI-NEB runs** (`PathDriver` + `NebCoupling` + `BandLBFGS`). No dimer,
  string, IRC. PIMD not yet a second IdentityReplication consumer.

### Optimization

- CPU: L-BFGS + Polak–Ribière⁺ CG. CUDA: those + steepest descent. **FIRE
  is ○** (good breadth candidate).
- **Atoms only** — no cell / variable-cell / joint optimization (good
  breadth candidate once stress/virial is wired into Minimize).
- Cartesian only; no internal coordinates.
