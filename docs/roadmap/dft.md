# Density-Functional Theory

Legend and caveats: [Roadmap index](index.md).

The molecular Kohn–Sham path shares the Gaussian SCF solver with
[Hartree–Fock](hartree-fock.md): the same 2×2 config product of
{restricted, unrestricted} × {HF, KS}, the same DIIS/level-shift machinery, the
same analytic-gradient assembly.

The molecular functional set includes pure LDA/GGA selectors (`LSDA`, `PBE`,
`PBEsol`, `BLYP`), global hybrids (`B3LYP`, `PBE0`), and double hybrid
`B2PLYP` — evaluated on Becke-partitioned atomic grids with Lebedev angular
quadrature. All functionals are hand-coded and gated against libxc to 1e-12 in
tests; libxc is never a runtime dependency.

**Periodic DFT is DFT**, so it lives here too: sections 10–11 cover the
plane-wave and GPW engines, their functionals, forces/stress, phonons, and the
remaining periodic-only gaps (hybrids, meta-GGA, finite-q DFPT, multi-k MP2 /
RPA / GW depth; Gaussian k-point RKS and Γ DFPT already ship).

## 1. Kohn–Sham Methods

| Method                                    | CPU | CUDA | Kokkos |
| ----------------------------------------- | :-: | :--: | :---: |
| Restricted Kohn–Sham DFT, RKS             |  ✓  |   ○  |   ○  |
| Unrestricted Kohn–Sham DFT, UKS           |  ✓  |   ○  |   ○  |
| Restricted open-shell Kohn–Sham DFT, ROKS |  ✓  |   ○  |   ○  |
| Generalized Kohn–Sham DFT                 |  ○  |   ○  |   ○  |
| Noncollinear spin DFT                     |  ○  |   ○  |   ○  |
| Spinor Kohn–Sham DFT                      |  ○  |   ○  |   ○  |
| Finite-temperature DFT                    |  ○  |   ○  |   ○  |
| Constrained DFT                           |  ○  |   ○  |   ○  |
| Delta-SCF DFT                             |  ○  |   ○  |   ○  |
| Maximum-overlap excited-state DFT         |  ○  |   ○  |   ○  |
| Ensemble DFT                              |  ○  |   ○  |   ○  |
| DFT+U                                     |  ○  |   ○  |   ○  |
| RKS analytic gradient                     |  ✓  |   ○  |   ○  |
| UKS analytic gradient                     |  ✓  |   ○  |   ○  |
| ROKS analytic gradient                    |  ✓  |   ○  |   ○  |
| RKS analytic Hessian                      |  ✓  |   ○  |   ○  |
| UKS analytic Hessian                      |  ✓  |   ○  |   ○  |

**ROKS** (`cpu::scf::ROKS`) is energy + **analytic nuclear gradient** on the
Roothaan effective-Fock path (Guest–Saunders / PySCF canonical;
`reg_roks-o2-sto3g` energy, `reg_roks-grad-o2-sto3g` forces — O₂ triplet
LSDA/STO-3G). Forces reuse the UKS electronic gradient on split-occupation
α/β densities with the ROHF-style energy-weighted density from **KS spin
Focks** (not effective-Fock eigenvalues). KS gradients (RKS/UKS/ROKS) include
the semi-local XC contribution (LDA, GGA, and the hybrid exact-exchange
admixture) and the Becke grid-response terms. **RKS nuclear Hessian**
(`RKS::hessian`; `reg_rks-hessian-h2o-sto3g`) and **UKS nuclear Hessian**
(`UKS::hessian`; `reg_uks-hessian-o2-sto3g`) are product footholds via FD of
analytic KS forces (same style as `ROHF::hessian`), gated against PySCF FD
twins; full CPKS / UCPKS closed-form `hessian.RKS` / `hessian.UKS` is a
follow-up. There is no Kohn–Sham path on CUDA at all — the device SCF solver
is restricted HF only.

## 2. LDA Functionals

| Functional           | CPU | CUDA | Kokkos |
| -------------------- | :-: | :--: | :---: |
| Slater exchange      |  ✓  |   ✓  |   ○  |
| VWN correlation      |  ✓  |   ○  |   ○  |
| SVWN                 |  ✓  |   ○  |   ○  |
| Perdew–Wang 1992 LDA |  ✓  |   ✓  |   ○  |
| Spin-polarized LDA   |  ✓  |   ○  |   ○  |

`LSDA` = Slater (Dirac) exchange + VWN5 correlation. VWN-RPA is also
implemented, as the LDA correlation component inside B3LYP. PW92 is present as
the LDA base of PBE. The CUDA column here is **plane-wave only**: Slater + PW92
on the PW grid (`cuda/planewave/internal/xc_lda.cu`), spin-unpolarized — see
[Periodic Kohn–Sham Methods](#10-periodic-kohnsham-methods) below. There is no
XC on the device in the Gaussian path.

## 3. GGA Functionals

| Functional | CPU | CUDA | Kokkos |
| ---------- | :-: | :--: | :---: |
| PBE        |  ✓  |   ○  |   ○  |
| PBEsol     |  ✓  |   ○  |   ○  |
| revPBE     |  ○  |   ○  |   ○  |
| RPBE       |  ○  |   ○  |   ○  |
| BLYP       |  ✓  |   ○  |   ○  |
| BP86       |  ○  |   ○  |   ○  |
| B97        |  ○  |   ○  |   ○  |
| OLYP       |  ○  |   ○  |   ○  |
| HCTH       |  ○  |   ○  |   ○  |
| HCTH/407   |  ○  |   ○  |   ○  |
| SOGGA      |  ○  |   ○  |   ○  |

`XCType::BLYP` is pure B88+LYP (unit weights, no HF admixture); regression
`reg_blyp-water-sto3g` pins H2O/STO-3G vs PySCF xc=blyp. **BLYP on GPW**
ships as `XcKind::BLYP` (`reg_gpw-blyp-water`; dual-gate vs PBE). **PBEsol**
(Perdew et al. 2008; solid-state μ/β) shares `evaluate_xc` and ships on the
molecular Gaussian path and **GPW** (`XcKind::PBEsol`; `reg_gpw-pbesol-si`).

## 4. Meta-GGA Functionals

| Functional | CPU | CUDA | Kokkos |
| ---------- | :-: | :--: | :---: |
| TPSS       |  ○  |   ○  |   ○  |
| revTPSS    |  ○  |   ○  |   ○  |
| r2SCAN     |  ○  |   ○  |   ○  |
| M06-L      |  ○  |   ○  |   ○  |
| MN12-L     |  ○  |   ○  |   ○  |
| MN15-L     |  ○  |   ○  |   ○  |
| B97M-V     |  ○  |   ○  |   ○  |
| TASK       |  ○  |   ○  |   ○  |

No kinetic-energy-density (τ) channel exists in the grid engine, so the
meta-GGA rung is structurally absent, not merely unparameterized.

## 5. Global Hybrid Functionals

| Functional | CPU | CUDA | Kokkos |
| ---------- | :-: | :--: | :---: |
| B3LYP      |  ✓  |   ○  |   ○  |
| B3PW91     |  ○  |   ○  |   ○  |
| PBE0       |  ✓  |   ○  |   ○  |
| BHLYP      |  ○  |   ○  |   ○  |
| BHandHLYP  |  ○  |   ○  |   ○  |
| X3LYP      |  ○  |   ○  |   ○  |
| TPSSh      |  ○  |   ○  |   ○  |
| M06        |  ○  |   ○  |   ○  |
| M06-2X     |  ○  |   ○  |   ○  |
| M06-HF     |  ○  |   ○  |   ○  |
| BMK        |  ○  |   ○  |   ○  |
| SCAN0      |  ○  |   ○  |   ○  |
| r2SCAN0    |  ○  |   ○  |   ○  |
| MN15       |  ○  |   ○  |   ○  |

B3LYP carries a 0.20 exact-exchange fraction, PBE0 0.25. The hybrid path builds
a real exchange matrix from the ERIs, so it is available anywhere the Gaussian
SCF runs — but not in GPW or plane waves.

## 6. Range-Separated and Screened Hybrid Functionals

| Functional                                 | CPU | CUDA | Kokkos |
| ------------------------------------------ | :-: | :--: | :---: |
| CAM-B3LYP                                  |  ○  |   ○  |   ○  |
| LC-ωPBE                                    |  ○  |   ○  |   ○  |
| LC-BLYP                                    |  ○  |   ○  |   ○  |
| ωB97                                       |  ○  |   ○  |   ○  |
| ωB97X                                      |  ○  |   ○  |   ○  |
| ωB97X-D                                    |  ○  |   ○  |   ○  |
| ωB97X-D3                                   |  ○  |   ○  |   ○  |
| ωB97M-V                                    |  ○  |   ○  |   ○  |
| HSE03                                      |  ○  |   ○  |   ○  |
| HSE06                                      |  ○  |   ○  |   ○  |
| Optimally tuned range-separated hybrid DFT |  ○  |   ○  |   ○  |

Range separation would need short/long-range-partitioned ERIs (an `erf`-attenuated
integral engine); the Rys kernel has no attenuation parameter.

## 7. Double-Hybrid Functionals

| Functional                              | CPU | CUDA | Kokkos |
| --------------------------------------- | :-: | :--: | :---: |
| B2PLYP                                  |  ✓  |   ○  |   ○  |
| B2GP-PLYP                               |  ○  |   ○  |   ○  |
| PWPB95                                  |  ○  |   ○  |   ○  |
| DSD-BLYP                                |  ○  |   ○  |   ○  |
| DSD-PBEP86                              |  ○  |   ○  |   ○  |
| revDSD-PBEP86-D4                        |  ○  |   ○  |   ○  |
| PBE0-DH                                 |  ○  |   ○  |   ○  |
| PBE-QIDH                                |  ○  |   ○  |   ○  |
| ωB97X-2                                 |  ○  |   ○  |   ○  |
| Spin-component-scaled double-hybrid DFT |  ○  |   ○  |   ○  |

**B2PLYP** (`XCType::B2PLYP` + `RKS::mp2` with `hybrid_pt2_fraction` a_c=0.27)
composes hybrid KS (a_x=0.53 HF + 0.47 B88 + 0.73 LYP) with scaled RMP2 on the
KS orbitals (`reg_b2plyp-h2o` vs offline PySCF-composed pin). Other double
hybrids and SCS-DH remain ○.

## 8. Dispersion and Nonlocal Correlation

| Method                           | CPU | CUDA | Kokkos |
| -------------------------------- | :-: | :--: | :---: |
| DFT-D3 zero damping              |  ✓  |   ○  |   ○  |
| DFT-D3 Becke–Johnson damping     |  ✓  |   ○  |   ○  |
| DFT-D4                           |  ✓  |   ○  |   ○  |
| Tkatchenko–Scheffler dispersion  |  ✓  |   ○  |   ○  |
| Many-body dispersion, MBD        |  ○  |   ○  |   ○  |
| Exchange-hole dipole moment, XDM |  ○  |   ○  |   ○  |
| VV10                             |  ○  |   ○  |   ○  |
| rVV10                            |  ○  |   ○  |   ○  |
| vdW-DF                           |  ○  |   ○  |   ○  |
| vdW-DF2                          |  ○  |   ○  |   ○  |
| optB88-vdW                       |  ○  |   ○  |   ○  |

DFT-D3 zero damping (PBE: s6 = 1.0, s8 = 0.722, sr6 = 1.217) is wired into
**GPW**, from in-tree reference tables. DFT-D3 **Becke–Johnson** damping
(PBE: s6 = 1, s8 = 0.7875, a1 = 0.4289, a2 = 4.4407 Bohr) is on the
**molecular Gaussian** path via `cpu::scf::D3BJ` / `MeanFieldMethod::d3bj`
(PySCF `disp='d3bj'` twin; regression `d3bj-molecular-pbe`) **and on GPW**
via `GpwDft::d3bj` / `GPWEvaluator::d3bj` (CP2K CN + BJ pin;
`reg_gpw-pbe-d3bj-ice`). **DFT-D4**
(Caldeweyher et al.; two-body BJ + EEQ/CN-dependent C6) is on the molecular
path via `cpu::scf::D4` / `MeanFieldMethod::d4` (PySCF `disp='d4'` twin;
regression `dft-d4-h2o-pbe`; three-body ATM not included) **and on GPW** via
`GpwDft::d4` / `GPWEvaluator::d4` (dftd4 3.7.0 offline pin;
`reg_gpw-d4-water`; three-body ATM off). **Tkatchenko–Scheffler (TS)** free-atom
pairwise dispersion ships on GPW (`GpwDft::ts` / `GPWEvaluator::ts`; libMBD
tables; `reg_gpw-ts-ice`). Hirshfeld density volumes, MBD, and vdW-DF remain ○.

## 9. Composite DFT Methods

| Method    | CPU | CUDA | Kokkos |
| --------- | :-: | :--: | :---: |
| HF-3c     |  ○  |   ○  |   ○  |
| PBEh-3c   |  ○  |   ○  |   ○  |
| B97-3c    |  ○  |   ○  |   ○  |
| r2SCAN-3c |  ○  |   ○  |   ○  |

---

## 10. Periodic Kohn–Sham Methods

Two periodic DFT engines exist, with disjoint strengths:

- **Plane-wave DFT** (`cpu::pw::PW`, `cuda::pw`) — k-points, norm-conserving
  UPF pseudopotentials, block Davidson. **Ionic forces via analytic local HF +
  Ewald + Kleinman–Bylander nonlocal** (`pw_forces_analytic`, default for
  `PWEvaluator::with_gradient(true)`); FD (`pw_forces_fd`) is the oracle.
  **Cell stress** is **Nielsen–Martin analytic** (`pw_stress` → `kStress`;
  multi-SCF strain-FD oracle `pw_stress_fd`; `reg_pw-analytic-stress-si`).
  XC is **true GGA PBE** on the FFT grid (density gradients enter Vxc; default
  `xc = PBE`) or LSDA; **collinear spin** (`spin_polarized`) including
  **k-mesh UKS** (`reg_pw-uks-kmesh-fe`), **metallic smearing**
  (Fermi–Dirac / Methfessel–Paxton; Al pin + **bulk Fe UKS+smear**
  `reg_pw-uks-fe-smearing`) ship on CPU. **Bands / DOS** post-process
  (`cpu::pw::Bands` / `Dos`; `reg_bands-dos-si`). **PW CELL_OPT**
  (`IsotropicCellRelax` + intensive `kStress`; `reg_pw-variable-cell-si`).
  CUDA path remains LDA / spin-unpolarized.
- **GPW**, Gaussian and plane waves (`cpu::gpw::GpwDft` / `GPWEvaluator`) —
  Gaussian basis collocated on a uniform grid, GTH pseudopotentials, LDA
  (Padé) or **PBE / PBEsol / BLYP**, DFT-D3 zero + **BJ** + **D4** + **TS**,
  **analytic forces**, **analytic cell stress** (`gpw_stress`;
  `reg_gpw-analytic-stress` / `reg_gpw-stress-water`),
  product **OT** minimizer (`use_ot`), **Dudarev DFT+U** (`plus_u`;
  `reg_gpw-dftu-nio`). **Γ-point** RKS and **Γ-point UKS** (PADE LSD +
  **PBE-GGA** `reg_gpw-uks-pbe`; `n_alpha`/`n_beta`), plus **Monkhorst–Pack
  k-point RKS** (`k_mesh` / explicit `k_points`; phased lattice sums, KB
  nonlocal, per-k DIIS; `reg_gpw-kpoints-si`). **Open GPW follow-ons**
  (queued): k-point UKS, analytic k-point forces (WIP `wf/gpw-kpoint-forces`),
  k-point stress, Γ-UKS forces, revPBE/RPBE, multi-k MP2. CPU only.

| Method                           | CPU | CUDA | Kokkos |
| -------------------------------- | :-: | :--: | :---: |
| Γ-point periodic RKS             |  ✓  |   ○  |   ○  |
| k-point periodic RKS             |  ✓  |   ✓  |   ○  |
| Γ-point periodic UKS             |  ✓  |   ○  |   ○  |
| k-point periodic UKS             |  ✓  |   ○  |   ○  |
| Plane-wave DFT                   |  ✓  |   ✓  |   ○  |
| Gaussian and plane-wave DFT, GPW |  ✓  |   ○  |   ○  |
| GAPW DFT                         |  ○  |   ○  |   ○  |
| Periodic Gaussian-basis DFT      |  ✓  |   ○  |   ○  |
| Periodic LDA                     |  ✓  |   ✓  |   ○  |
| Periodic GGA                     |  ✓  |   ○  |   ○  |
| Periodic meta-GGA                |  ○  |   ○  |   ○  |
| Periodic hybrid DFT              |  ○  |   ○  |   ○  |
| Periodic range-separated hybrid  |  ○  |   ○  |   ○  |
| Periodic DFT+U                   |  ✓  |   ○  |   ○  |
| Periodic dispersion-corrected DFT |  ✓  |  ○  |   ○  |
| Periodic MP2                     |  ✓  |   ○  |   ○  |
| Periodic RPA                     |  ○  |   ○  |   ○  |
| Periodic GW                      |  ○  |   ○  |   ○  |

Which engine backs which ✓:

- **Γ-point periodic RKS** is GPW (`reg_gpw-water-pbe-mt`). **Γ-point periodic
  UKS** is GPW PADE LSD (`reg_gpw-uks-o2-gamma`). **k-point periodic RKS** is
  the plane-wave solver (Monkhorst–Pack mesh) on both CPU and CUDA, **and**
  GPW on CPU — the Gaussian-basis k-point path (phased lattice sums, KB
  nonlocal, per-k DIIS; `reg_gpw-kpoints-si`, diamond Si vs CP2K Quickstep,
  2026-07-25).
  **k-point periodic UKS** is CPU PW collinear spin on a multi-k mesh
  (`reg_pw-uks-kmesh-fe`; O₂ twin) plus **bulk Fe UKS + metallic smearing**
  (`reg_pw-uks-fe-smearing`). **Γ UKS PBE-GGA** on GPW: `reg_gpw-uks-pbe`.
- **Periodic LDA** is provided twice: Padé-LDA in GPW, and Slater + PW92 on the
  plane-wave grid (the latter is the CUDA ✓; CPU PW also has LDA). **Periodic
  GGA** is GPW PBE / **PBEsol** (`reg_gpw-pbesol-si`) / **BLYP**
  (`reg_gpw-blyp-water`) **and** CPU PW true PBE (`reg_pw-pbe-gga-si`).
- **Periodic DFT+U** is GPW Dudarev PLUS_U (`reg_gpw-dftu-nio`; NiO mol pin).
  Molecular Gaussian DFT+U remains ○.
- **Periodic dispersion-corrected DFT** is GPW + DFT-D3 zero damping,
  **D3(BJ)** (§8; `reg_gpw-pbe-d3bj-ice`), **D4** (`reg_gpw-d4-water`), and
  **TS** (`reg_gpw-ts-ice`).
- **Periodic MP2** is **Γ-point RI-GPW product** (`GpwDft::mp2` with
  `mp2_ri` default true; Del Ben et al. 2012 layout;
  `reg_periodic-ri-mp2-si` vs CP2K 2026.1; Poisson foothold
  `reg_periodic-mp2-si`). Multi-k MP2 residual throws; **periodic RPA** and
  **periodic GW** remain ○.
- Periodic hybrids are ○ for the same reason periodic HF is (see
  [Hartree–Fock](hartree-fock.md#2-periodic-hartreefock)): no periodic exact
  exchange.

## 11. Periodic Response and Analysis

| Method                                 | CPU | CUDA | Kokkos |
| -------------------------------------- | :-: | :--: | :---: |
| Periodic analytic forces               |  ✓  |   ○  |   ○  |
| Periodic analytic stress               |  ✓  |   ○  |   ○  |
| Density-functional perturbation theory |  ✓  |   ○  |   ○  |
| Finite-displacement phonons            |  ✓  |   ○  |   ○  |
| Berry-phase polarization               |  ✓  |   ○  |   ○  |
| Maximally localized Wannier functions  |  ✓  |   ○  |   ○  |

- **Periodic analytic forces** are GPW (Hellmann–Feynman + Pulay) and
  plane-wave **local HF + Ewald + KB nonlocal** (`pw_forces_analytic`). PW FD
  remains as oracle.
- **PW stress** is **Nielsen–Martin analytic** (`pw_stress` → Frame `kStress`;
  multi-SCF strain-FD oracle `pw_stress_fd`; `reg_pw-analytic-stress-si`).
  **GPW stress** is **analytic fixed-\(P\)+Pulay** product
  (`gpw_stress` / `GPWEvaluator::compute_stress`, `stress_use_fd=false`;
  multi-SCF oracle `gpw_stress_fd`; `reg_gpw-analytic-stress` dual-gate +
  `reg_gpw-stress-water` vs CP2K NUMERICAL). **GPW variable-cell relax ships**
  (evaluator re-target, stress axes mask, `StressConverged`;
  `reg_gpw-variable-cell-ice` vs CP2K CELL_OPT — see
  [Optimization](optimization.md)). **PW variable-cell relax** ships
  (`reg_pw-variable-cell-si`).
- **Finite-displacement phonons** at Γ: `cpu::vib::Phonon` expands a supercell,
  runs 6N force evaluations, and diagonalizes the mass-weighted Hessian
  (`reg_fd-phonon-supercell`; SW 64-atom Si twin). Molecular
  `cpu::vib::Hessian` remains the non-periodic path. **Phonon dispersions
  ω(q)** productize supercell IFC Fourier on the same SW twin
  (`cpu::vib::PhononDispersion`; `reg_phonon-dispersion-si` Γ/X/L windows).
- **DFPT** ships the **Γ-point (q = 0) foothold** for the PW method:
  `cpu::pw::DfptPhonon` — self-consistent Sternheimer response (conduction
  CG, Hartree + LDA kernel), analytic Ewald/KB second derivatives
  (`regressions/response/dfpt-phonon-si`; energy-second-difference twin to
  0.07 %). Restricted LSDA / NC-UPF / fixed occupations only; **finite-q
  DFPT Sternheimer** remains a follow-on (IFC-Fourier dispersions cover the
  classical q-path surface).
- **Berry-phase polarization**: King-Smith–Vanderbilt (KSV) electronic
  polarization on the PW method (`cpu::pw::BerryPolarization`; string-averaged
  Berry phases over occupied Bloch manifolds; symmetry-exact Si pins,
  `regressions/pw/berry-polarization`).
- **Maximally localized Wannier functions**: molecular **Foster–Boys** on
  closed-shell Gaussian RHF/RKS (`cpu::scf::WannierCenters`;
  `regressions/scf/wannier-centers`) **and** periodic multi-k
  **Marzari–Vanderbilt** on PW (`cpu::pw::PeriodicWannier`; bond-midpoint Si
  pin; `reg_periodic-wannier-si`). GPW-grid MLWF residual stays unqueued.

## 12. Periodic infrastructure

| Capability | CPU | CUDA | Kokkos |
| ---------- | :-: | :--: | :---: |
| Monkhorst–Pack k-mesh | ✓ | ✓ |   ○  |
| Block Davidson band solver | ✓ | ✓ |   ○  |
| Norm-conserving UPF pseudopotentials (KB projectors) | ✓ | ✓ |   ○  |
| GTH pseudopotentials (local + nonlocal) | ✓ | ○ |   ○  |
| FFT Poisson solver | ✓ | ✓ |   ○  |
| Martyna–Tuckerman non-periodic Poisson | ✓ | ○ |   ○  |
| Ewald ion–ion energy | ✓ | ✓ |   ○  |
| Pulay/DIIS density mixing | ✓ | ○ |   ○  |
| Orbital-transformation (OT) minimizer | ✓ | ○ |   ○  |
| Real-space multigrid collocation | ✓ | ○ |   ○  |
| Ultrasoft / PAW | ○ | ○ |   ○  |
| Smearing for metals | ✓ | ○ |   ○  |
| Band structure / DOS post-processing | ✓ | ○ |   ○  |

**Smearing for metals** on CPU PW: `PW::smearing` (Fermi–Dirac / Methfessel–Paxton)
+ `degauss` [Ry] + `extra_bands` (`reg_pw-smearing-al`). OT on GPW is a product
member (`use_ot`; `reg_ot-minimizer-water`). **Band structure / DOS** on CPU PW:
`cpu::pw::Bands` (high-symmetry path) + `cpu::pw::Dos` (MP mesh + Gaussian
smearing; `reg_bands-dos-si`).

On CUDA the plane-wave hot path (Davidson H·ψ, cuFFT, Poisson, V_eff, LDA,
density accumulation) is device-resident; Rayleigh–Ritz subspace bookkeeping and
the SCF control loop stay on the host.
