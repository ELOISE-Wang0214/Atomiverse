# Density-Functional Theory

Legend and caveats: [Roadmap index](index.md).

The molecular Kohn–Sham path shares the Gaussian SCF solver with
[Hartree–Fock](hartree-fock.md): the same 2×2 config product of
{restricted, unrestricted} × {HF, KS}, the same DIIS/level-shift machinery, the
same analytic-gradient assembly.

The molecular functional set is exactly four selectors — `LSDA`, `PBE`, `B3LYP`,
`PBE0` — evaluated on Becke-partitioned atomic grids with Lebedev angular
quadrature. All functionals are hand-coded and gated against libxc to 1e-12 in
tests; libxc is never a runtime dependency.

**Periodic DFT is DFT**, so it lives here too: sections 10–11 cover the
plane-wave and GPW engines, their functionals, forces, and the periodic-only
gaps (stress, phonons, hybrids).

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
| PBEsol     |  ○  |   ○  |   ○  |
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
`reg_blyp-water-sto3g` pins H2O/STO-3G vs PySCF xc=blyp.

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
| Tkatchenko–Scheffler dispersion  |  ○  |   ○  |   ○  |
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
regression `dft-d4-h2o-pbe`; three-body ATM not included). Periodic D4 remains ○.

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
  **Cell stress via strain-FD** (`pw_stress` → `kStress`). XC is **true GGA
  PBE** on the FFT grid (density gradients enter Vxc; default `xc = PBE`) or
  LSDA; **collinear spin** (`spin_polarized`) and **metallic smearing**
  (Fermi–Dirac / Methfessel–Paxton) ship on CPU. CUDA path remains LDA /
  spin-unpolarized.
- **GPW**, Gaussian and plane waves (`cpu::gpw::GpwDft` / `GPWEvaluator`) —
  Gaussian basis collocated on a uniform grid, GTH pseudopotentials, LDA
  (Padé) or PBE, DFT-D3 zero + **BJ**, **analytic forces**, product **OT**
  minimizer (`use_ot`). **Γ-point** RKS and **Γ-point UKS** (PADE LSD;
  `n_alpha`/`n_beta`). No k-points. CPU only.

| Method                           | CPU | CUDA | Kokkos |
| -------------------------------- | :-: | :--: | :---: |
| Γ-point periodic RKS             |  ✓  |   ○  |   ○  |
| k-point periodic RKS             |  ✓  |   ✓  |   ○  |
| Γ-point periodic UKS             |  ✓  |   ○  |   ○  |
| k-point periodic UKS             |  ○  |   ○  |   ○  |
| Plane-wave DFT                   |  ✓  |   ✓  |   ○  |
| Gaussian and plane-wave DFT, GPW |  ✓  |   ○  |   ○  |
| GAPW DFT                         |  ○  |   ○  |   ○  |
| Periodic Gaussian-basis DFT      |  ✓  |   ○  |   ○  |
| Periodic LDA                     |  ✓  |   ✓  |   ○  |
| Periodic GGA                     |  ✓  |   ○  |   ○  |
| Periodic meta-GGA                |  ○  |   ○  |   ○  |
| Periodic hybrid DFT              |  ○  |   ○  |   ○  |
| Periodic range-separated hybrid  |  ○  |   ○  |   ○  |
| Periodic DFT+U                   |  ○  |   ○  |   ○  |
| Periodic dispersion-corrected DFT |  ✓  |  ○  |   ○  |
| Periodic MP2                     |  ○  |   ○  |   ○  |
| Periodic RPA                     |  ○  |   ○  |   ○  |
| Periodic GW                      |  ○  |   ○  |   ○  |

Which engine backs which ✓:

- **Γ-point periodic RKS** is GPW (`reg_gpw-water-pbe-mt`). **Γ-point periodic
  UKS** is GPW PADE LSD (`reg_gpw-uks-o2-gamma`). **k-point periodic RKS** is
  the plane-wave solver (Monkhorst–Pack mesh) — on both CPU and CUDA. There is
  no periodic *Gaussian-basis* k-point path. **k-point periodic UKS** remains ○.
- **Periodic LDA** is provided twice: Padé-LDA in GPW, and Slater + PW92 on the
  plane-wave grid (the latter is the CUDA ✓; CPU PW also has LDA). **Periodic
  GGA** is GPW/PBE **and** CPU PW true PBE (`reg_pw-pbe-gga-si`).
- **Periodic dispersion-corrected DFT** is GPW + DFT-D3 zero damping and
  **D3(BJ)** (§8; `reg_gpw-pbe-d3bj-ice`).
- Periodic hybrids are ○ for the same reason periodic HF is (see
  [Hartree–Fock](hartree-fock.md#2-periodic-hartreefock)): no periodic exact
  exchange.

## 11. Periodic Response and Analysis

| Method                                 | CPU | CUDA | Kokkos |
| -------------------------------------- | :-: | :--: | :---: |
| Periodic analytic forces               |  ✓  |   ○  |   ○  |
| Periodic analytic stress               |  ✓  |   ○  |   ○  |
| Density-functional perturbation theory |  ○  |   ○  |   ○  |
| Finite-displacement phonons            |  ○  |   ○  |   ○  |
| Berry-phase polarization               |  ○  |   ○  |   ○  |
| Maximally localized Wannier functions  |  ○  |   ○  |   ○  |

- **Periodic analytic forces** are GPW (Hellmann–Feynman + Pulay) and
  plane-wave **local HF + Ewald + KB nonlocal** (`pw_forces_analytic`). PW FD
  remains as oracle. **PW stress** ships as strain-FD of energy (`pw_stress`);
  GPW stress remains open.
- **PW stress** is available via strain-FD (`pw_stress` → Frame `kStress`).
  GPW still has no stress tensor; full variable-cell drivers remain productization
  work (see [Optimization](optimization.md)).
- **Finite-displacement phonons** is ○ even though a finite-difference Hessian
  exists: `cpu::vib::Hessian` is molecular (3N Cartesian displacements, no
  supercell or q-point machinery), so it yields no phonon band structure.

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
| Band structure / DOS post-processing | ○ | ○ |   ○  |

**Smearing for metals** on CPU PW: `PW::smearing` (Fermi–Dirac / Methfessel–Paxton)
+ `degauss` [Ry] + `extra_bands` (`reg_pw-smearing-al`). OT on GPW is a product
member (`use_ot`; `reg_ot-minimizer-water`). Band structure / DOS post-processing
remains ○ (`bands-dos-si`).

On CUDA the plane-wave hot path (Davidson H·ψ, cuFFT, Poisson, V_eff, LDA,
density accumulation) is device-resident; Rayleigh–Ritz subspace bookkeeping and
the SCF control loop stay on the host.
