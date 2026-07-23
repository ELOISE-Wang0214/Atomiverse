# PySCF method queue

Durable worklist for workflow `pyscf-method-loop`. Do not delete done rows.
The census phase rewrites **pending** from `docs/roadmap/` ○ rows that PySCF
can oracle offline; **done** / **blocked** history is preserved across runs.

Census date: 2026-07-23. Roadmap CPU ○ ∩ PySCF oracle; excludes classical MD /
FF / sampling / CUDA-Kokkos / ML / PW-GPW-without-Gaussian-twin / parked
hf-parity product-depth / SE tables without PySCF twin / meta-GGA (no τ) /
range-separated (no erf ERI). Axis extensions first after closed UMP2 + hybrid
TDA: frozen-core MP2, CIS(D), SCS/SOS-MP2, MP2 grad, full LR-TDDFT, BLYP,
GHF grad, ROMP2, ROKS; CPHF polarizability / open-shell Hessians late.

## pending

- [ ] `sos-mp2-h2o-sto3g` | SOS-MP2 | post-hf | water | sto-3g | e_corr | opposite-spin only scaling of RMP2
- [ ] `mp2-grad-h2o-sto3g` | MP2 analytic gradient | post-hf | water | sto-3g | forces Ha/Bohr | RMP2 nuclear gradient; PySCF mp.MP2.nuc_grad_method
- [ ] `romp2-o2-sto3g` | ROMP2 | post-hf | o2 | sto-3g | e_corr | ROHF ref ✓ + analytic ROHF grad ✓; MP2 on ROHF ○
- [ ] `adc2-h2co-sto3g` | ADC(2) | excited-states | h2co | sto-3g | S1 eV | algebraic-diagrammatic foothold; PySCF adc.ADC
- [ ] `spin-flip-tda-o2` | Spin-flip TDA | excited-states | o2 | sto-3g | excitation eV | SF-TDA on UHF/ROHF; PySCF TDA spin-flip
- [ ] `uks-tda-o2-lsda` | UKS TDA | excited-states | o2 | sto-3g | S/T roots | open-shell TDDFT-TDA; RKS::tda ✓ only
- [ ] `uhf-hessian-o2` | UHF analytic Hessian | hartree-fock | o2 | sto-3g | Hessian eV/Å² | RHF CPHF Hessian ✓ (reg_h2o_rhf_freq); UHF Hessian ○
- [ ] `rohf-hessian-o2` | ROHF analytic Hessian | hartree-fock | o2 | sto-3g | Hessian eV/Å² | ROHF energy+grad ✓; analytic Hessian ○; PySCF hessian.ROHF
- [ ] `mp3-h2-sto3g` | MP3 | post-hf | h2 | sto-3g | e_corr | third-order PT; PySCF mp.MP3
- [ ] `ccsd-h2-sto3g` | CCSD | post-hf | h2 | sto-3g | e_corr+e_tot | coupled-cluster doubles foothold; PySCF cc.CCSD
- [ ] `ccsd-t-h2o-sto3g` | CCSD(T) | post-hf | water | sto-3g | e_corr | needs CCSD; PySCF cc.CCSD.ccsd_t
- [ ] `casci-h2-sto3g` | CASCI | post-hf | h2 | sto-3g | e_tot | active-space on existing FCI kernel; PySCF mcscf.CASCI
- [ ] `g0w0-h2-sto3g` | G0W0 | excited-states | h2 | sto-3g | QP gaps | Green's-function foothold; PySCF gw.G0W0
- [ ] `d3bj-molecular-pbe` | DFT-D3 BJ damping | dft | water | sto-3g | E_disp | molecular Gaussian D3 ○ (GPW zero-damping only); PySCF dftd3
- [ ] `b2plyp-h2o` | B2PLYP double hybrid | dft | water | sto-3g | e_tot | hybrid KS + MP2 composition; ingredients exist separately
- [ ] `delta-scf-n2` | Delta-SCF | excited-states | n2 | sto-3g | ΔE | MOM/IMOM foothold; PySCF mom_occ_; not hf-parity parking product depth if single textbook case

## in_progress

- [~] `scs-mp2-h2o-sto3g` | SCS-MP2 | post-hf | water | sto-3g | e_corr | Grimme OS/SS scale of existing RMP2 spin components
- [~] `full-lr-tddft-h2-pbe` | Full linear-response TDDFT | excited-states | h2 | sto-3g | S1 eV | RKS beyond TDA (Casida A/B); RHF::tdhf ✓, pure/hybrid KS full TDDFT ○
- [~] `roks-o2-sto3g` | ROKS | dft | o2 | sto-3g | e_tot | restricted open-shell KS; PySCF dft.ROKS
- [~] `ief-pcm-h2o` | IEF-PCM | multiscale | water | 6-31g* | e_solv | conductor-like only today; PySCF solvent IEF-PCM


## blocked

- [!] `eom-ee-ccsd` | blocked upstream — no CCSD product path yet (queue `ccsd-h2-sto3g` first)
- [!] `cam-b3lyp` | blocked — range-separated / erf-attenuated ERI engine absent (Rys has no attenuation)
- [!] `periodic-rhf-gamma` | blocked — no periodic exact exchange; PW/GPW are semi-local DFT only
- [!] `td-dftb` | blocked upstream — DFTB2 has no response path; not a PySCF twin
- [!] `gfn2-xtb` | blocked — no PySCF oracle; xTB is external (GFN0 toy already ✓)
- [!] `pm6-pm7` | blocked — NDDO tables / MOPAC twin, not PySCF

## done
- [x] `static-polariz-h2o` | RHF static polarizability (CPHF) | regressions/scf/static-polariz-h2o | H2O/STO-3G α_ii=(0.04190863, 5.20973255, 2.12715820) a.u. (PySCF 2.14.0 CPHF/FD); `RHF::polarizability`; tol 1e-5 a.u.
- [x] `pcm-grad-h2o` | Implicit-solvent analytic gradient | regressions/scf/pcm-grad-h2o | FD of free energy ≤1e-5 Ha/Bohr; PySCF pin in REFERENCE
- [x] `ghf-grad-h3-sto3g` | GHF analytic gradient | regressions/scf/ghf-grad-h3-sto3g | noncollinear H3 forces 1e-5 Ha/Bohr
- [x] `blyp-water-sto3g` | BLYP GGA | regressions/scf/blyp-water-sto3g | H2O/STO-3G cart e_tot=-75.277447482665 Ha (PySCF 2.14.0 xc=blyp); tol 1e-6 Ha
- [x] `cisd-h2co-sto3g` | CIS(D) | regressions/scf/cisd-h2co-sto3g | H2CO/STO-3G S1 CIS(D)=4.54457685 eV (PySCF+Head-Gordon U+V); tol 1e-4 eV
- [x] `frozen-core-mp2-h2o` | Frozen-core RMP2 | regressions/scf/frozen-core-mp2-h2o | H2O/6-31G* cart frozen=1 (O 1s); PySCF 2.14.0 e_corr=-0.186119540118 Ha, e_tot=-76.196624528375 Ha; `rhf.frozen=1`; tol 1e-6 Ha

- [x] `rhf-analytic-hessian-h2o` | RHF analytic Hessian (CPHF) | regressions/scf/h2o_rhf_freq | PySCF analytic Hessian freq ≤1 cm⁻¹
- [x] `lowdin-mk-chelpg-h2o` | Löwdin + MK + CHELPG charges | regressions/scf/h2o_charges_pyscf | PySCF 2.14.0 pins in REFERENCE.md
- [x] `rhf-dipole-h2o` | RHF dipole | regressions/scf/h2o_rhf_dipole | PySCF dipole pin
- [x] `rmp2-h2o-sto3g` | RMP2 | regressions/scf/h2o_mp2_sto3g | PySCF mp.MP2
- [x] `ri-mp2-h2o` | RI-MP2 | regressions/scf/h2o_ri_mp2 | density-fitted MP2
- [x] `ump2-o2-sto3g` | UMP2 | regressions/scf/ump2-o2-sto3g | UMP2/O2/STO-3G triplet cart conv_tol=1e-12 PySCF 2.14.0: e_corr=-0.096393144130 Ha, e_tot=-147.730339929580 Ha (UHF e_hf=-147.633946785450 Ha); tol 1e-6 Ha
- [x] `rhf-tda-tdhf-ch2o` | RHF TDA/TDHF | regressions/scf/ch2o_tdhf | CIS/TDHF vs PySCF
- [x] `rks-tda-alda` | RKS TDA full-grid ALDA | regressions/scf/h2_full_grid_alda | LSDA ALDA kernel
- [x] `hybrid-tda-b3lyp-h2` | Hybrid TDDFT-TDA | regressions/scf/h2-b3lyp-hybrid-tda | B3LYP exact-exchange kernel; S1=25.928364822972 eV (PySCF 2.14.0)
- [x] `analytic-tda-force-h2` | CIS/TDA analytic forces | regressions/scf/h2_analytic_tda_force | analytic default
- [x] `tda-nac-h2` | TDA NAC vectors | regressions/scf/h2_nac_vectors | true NAC
- [x] `rohf-o2-energy-grad` | ROHF energy + gradient | regressions/scf/o2_rohf_sto3g, o2_rohf_grad | PySCF ROHF
- [x] `ghf-h3-energy` | GHF energy | regressions/scf/h3_ghf_sto3g | energy-only
- [x] `cpcm-h2o` | C-PCM/COSMO energy | regressions/scf/h2o_rhf_631gs_cpcm | conductor-like solvent
