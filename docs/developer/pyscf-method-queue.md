# PySCF method queue

Durable worklist for workflow `pyscf-method-loop`. Do not delete done rows.
The census phase rewrites **pending** from `docs/roadmap/` ○ rows that PySCF
can oracle offline; **done** / **blocked** history is preserved across runs.

Census date: 2026-07-23 (re-census after wave-6 exhaustion). Roadmap CPU ○ ∩
PySCF oracle; excludes classical MD / FF / sampling / CUDA-Kokkos / ML /
PW-GPW-without-Gaussian-twin / parked hf-parity product-depth / SE tables
without PySCF twin / meta-GGA (no τ) / range-separated (no erf ERI) /
CASPT2 (no stock PySCF). Prefer new family or scale over XC variants
(BP86/B3PW91/B2GP-PLYP/…). Next axis after closed EOM-EE / CASCI / CCSD(T) /
G0W0 / B2PLYP: EOM-IP/EA, UKS full-LR, ADC(2)-x/IP-ADC, CASSCF/NEVPT2, CC
gradients, KS Hessians, open-shell MP2 forces, SF-TDDFT, MP4, CC2/EOM-CC2,
dRPA/BSE, D4, EOM-SF, GF2.

## pending

- [ ] `eom-ip-ccsd` | EOM-IP-CCSD | excited-states | h2o | sto-3g | IP eV | EOM-EE ✓; ionization channel; PySCF cc.eom_rccsd IP
- [ ] `eom-ea-ccsd` | EOM-EA-CCSD | excited-states | h2 | sto-3g | EA eV | EOM-EE ✓; attachment channel; PySCF cc.eom_rccsd EA
- [ ] `uks-full-lr-o2-lsda` | UKS full-LR TDDFT | excited-states | o2 | sto-3g | roots eV | UKS-TDA ✓; RKS full-LR ✓; open-shell Casida A/B; PySCF tdscf.TDDFT(UKS)
- [ ] `adc2x-h2co-sto3g` | ADC(2)-x | excited-states | h2co | sto-3g | S1 eV | strict ADC(2) ✓; extended ADC(2)-x; PySCF adc method=adc(2)-x
- [ ] `ip-adc2-h2o-sto3g` | IP-ADC(2) | excited-states | h2o | sto-3g | IP eV | ADC(2) EE ✓; ionization ADC; PySCF adc IP
- [ ] `ea-adc2-h2-sto3g` | EA-ADC(2) | excited-states | h2 | sto-3g | EA eV | ADC(2) EE ✓; attachment ADC; PySCF adc EA
- [ ] `casscf-h2-sto3g` | CASSCF | post-hf | h2 | sto-3g | e_tot | CASCI ✓; orbital optimization; PySCF mcscf.CASSCF
- [ ] `nevpt2-h2-sto3g` | NEVPT2 | post-hf | h2 | sto-3g | e_corr | needs CASSCF; partially/strongly contracted; PySCF mrpt.NEVPT
- [ ] `ccsd-grad-h2o-sto3g` | CCSD analytic gradient | post-hf | water | sto-3g | forces Ha/Bohr | CCSD energy ✓; PySCF cc.CCSD.nuc_grad_method
- [ ] `rks-hessian-h2o-sto3g` | RKS analytic Hessian | dft | water | sto-3g | Hessian eV/Å² | RHF/UHF CPHF Hessian ✓; KS Hessian ○; PySCF hessian.RKS
- [ ] `uks-hessian-o2-sto3g` | UKS analytic Hessian | dft | o2 | sto-3g | Hessian eV/Å² | UKS energy+TDA ✓; analytic Hessian ○; PySCF hessian.UKS
- [ ] `ump2-grad-o2-sto3g` | UMP2 analytic gradient | post-hf | o2 | sto-3g | forces Ha/Bohr | RMP2 grad ✓; UMP2 energy ✓; open-shell MP2 forces; PySCF UMP2.nuc_grad_method
- [ ] `sf-tddft-o2` | Spin-flip TDDFT | excited-states | o2 | sto-3g | excitation eV | SF-TDA ✓; full SF-RPA; PySCF TDDFT spin-flip
- [ ] `mp4-h2-sto3g` | MP4(SDTQ) | post-hf | h2 | sto-3g | e_corr | MP3 ✓; fourth-order PT; PySCF mp.MP4
- [ ] `cc2-h2-sto3g` | CC2 | post-hf | h2 | sto-3g | e_corr | approximate CC doubles foothold; unlocks EOM-CC2; PySCF cc.CC2
- [ ] `eom-cc2-h2` | EOM-CC2 | excited-states | h2 | sto-3g | ω eV | needs CC2 energy path; PySCF EOM-CC2
- [ ] `direct-rpa-h2-sto3g` | Direct RPA | excited-states | h2 | sto-3g | e_corr | G0W0 screening twin as standalone corr; PySCF rpa.RPA
- [ ] `bse-h2-sto3g` | Bethe–Salpeter BSE | excited-states | h2 | sto-3g | S1 eV | after G0W0 ✓; PySCF tdscf.BSE
- [ ] `dft-d4-h2o-pbe` | DFT-D4 | dft | water | sto-3g | E_disp | D3 zero + BJ ✓; D4 open; PySCF/dftd4 offline pin
- [ ] `eom-sf-ccsd-o2` | EOM-SF-CCSD | excited-states | o2 | sto-3g | ω eV | EOM-EE ✓; spin-flip EOM channel; PySCF eom_uccsd SF
- [ ] `gf2-h2-sto3g` | GF2 | excited-states | h2 | sto-3g | QP | second-order Green's function; PySCF gf2

## in_progress

(empty)

## blocked

- [!] `cam-b3lyp` | blocked — range-separated / erf-attenuated ERI engine absent (Rys has no attenuation)
- [!] `periodic-rhf-gamma` | blocked — no periodic exact exchange; PW/GPW are semi-local DFT only
- [!] `td-dftb` | blocked upstream — DFTB2 has no response path; not a PySCF twin
- [!] `gfn2-xtb` | blocked — no PySCF oracle; xTB is external (GFN0 toy already ✓)
- [!] `pm6-pm7` | blocked — NDDO tables / MOPAC twin, not PySCF
- [!] `caspt2` | blocked — no stock PySCF oracle (OpenMolcas / external); NEVPT2 is the in-tree PySCF MRPT path

## done

- [x] `eom-ee-ccsd` | EOM-EE-CCSD | regressions/scf/eom-ee-ccsd | H2/STO-3G cart singlet EOM-EE; PySCF 2.14.0 S1=26.365966612792 eV; `RHF::eom_ee_ccsd`; tol 1e-4 eV; ratio 0.193
- [x] `delta-scf-n2` | Delta-SCF | regressions/scf/delta-scf-n2 | N2/STO-3G α HOMO→LUMO UHF-IMOM; PySCF 2.14.0 mom_occ_ ΔE=8.24694070 eV; `UHF::delta_scf` + `MomStyle::IMOM`; tol 1e-4 eV; ratio 7.83
- [x] `casci-h2-sto3g` | CASCI | regressions/scf/casci-h2-sto3g | H2 CASCI(2,2)/STO-3G cart; PySCF 2.14.0 e_tot=-1.137283834489 Ha; `RHF::casci(frame,2,2)`; tol 1e-6 Ha; ratio 0.076
- [x] `adc2-h2co-sto3g` | ADC(2) | regressions/scf/adc2-h2co-sto3g | H2CO/STO-3G cart strict ADC(2) EE; PySCF 2.14.0 S1=4.535317180049 eV; `RHF::adc2`; tol 1e-4 eV; ratio 0.267
- [x] `mp3-h2-sto3g` | MP3 | regressions/scf/mp3-h2-sto3g | H2/STO-3G cart RMP3; offline pin e_corr=-0.017974146227 Ha (Szabo–Ostlund on PySCF RHF orbs); `RHF::mp3`; tol 1e-6 Ha; ratio 0.774
- [x] `ccsd-t-h2o-sto3g` | CCSD(T) | regressions/scf/ccsd-t-h2o-sto3g | H2O/STO-3G cart CCSD+(T); PySCF 2.14.0 e_corr=-0.049505972715 Ha, e_t=-6.740968e-05 Ha; `RHF::ccsd_t`; tol 1e-6 Ha; ratio 0.333
- [x] `ccsd-h2-sto3g` | CCSD | regressions/scf/ccsd-h2-sto3g | H2/STO-3G cart CCSD; PySCF 2.14.0 e_corr=-0.020524527092 Ha, e_tot=-1.137283834489 Ha; `RHF::ccsd`; tol 1e-6 Ha; ratio 0.317
- [x] `g0w0-h2-sto3g` | G0W0 | regressions/scf/g0w0-h2-sto3g | H2 RHF/STO-3G spectral G0W0 QP gap=35.007583892680 eV (PySCF 2.14.0 GWExact+dRPA); `RHF::g0w0`; tol 1e-4 eV; ratio 0.109
- [x] `b2plyp-h2o` | B2PLYP double hybrid | regressions/scf/b2plyp-h2o | H2O/STO-3G cart; PySCF 2.14.0 composed E_KS+a_c·MP2 e_tot=-75.210908483134 Ha; `XCType::B2PLYP` a_x=0.53 a_c=0.27 + `RKS::mp2`; tol 1e-6 Ha; ratio 1.75
- [x] `romp2-o2-sto3g` | ROMP2 | regressions/scf/romp2-o2-sto3g | O2 triplet ROHF→ROMP2/STO-3G cart; PySCF 2.14.0 UMP2(ROHF) e_corr=-0.096718716392 Ha, e_tot=-147.728885707074 Ha; `ROHF::mp2`; tol 1e-6 Ha; ratio 1.37
- [x] `roks-grad-o2-sto3g` | ROKS analytic gradient | regressions/scf/roks-grad-o2-sto3g | O2 ROKS-LSDA/STO-3G cart forces; PySCF 2.14.0 g_z=±0.162926534046 Ha/Bohr; `ROKS::forces` spin-Fock W; tol 2e-5 Ha/Bohr; ratio 1.117
- [x] `h2-b3lyp-rpa-tddft` | Hybrid full-LR TDDFT (B3LYP RPA) | regressions/scf/h2-b3lyp-rpa-tddft | H2 RKS-B3LYP/STO-3G Casida A/B; PySCF 2.14.0 S1=25.307797345333 eV; `RKS::tddft`; tol 1e-4 eV; ratio 0.430
- [x] `spin-flip-tda-o2` | Spin-flip TDA | regressions/scf/spin-flip-tda-o2 | O2 UHF/STO-3G collinear SF-CIS; offline pin roots 0.0451/1.5677/1.5677/3.0384/7.8486 eV; `UHF::sf_tda`; tol 1e-4 eV; ratio 9.56
- [x] `rohf-hessian-o2` | ROHF nuclear Hessian | regressions/scf/rohf-hessian-o2 | O2 ROHF/STO-3G FD of analytic forces; stretch H_zz~135.57 eV/Å²; `ROHF::hessian`; tol 1e-5 eV/Å²; perf_exempt (no PySCF hessian.ROHF)
- [x] `mp2-grad-h2o-sto3g` | MP2 analytic gradient | regressions/scf/mp2-grad-h2o-sto3g | H2O/STO-3G cart RMP2; PySCF 2.14.0 g_O,z=-9.528798e-02 Ha/Bohr; `RHF::mp2_forces`; tol 1e-5 Ha/Bohr
- [x] `sos-mp2-h2o-sto3g` | SOS-MP2 | regressions/scf/sos-mp2-h2o-sto3g | H2O/STO-3G cart Jung c_OS=1.3 c_SS=0; PySCF 2.14.0 e_corr=-0.043570566111 Ha, e_tot=-75.006593704574 Ha; `mp2_os_scale`/`mp2_ss_scale` members; tol 1e-6 Ha
- [x] `uks-tda-o2-lsda` | UKS TDA | regressions/scf/uks-tda-o2-lsda | O2 UKS-LSDA/STO-3G cart spin-conserving TDA; PySCF 2.14.0 roots 7.6209/7.6424/7.6425/8.4033/8.4033 eV; `UKS::tda`; tol 1e-4 eV
- [x] `d3bj-molecular-pbe` | DFT-D3 BJ damping | regressions/scf/d3bj-molecular-pbe | H2O molecular PBE-D3(BJ) E_disp=-0.000359468656 Ha (PySCF 2.14.0+pyscf-dispersion); `cpu::scf::D3BJ` + `MeanFieldMethod::d3bj`; tol 1e-6 Ha
- [x] `uhf-hessian-o2` | UHF analytic Hessian | regressions/scf/uhf-hessian-o2 | O2 UHF/STO-3G cart CPHF Hessian; PySCF 2.14.0 `hessian.UHF`; `UHF::hessian`; tol 1e-5 eV/Å²
- [x] `scs-mp2-h2o-sto3g` | SCS-MP2 | regressions/scf/scs-mp2-h2o-sto3g | H2O/STO-3G cart Grimme c_OS=1.2 c_SS=1/3; PySCF 2.14.0 e_corr=-0.040895594623 Ha, e_tot=-75.003918733086 Ha; `mp2_os_scale`/`mp2_ss_scale` members; tol 1e-6 Ha
- [x] `full-lr-tddft-h2-pbe` | Full linear-response TDDFT | regressions/scf/full-lr-tddft-h2-pbe | RKS::tddft Casida A/B PBE; S1=25.519450223124 eV (PySCF 2.14.0)
- [x] `roks-o2-sto3g` | ROKS | regressions/scf/roks-o2-sto3g | O2 triplet ROKS/LSDA/STO-3G cart e_tot=-147.193715448518 Ha (PySCF 2.14.0 dft.ROKS xc=LDA,VWN); tol 1e-6 Ha; spin purity ⟨S²⟩=2
- [x] `static-polariz-h2o` | RHF static polarizability (CPHF) | regressions/scf/static-polariz-h2o | H2O/STO-3G α_ii=(0.04190863, 5.20973255, 2.12715820) a.u. (PySCF 2.14.0 CPHF/FD); `RHF::polarizability`; tol 1e-5 a.u.
- [x] `ief-pcm-h2o` | IEF-PCM | regressions/scf/ief-pcm-h2o | H2O RHF/6-31G* cart ε=78.4 Bondi×1.2 Lebedev 302; PySCF 2.14.0 e_solv=-0.01265803490853 Ha; `solvent_model=IEFPCM`; tol 1e-4 Ha (hard bury-delete vs SWIG)
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
