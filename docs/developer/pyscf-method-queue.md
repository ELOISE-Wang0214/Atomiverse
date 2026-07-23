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

(empty)

## in_progress

- [~] `ccsd-t-h2o-sto3g` | CCSD(T) | post-hf | water | sto-3g | e_corr | CCSD energy foothold lands (`RHF::ccsd`); triples (T) next; PySCF cc.CCSD.ccsd_t
- [~] `casci-h2-sto3g` | CASCI | post-hf | h2 | sto-3g | e_tot | active-space on existing FCI kernel; PySCF mcscf.CASCI

- [~] `eom-ee-ccsd` | EOM-EE-CCSD | excited-states | h2 | sto-3g | ω | CCSD energy foothold now ships; EOM layer still ○; PySCF cc.eom_rccsd

## blocked

- [!] `cam-b3lyp` | blocked — range-separated / erf-attenuated ERI engine absent (Rys has no attenuation)
- [!] `periodic-rhf-gamma` | blocked — no periodic exact exchange; PW/GPW are semi-local DFT only
- [!] `td-dftb` | blocked upstream — DFTB2 has no response path; not a PySCF twin
- [!] `gfn2-xtb` | blocked — no PySCF oracle; xTB is external (GFN0 toy already ✓)
- [!] `pm6-pm7` | blocked — NDDO tables / MOPAC twin, not PySCF

## done

- [x] `delta-scf-n2` | Delta-SCF | regressions/scf/delta-scf-n2 | N2/STO-3G α HOMO→LUMO UHF-IMOM; PySCF 2.14.0 mom_occ_ ΔE=8.24694070 eV; `UHF::delta_scf` + `MomStyle::IMOM`; tol 1e-4 eV; ratio 7.83
- [x] `adc2-h2co-sto3g` | ADC(2) | regressions/scf/adc2-h2co-sto3g | H2CO/STO-3G cart strict ADC(2) EE; PySCF 2.14.0 S1=4.535317180049 eV; `RHF::adc2`; tol 1e-4 eV; ratio 0.267
- [x] `mp3-h2-sto3g` | MP3 | regressions/scf/mp3-h2-sto3g | H2/STO-3G cart RMP3; offline pin e_corr=-0.017974146227 Ha (Szabo–Ostlund on PySCF RHF orbs); `RHF::mp3`; tol 1e-6 Ha; ratio 0.774
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
