# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

## Live — esmm campaign: ES × MM connection layer (2026-07-26, 16 specs)

Operator-approved three-wave campaign (plan-mode design review + Q&A 2026-07-26):
unify the force-pipeline composition contract, redesign QM/MM around an engine
slot + `qm_z` mask column, generalize CPMD into a shared `cpu/electronic/`
propagation layer. Execution order (operator, 2026-07-26): **esmm-cpmd wave
first** (`/mol:impl-all esmm-cpmd`), then contract + qmmm waves;
esmm-contract-05 rebases its gpw/scf `=`→`+=` sites over the cpmd wave.
The libcint/C++20 sweep gate was cleared by committing it (2026-07-26).

| # | Spec | Deps | One-line |
| - | ---- | ---- | -------- |
| W1-1 | esmm-contract-01-zeroing | — | Zero* unification (Hessian/Path/MC/Minimize/Evaluate) + k/2 oracle fix |
| W1-2 | esmm-contract-02-validation-walks | — | McBase + PathDriver node-aware energy-source validation |
| W1-3 | esmm-contract-03-runslot-nested | — | run_slot_ dispatches nested drivers (line-search zero-force fix) + R2 |
| W1-4 | esmm-contract-04-forces-requested | W1-2 | phase_wants_forces plumbing + ASPC-under-MC guard + R3 |
| W1-5 | esmm-contract-05-accumulate-switch | W1-1 | all ES sources `=`→`+=` + kEnergyElectronic + R1 |
| W2-1 | esmm-qmmm-01-embed-forces | — | analytic embedding gradient + MM-site forces (third law) |
| W2-2 | esmm-qmmm-02-qm-engine-protocol | W2-1 | QmEngine/EmbeddingEngine protocols; MFM×7 + AM1 + DFTB2 |
| W2-3 | esmm-qmmm-03-electrostatic-rewrite | W2-2, W1-4, W1-5 | qm_z mask column + engine slot + real MM coupling |
| W2-4 | esmm-qmmm-04-mech-pol-port | W2-3 | Mechanical/Polarizable onto QmSubsystem; quantitative re-pins |
| W2-5 | esmm-qmmm-05-md-closure | W2-4 | QM/MM-MD conservation regression; roadmap ○→✓; tutorials |
| W3-1 | esmm-cpmd-01-electronic-home | — | DONE 2026-07-26 (closed; see git c40acdb3) — cpu/electronic/ shared propagators + linalg orbital_tangent extraction |
| W3-2 | esmm-cpmd-02-timestep-protocol | W3-1 | DONE 2026-07-27 (closed; see git) — TimestepProvider/Consumer; cpmd_dt deleted; cpmd×stress gate |
| W3-3 | esmm-cpmd-03-observability | W3-2 | kEnergyElectronKinetic/Thermostat + CpmdConserved → Thermo drift |
| W3-4 | esmm-cpmd-04-ao-cpmd | W3-1..3 | Gaussian AO-CPMD on MeanFieldMethod (fixed-density eval + W=occ·C(CᵀFC)Cᵀ) |
| W3-5 | esmm-cpmd-05-scf-xlbomd | W3-4 | XL-BOMD on the SCF methods (guess ladder) |
| W3-6 | esmm-cpmd-06-single-eval | W3-3 | GpwDft::overlap(); CPMD single KS eval per step (~2×) |

## Live — deferred PySCF engine gaps (2026-07-24)

Approved for **later** implementation (do **not** auto `/mol:impl` until
scheduled). Unblock `docs/developer/pyscf-method-queue.md` blocked rows.

| Spec | Status | Unblocks | One-line |
| ---- | ------ | -------- | -------- |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/spec.md) | approved, deferred | `cam-b3lyp` | erf-attenuated ERI + CAM-B3LYP RKS energy foothold |
| [periodic-rhf-gamma](periodic-rhf-gamma/spec.md) | approved, deferred | `periodic-rhf-gamma` | Γ-point periodic RHF energy (lattice HF-X) foothold |

Acceptance: `rsh-eri-cam-b3lyp/acceptance.md`, `periodic-rhf-gamma/acceptance.md`.
Independent — no cross-deps. Roadmap: `docs/roadmap/dft.md` §6, `hartree-fock.md` §2.

## Closed — SCF review remediation chain (2026-07-24, all 8 done)

Fixed all 27 findings of the 2026-07-24 /mol:review of `cpu/scf` (verdict
was REQUEST CHANGES: 11 🔴 / 9 🟡 / 7 🟢). Both architecture scrubs restored
to green. Specs deleted on close; history in git (commits f9a6443a..HEAD).

| # | Spec | Deps | One-line |
| - | ---- | ---- | -------- |
| 01 | scf-review-01-linalg | — | DONE 2026-07-24 (closed; see git) — sym_power/eig_nonsymmetric/solve_pcg in cpu/linalg + OrbitalHessian::direct_rpa |
| 02 | scf-review-02-eri-hoists | 01 | DONE 2026-07-24 (closed; see git) — internal/ao2mo.h + internal/spin_orbital_eri.h + mp_denominator |
| 03 | scf-review-03-posthf-dedup | 01, 02 | DONE 2026-07-24 (closed; see git) — eom_ccsd_common (EomImds ctor) + ci_engine; build_imds scrub debt cleared |
| 04 | scf-review-04-meanfield-dedup | — | DONE 2026-07-24 (closed; see git) — FD-Hessian CRTP template + check_hessian_preconditions + GradientResult tag ctor |
| 05 | scf-review-05-api-hygiene | 01 | DONE 2026-07-24 (closed; see git) — grid ctor + throws + units::pi + .cpp moves; BOTH scrubs restored green |
| 06 | scf-review-06-scf-perf | 05 | DONE 2026-07-24 (closed; see git) — XCGridBasis cache + jk_rhf one-sweep + ERI micro-opts (~14% on fixture) |
| 07 | scf-review-07-posthf-perf | 02 | DONE 2026-07-24 (closed; see git) — block eri_dense (10x/54x) + RI-MP2 GEMM + copy elision |
| 08 | scf-review-08-ccsd-spatial-diis | 02, 03, 07 | DONE 2026-07-24 (closed; see git) — CcAmplitudeDiis + spatial CCSD (~87x, oracle parity 1e-9) |

## Live — deferred PySCF only (after harness wave 12)

| Spec | Branch / track | Status |
| ---- | -------------- | ------ |
| [periodic-rhf-gamma](periodic-rhf-gamma/) | pyscf / HF periodic | open (deferred) |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/) | pyscf / RSH | open (deferred) |

Open CP2K specs on `dev` = **0**. Stock VV10 is **blocked**. Incomplete
`[~]` = **0** after close of `gpw-fd-phonon-si` (2026-07-26; reg registered
with 1D-curvature dual-gate + UKS density warm-start). Pending residual
**12** (GapCensus wave 13/14 list); **not** queue exhausted.

## Closed — CP2K harness wave 12 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 12 | xl-bomd-linearized-forces, second-gen-cpmd-kpoints, gpw-d3-atm, gpw-d4-atm, bands-spin-metallic-si, dos-spin-metallic-si | zero-SCF linearized XL; multi-k SGCP; D3/D4 ATM; spin+metallic bands/DOS |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 11 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 11 | ugw-uks-rpbe, gpw-ts-forces, gpw-xdm-forces, gpw-kpoint-uks-stress-gga, cpmd-electronic-thermostat, xl-bomd-uks-kpoints, bands-metallic-si, dos-spin-si, dfpt-finite-q-spin-smearing | UGW@UKS-RPBE; TS/XDM force pins; multi-k UKS stress GGA; CPMD e-T; multi-k UKS XL; metallic bands; spin DOS; finite-q DFPT spin+smear |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 10 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 10 | dfpt-spin-smearing, ot-uks-dftu-smearing, second-gen-cpmd-uks-ot, xl-bomd-kpoints, bands-spin-si, dos-metallic-si, periodic-wannier-uks, berry-spin-si, gpw-mbd-kpoint-forces | DFPT spin+smear; multi-k UKS OT+U+smear; SGCP UKS OT; multi-k XL-BOMD; spin bands; metallic DOS; spin MLWF; spin Berry; multi-k MBD forces |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 9 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 9 | ot-uks-smearing, ot-uks-dftu, ot-dftu-smearing, gpw-d4-forces, gpw-mbd-forces, gpw-uks-multigrid, dfpt-finite-q-smearing, second-gen-cpmd-ot | UKS OT+smear; multi-k UKS OT+U; OT+U+smear; D4/MBD forces; UKS multigrid; finite-q DFPT smear; SGCP OT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 8 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 8 | periodic-urpa-kpoints, ot-smearing, ot-uks-kpoints, ot-dftu-kpoints, gpw-d4-stress, gpw-vdw-scf-potential, dfpt-metallic-si, dfpt-spin-incommensurate, gpw-optb88-vdw | multi-k URPA; OT+smear; multi-k UKS OT; multi-k OT+DFT+U; D4 analytic stress; nonlocal vdW SCF pot/forces; metallic DFPT; spin-incommensurate DFPT; optB88-vdW |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 7 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 7 | periodic-urpa-si, periodic-ugw-kpoints, gpw-uks-smearing, gpw-hirshfeld-uks, cpmd-uks, second-gen-cpmd-uks, ot-kpoints, ot-dftu, dfpt-spin-si | URPA@UKS; multi-k UGW; GPW UKS FD/MP smear; UKS Hirshfeld; CPMD UKS; SGCP UKS; OT+k; OT+DFT+U; collinear-spin DFPT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 6 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 6 | gpw-metallic-smearing, gpw-xdm-ice, gpw-uks-rpbe, gpw-dftu-kpoint-stress, periodic-ump2-kpoints, periodic-ugw-si, gpw-hirshfeld-kpoints, dfpt-gga-si, xl-bomd-uks | GPW FD/MP smear; free-atom XDM; Γ-UKS RPBE; multi-k DFT+U stress; multi-k UMP2; UGW@UKS; multi-k Hirshfeld; DFPT-GGA; XL-BOMD UKS |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 5 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 5 | gpw-dftu-kpoints, periodic-gw-finite-q, gpw-uks-revpbe, cpmd-rattle-water, periodic-ump2-si, ot-uks-o2 | multi-k PLUS_U; finite-q G0W0 \(W\); Γ-UKS revPBE; CPMD RATTLE; UMP2@UKS; OT-UKS O₂ |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 4 (2026-07-26)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 4 | gpw-kpoint-uks-stress, periodic-gw-kpoints, gpw-vdw-df2-ice, gpw-rvv10-ice, gpw-mbd-periodic, pw-spin-gga, gpw-uks-blyp, periodic-rpa-finite-q, gpw-dftu-stress | multi-k UKS stress; multi-k G0W0; vdW-DF2 LMKLL; rVV10; periodic Ewald MBD; PW collinear spin-GGA; Γ-UKS BLYP; finite-q dRPA; DFT+U analytic stress |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 3 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 3 | gpw-rpbe-water, gpw-kpoint-uks-forces, gpw-uks-stress, gpw-vdw-df-ice, gpw-mbd-hirshfeld, periodic-rpa-kpoints, periodic-gw-si, dfpt-incommensurate-si | RPBE GGA; multi-k UKS forces; Γ-UKS stress; vdW-DF DRSLL; Hirshfeld MBD/TS; multi-k dRPA; Γ spectral G0W0; incommensurate DFPT |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K harness wave 2 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 2 | gpw-kpoint-forces, gpw-revpbe-si, gpw-kpoints-uks, gpw-uks-forces, gpw-kpoint-stress, gpw-mbd-ice, periodic-mp2-kpoints, periodic-rpa-si, dfpt-finite-q-si, gpw-variable-cell-molopt | MP k-forces; revPBE GGA; k-UKS; Γ-UKS forces; k-stress; free-atom MBD; multi-k RI-MP2; Γ dRPA; finite-q DFPT X; MOLOPT CELL_OPT |

Do not re-create closed slugs without superseding intent.

MetalWalls wave 3 closed 2026-07-24 (specs deleted). **Wave-4 assign**
(2026-07-24) staffed 10 production-depth / residual map rows to
`## in_progress` (still unstarted as of **wave-5 gap-census** 2026-07-24:
MaZe, NeighborList salt/Gaussian, electrode forces, CV-SMD, aqueous EDL,
D-field, charge tools, …). Wave-5 census refilled **pending=10** (per-electrode
\(\eta\), const-charge control, MPI charge-halo, XL dipoles, QEq, PIM-NPT IL,
Voronoi TF, IL EDL box, nanoporous product, TI capacitance). No open
MD/electrochem specs on `dev` — workers write `.claude/specs/<slug>/` on
`wf/<slug>` before impl.

CP2K claude wave 1 closed 2026-07-25 (5 merged). **`xl-bomd-water` closed
2026-07-25** (resume of WIP — verified green, merged, spec deleted).


## Closed — CP2K harness wave 1 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop harness wave 1 | gpw-blyp-water, gpw-uks-pbe, pw-uks-fe-smearing, pw-variable-cell-si, gpw-ts-ice, phonon-dispersion-si, periodic-wannier-si, second-gen-cpmd-water, periodic-ri-mp2-si | GPW BLYP; Γ-UKS PBE O₂; bulk Fe UKS+smear; PW CELL_OPT Si; free-atom TS; phonon ω(q) IFC Fourier; periodic MV MLWF; 2nd-gen CP water; Γ RI-GPW MP2 product |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K parallel wave (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped |
| ----- | ------------------------ | ------- |
| cp2k-method-loop parallel | mpi-pw-scaling-4rank, gpw-analytic-stress, wannier-centers, mpi-gpw-scaling-4rank, cpmd-water | PW/GPW MPI k-RR gates; analytic GPW stress; Foster–Boys Wannier; CPMD H₂O |

## Closed — CP2K xl-bomd-water (2026-07-25)

| Queue | Spec (deleted; see git) | Shipped |
| ----- | ----------------------- | ------- |
| cp2k-method-loop | xl-bomd-water | `XlBomd` + `guess_P` + `GPWEvaluator::xl_bomd`; `reg_xl-bomd-water` NVE drift 1.2e-5 |

## Closed — CP2K capability claude wave 1 (2026-07-25)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop claude wave 1 | gpw-kpoints-si, gpw-variable-cell-ice, dfpt-phonon-si, berry-polarization, periodic-mp2-si | MP k-point GPW RKS; GPW variable-cell relax vs CP2K CELL_OPT; Γ-point DFPT phonons (PW); KSV Berry-phase polarization (PW); Γ-point canonical GPW MP2 foothold |

Do not re-create closed slugs without superseding intent.

## Closed — MetalWalls capability wave 3 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 3 | constant-force-pull-ion, oscillating-efield-md, potential-profile-observer, surface-tension-coupling, virtual-site-multi-topology, dipole-electrode-joint | FixAddForce; FixOscillatingEfield; Ψ(z) PotentialProfile; SurfaceTensionNPT; VirtualSites Weighted3+OutOfPlane3; DipoleElectrodeJoint |

Do not re-create closed slugs without superseding intent.

## Closed — MetalWalls capability wave 2 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 2 | const-pot-cg-large-electrode, electrode-charge-predictor, finite-field-single-electrode, pim-induced-dipole-bulk, piston-electrode-pz, thomas-fermi-electrode, fumi-tosi-nacl, steele-wall-argon, born-repulsion-foothold, capacitance-charge-fluctuation | ConstPotCG; ChargeASPC predictor; Dufils single-electrode finite-field; bulk PIM+Thole; PistonElectrode; TF metallicity; PairBorn/Fumi–Tosi; FixWallLJ1043 Steele; CapacitanceChargeFluctuation |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K capability wave 2 (2026-07-24)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop wave 2 | gpw-stress-water, gpw-bomd-water, pw-analytic-stress-si, fd-phonon-supercell, gpw-dftu-nio, bands-dos-si, pw-uks-kmesh-fe, gpw-pbesol-si, gpw-d4-water | GPW strain-FD stress; GPW BOMD NVE; PW Nielsen–Martin stress; FD supercell phonons; Dudarev DFT+U NiO; bands/DOS Si; k-mesh UKS O₂; PBEsol GPW; periodic D4 |

Do not re-create closed slugs without superseding intent.

## Closed — CP2K capability wave 1 (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cp2k-method-loop wave 1 | gpw-water-pbe-mt, gpw-forces-water, gpw-uks-o2-gamma, gpw-pbe-d3bj-ice, pw-lda-si, pw-pbe-gga-si, pw-smearing-al, pw-spin-polarized-o2, pw-bomd-h, ot-minimizer-water | MT water PBE pin; GPW forces; Γ-UKS O₂; GPW D3(BJ); PW LDA Si; true PBE GGA PW; Fermi/MP Al; PW spin O₂; PW BOMD H; OT water pin |

## Closed — MetalWalls capability wave 1 (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 1 | const-pot-bo-parallel-plate, freeze-electrode-atoms, const-pot-charge-observer, charging-work-energy, const-pot-md-nvt-short, gaussian-electrode-charge, slab-2d-ewald-foothold, external-efield-md, buckingham-pair-lammps, four-site-tip4p-foothold | ConstPotBO + FixFreeze + electrode charge observer + charging work + short const-pot NVT; Gaussian electrode Coulomb; slab Ewald (`slab_factor`); FixEfield; PairBuck; TIP4P `lj/cut/tip4p/cut` |

## Closed — CPU next campaign (2026-07-23, parallel worktrees)

| Branch (deleted) | Specs | Highlights |
| ------ | ----- | ---------- |
| `feat/cpu-next-md` | parrinello-rahman-npt, nph-ensemble | PR full-cell NPT, isotropic NPH |
| `feat/cpu-next-tb-arch` | tb-sk-product, concrete-pw-method, concrete-gpw-method, concrete-nddo-method | Multi-element DFTB2; concrete PW/GPW/AM1-RM1 |
| `feat/cpu-next-msamp` | qmmm-link-atoms, polarizable-qmmm, multi-cv-metad, multi-window-wham | Link atoms, pol QM/MM, ND metaD, WHAM (molrs) |
| `feat/cpu-next-response` | full-grid-alda, analytic-tda-forces, true-nac-vectors, ri-mp2 | Full MO-product ALDA, CIS TDA forces, NAC vectors, RI-MP2 |

## Closed earlier on `dev` (2026-07-23)

| depth-next | grid-alda-tda, qmmm-electrostatic, tda-nac-skeleton, triclinic-cell-relax | ALDA diagonal, ES QM/MM, NAC strength, triclinic cell |
| open scales | alda-tddft-tda … ri-mp2-skeleton (8) | GFN0, aniso cell, mech QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| method ladder | rks-tda … hamiltonian-remd-lambda (8) | RKS-TDA, DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD |
| breadth-next | soft-core-coulomb … brownian-dynamics (9) | soft Coul, AngleCV, TI, joint cell, string, IRC, Rg², SETTLE, BD |
| breadth-first | FIRE … dimer (6) | soft LJ, Widom, PIMD, isotropic cell, dimer |
| architecture | backend-mirror, cuda-resurrect, cuda-pme-accuracy | layout rewrite, CUDA live, PME Madelung |

## Closed — CUDA OpenMM MD engine (2026-07-23, from origin/dev)

| cuda-openmm-md 01–05 | contract, rebuild, tile-force, multigpu, delete | CS=32 clusters, real masks, tile LJ/Coul, multi-GPU classify; Verlet default removed on CUDA |
