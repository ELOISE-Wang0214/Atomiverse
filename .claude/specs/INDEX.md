# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

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

## Live — deferred PySCF only (CP2K clean after bootstrap 2026-07-25)

| Spec | Branch / track | Status |
| ---- | -------------- | ------ |
| [periodic-rhf-gamma](periodic-rhf-gamma/) | pyscf / HF periodic | open (deferred) |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/) | pyscf / RSH | open (deferred) |

No open CP2K/GPW/PW specs on `dev` — workers write `.claude/specs/<slug>/` on
`wf/<slug>` before impl. Parent bootstrap GapCensus 2026-07-25: pending=14,
in_progress=0.

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
