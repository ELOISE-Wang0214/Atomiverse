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

## Live

| Spec | Branch / track | Status |
| ---- | -------------- | ------ |
| [periodic-rhf-gamma](periodic-rhf-gamma/) | pyscf / HF periodic | open (deferred) |
| [rsh-eri-cam-b3lyp](rsh-eri-cam-b3lyp/) | pyscf / RSH | open (deferred) |

MetalWalls wave 3 closed 2026-07-24 (specs deleted). No open MD/electrochem
specs on `dev`. Classical MW staff queue exhausted (pending empty).

CP2K wave 2 closed 2026-07-24 (specs deleted). Incomplete leftover
`gpw-kpoints-si` has **no** open spec dir — write-first on resume.

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
