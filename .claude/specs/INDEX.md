# Specs index
Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.
## Live
*(none — MetalWalls wave 1 closed 2026-07-23; next workers re-create specs on
`wf/<slug>` before code)*
## Closed — MetalWalls capability wave 1 (2026-07-23)
| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| metalwalls-method-loop wave 1 | const-pot-bo-parallel-plate, freeze-electrode-atoms, const-pot-charge-observer, charging-work-energy, const-pot-md-nvt-short, gaussian-electrode-charge, slab-2d-ewald-foothold, external-efield-md, buckingham-pair-lammps, four-site-tip4p-foothold | ConstPotBO + FixFreeze + electrode charge observer + charging work + short const-pot NVT; Gaussian electrode Coulomb; slab Ewald (`slab_factor`); FixEfield; PairBuck; TIP4P `lj/cut/tip4p/cut` |
Do not re-create closed slugs without superseding intent.
## Closed — CP2K capability wave 1 (2026-07-23)
| cp2k-method-loop wave 1 | gpw-water-pbe-mt, gpw-forces-water, gpw-uks-o2-gamma, gpw-pbe-d3bj-ice, pw-lda-si, pw-pbe-gga-si, pw-smearing-al, pw-spin-polarized-o2, pw-bomd-h, ot-minimizer-water | MT water PBE pin; GPW forces; Γ-UKS O₂; GPW D3(BJ); PW LDA Si; true PBE GGA PW; Fermi/MP Al; PW spin O₂; PW BOMD H; OT water pin |
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
| Spec | Branch | Status |
| ---- | ------ | ------ |
| [gpw-stress-water](gpw-stress-water/) | `wf/gpw-stress-water` | open — GPW cell stress (water) |
| [buckingham-pair-lammps](buckingham-pair-lammps/) | `wf/buckingham-pair-lammps` | open — buckingham-pair-lammps |
| [charging-work-energy](charging-work-energy/) | `wf/charging-work-energy` | open — charging-work-energy |
| [const-pot-bo-parallel-plate](const-pot-bo-parallel-plate/) | `wf/const-pot-bo-parallel-plate` | open — Spec: const-pot-bo-parallel-plate |
| [const-pot-charge-observer](const-pot-charge-observer/) | `wf/const-pot-charge-observer` | open — Spec: `const-pot-charge-observer` |
| [const-pot-md-nvt-short](const-pot-md-nvt-short/) | `wf/const-pot-md-nvt-short` | open — Spec: `const-pot-md-nvt-short` |
| [external-efield-md](external-efield-md/) | `wf/external-efield-md` | open — Spec: `external-efield-md` |
| [four-site-tip4p-foothold](four-site-tip4p-foothold/) | `wf/four-site-tip4p-foothold` | open — four-site-tip4p-foothold |
| [freeze-electrode-atoms](freeze-electrode-atoms/) | `wf/freeze-electrode-atoms` | open — Spec: `freeze-electrode-atoms` |
| [gaussian-electrode-charge](gaussian-electrode-charge/) | `wf/gaussian-electrode-charge` | open — gaussian-electrode-charge |
| [slab-2d-ewald-foothold](slab-2d-ewald-foothold/) | `wf/slab-2d-ewald-foothold` | open — Spec: slab-2d-ewald-foothold |
| [gpw-bomd-water](gpw-bomd-water/) | `wf/gpw-bomd-water` | open — GPW BOMD water NVE energy-drift pin |
| [pw-analytic-stress-si](pw-analytic-stress-si/) | `wf/pw-analytic-stress-si` | open — PW Nielsen–Martin analytic stress (Si NC-UPF) |
| [fd-phonon-supercell](fd-phonon-supercell/) | `wf/fd-phonon-supercell` | open — FD supercell phonons (Γ) |
| [bands-dos-si](bands-dos-si/) | `wf/bands-dos-si` | open — Band structure / DOS PW Si |
| [gpw-dftu-nio](gpw-dftu-nio/) | `wf/gpw-dftu-nio` | open — Periodic DFT+U (NiO mol PLUS_U pin) |
