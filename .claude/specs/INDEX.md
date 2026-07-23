# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

## Live

| Spec | Branch | Status |
| ---- | ------ | ------ |
| [gpw-forces-water](gpw-forces-water/) | `wf/gpw-forces-water` | open — CP2K-loop wave 1 merge |
| [gpw-pbe-d3bj-ice](gpw-pbe-d3bj-ice/) | `wf/gpw-pbe-d3bj-ice` | open — CP2K-loop wave 1 merge |
| [gpw-uks-o2-gamma](gpw-uks-o2-gamma/) | `wf/gpw-uks-o2-gamma` | open — CP2K-loop wave 1 merge |
| [gpw-water-pbe-mt](gpw-water-pbe-mt/) | `wf/gpw-water-pbe-mt` | open — CP2K-loop wave 1 merge |
| [ot-minimizer-water](ot-minimizer-water/) | `wf/ot-minimizer-water` | open — CP2K-loop wave 1 merge |
| [pw-bomd-h](pw-bomd-h/) | `wf/pw-bomd-h` | open — CP2K-loop wave 1 merge |
| [pw-lda-si](pw-lda-si/) | `wf/pw-lda-si` | open — CP2K-loop wave 1 merge |
| [pw-pbe-gga-si](pw-pbe-gga-si/) | `wf/pw-pbe-gga-si` | open — CP2K-loop wave 1 merge |
| [pw-spin-polarized-o2](pw-spin-polarized-o2/) | `wf/pw-spin-polarized-o2` | open — CP2K-loop wave 1 merge |

## Closed — CPU next campaign (2026-07-23, parallel worktrees)

| Branch (deleted) | Specs | Highlights |
| ------ | ----- | ---------- |
| `feat/cpu-next-md` | parrinello-rahman-npt, nph-ensemble | PR full-cell NPT, isotropic NPH |
| `feat/cpu-next-tb-arch` | tb-sk-product, concrete-pw-method, concrete-gpw-method, concrete-nddo-method | Multi-element DFTB2; concrete PW/GPW/AM1-RM1 |
| `feat/cpu-next-msamp` | qmmm-link-atoms, polarizable-qmmm, multi-cv-metad, multi-window-wham | Link atoms, pol QM/MM, ND metaD, WHAM (molrs) |
| `feat/cpu-next-response` | full-grid-alda, analytic-tda-forces, true-nac-vectors, ri-mp2 | Full MO-product ALDA, CIS TDA forces, NAC vectors, RI-MP2 |

**Hygiene (2026-07-23):** stale Grok clones under
`~/.grok/worktrees/molcrafts-atomiverse/{cpu-deps,fix-cuda-pme,impl-kokkos}`
removed. Unique Kokkos WIP preserved only as local branch
`archive/impl-kokkos-wip` (not merged to `dev`). Roadmap resynced in
`docs/roadmap/`.

## Closed earlier on `dev` (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| depth-next | grid-alda-tda, qmmm-electrostatic, tda-nac-skeleton, triclinic-cell-relax | ALDA diagonal, ES QM/MM, NAC strength, triclinic cell |
| open scales | alda-tddft-tda … ri-mp2-skeleton (8) | GFN0, aniso cell, mech QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| method ladder | rks-tda … hamiltonian-remd-lambda (8) | RKS-TDA, DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD |
| breadth-next | soft-core-coulomb … brownian-dynamics (9) | soft Coul, AngleCV, TI, joint cell, string, IRC, Rg², SETTLE, BD |
| breadth-first | FIRE … dimer (6) | soft LJ, Widom, PIMD, isotropic cell, dimer |
| architecture | backend-mirror, cuda-resurrect, cuda-pme-accuracy | layout rewrite, CUDA live, PME Madelung |

Do not re-create closed slugs without superseding intent.

## Closed — CUDA OpenMM MD engine (2026-07-23, from origin/dev)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| cuda-openmm-md 01–05 | contract, rebuild, tile-force, multigpu, delete | CS=32 clusters, real masks, tile LJ/Coul, multi-GPU classify; Verlet default removed on CUDA |

Do not re-create closed slugs without superseding intent.
