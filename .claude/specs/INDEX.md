# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

## Live — CPU next campaign (2026-07-23)

### A/B — Response (merging next: `feat/cpu-next-response`)

| Spec | Summary |
| ---- | ------- |
| [full-grid-alda](full-grid-alda.md) | Full grid MO-product ALDA for RKS TDA [approved] |
| [analytic-tda-forces](analytic-tda-forces.md) | Analytic CIS/TDA nuclear forces on RHF [approved] |
| [true-nac-vectors](true-nac-vectors.md) | True TDA NAC vectors (not strength-only) [approved] |
| [ri-mp2](ri-mp2.md) | Density-fitted RI-MP2 energy on RHF [approved] |

### Closed this campaign (worktrees)

| Branch | Specs | Highlights |
| ------ | ----- | ---------- |
| `feat/cpu-next-md` | parrinello-rahman-npt, nph-ensemble | PR full-cell NPT, isotropic NPH |
| `feat/cpu-next-tb-arch` | tb-sk-product, concrete-pw/gpw/nddo | Multi-element DFTB2; concrete methods |
| `feat/cpu-next-msamp` | qmmm-link-atoms, polarizable-qmmm, multi-cv-metad, multi-window-wham | Link atoms, pol QM/MM, ND metaD, WHAM |

## Closed on `dev` (earlier 2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| depth-next | grid-alda-tda, qmmm-electrostatic, tda-nac-skeleton, triclinic-cell-relax | ALDA diagonal, ES QM/MM, NAC strength, triclinic cell |
| open scales | alda-tddft-tda … ri-mp2-skeleton (8) | GFN0, aniso cell, mech QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| method ladder | rks-tda … hamiltonian-remd-lambda (8) | RKS-TDA, DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD |
| breadth-next | soft-core-coulomb … brownian-dynamics (9) | soft Coul, AngleCV, TI, joint cell, string, IRC, Rg², SETTLE, BD |
| breadth-first | FIRE … dimer (6) | soft LJ, Widom, PIMD, isotropic cell, dimer |
| architecture | backend-mirror, cuda-resurrect, cuda-pme-accuracy | layout rewrite, CUDA live, PME Madelung |

Do not re-create closed slugs without superseding intent.
