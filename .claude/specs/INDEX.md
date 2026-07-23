# Specs index

Live specs under `.claude/specs/`. Completed work is **deleted**; history is
the git log. Method census: `docs/roadmap/`.

## Live — CPU next campaign (2026-07-23)

Parallel tracks **A** (roadmap recommended next) + **B** (high-leverage gaps,
folded into the same slugs) + **C** (architecture debt). All `approved`;
impl independently (no single ordered chain unless noted).

### A/B — Response

| Spec | Summary |
| ---- | ------- |
| [full-grid-alda](full-grid-alda.md) | Full grid MO-product ALDA for RKS TDA [approved] |
| [analytic-tda-forces](analytic-tda-forces.md) | Analytic CIS/TDA nuclear forces on RHF [approved] |
| [true-nac-vectors](true-nac-vectors.md) | True TDA NAC vectors (not strength-only) [approved] |

### A/B — Cell / MD

| Spec | Summary |
| ---- | ------- |
| [parrinello-rahman-npt](parrinello-rahman-npt.md) | Full-cell PR NPT integrator [approved] |
| [nph-ensemble](nph-ensemble.md) | Thermostat-free NPH path [approved] |

### A/B — Multiscale

| Spec | Summary |
| ---- | ------- |
| [qmmm-link-atoms](qmmm-link-atoms.md) | H link-atom QM/MM boundary [approved] |
| [polarizable-qmmm](polarizable-qmmm.md) | Induced-dipole polarizable QM/MM [approved] |

### A/B — Post-HF / TB / Sampling

| Spec | Summary |
| ---- | ------- |
| [ri-mp2](ri-mp2.md) | Density-fitted RI-MP2 energy on RHF [approved] |
| [multi-cv-metad](multi-cv-metad.md) | Multi-CV well-tempered metadynamics [approved] |
| [multi-window-wham](multi-window-wham.md) | Multi-window WHAM (math in molrs) [approved] |

### C — Architecture (concrete methods)

| Spec | Summary |
| ---- | ------- |

**Suggested parallel impl batches** (no hard deps inside a batch):

3. polarizable-qmmm · multi-window-wham · concrete-nddo · true-nac (after analytic forces if sharing Z-vector)

Soft deps: `true-nac-vectors` benefits from `analytic-tda-forces`; `nph-ensemble` may reuse PR cell math from `parrinello-rahman-npt` but isotropic v1 is standalone.

## Closed on `dev` (2026-07-23)

| Queue | Specs (deleted; see git) | Shipped highlights |
| ----- | ------------------------ | ------------------ |
| depth-next | grid-alda-tda, qmmm-electrostatic, tda-nac-skeleton, triclinic-cell-relax | ALDA diagonal, ES QM/MM, NAC strength, triclinic cell |
| open scales | alda-tddft-tda … ri-mp2-skeleton (8) | GFN0, aniso cell, mech QM/MM, metaD, multi-DPD, TDA forces, MP2 |
| method ladder | rks-tda … hamiltonian-remd-lambda (8) | RKS-TDA, DFTB2, PW-cell, DihedralCV, FENE+WCA, DPD, Yukawa, H-REMD |
| breadth-next | soft-core-coulomb … brownian-dynamics (9) | soft Coul, AngleCV, TI, joint cell, string, IRC, Rg², SETTLE, BD |
| breadth-first | FIRE … dimer (6) | soft LJ, Widom, PIMD, isotropic cell, dimer |
| architecture | backend-mirror, cuda-resurrect, cuda-pme-accuracy | layout rewrite, CUDA live, PME Madelung |

Do not re-create closed slugs without superseding intent.
