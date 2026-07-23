# Specs index

Live specs under `.claude/specs/`. Completed work is deleted; history is the
git log. Method census lives in `docs/roadmap/`.

## Live

*(none — CPU open-scale queue closed 2026-07-23)*

## Closed (2026-07-23) — CPU open scales

| # | Spec | Shipped as |
| - | ---- | ---------- |
| 1 | ~~alda-tddft-tda~~ | `RKS::tda` `alda` flag (LSDA path), `reg_h2_rks_alda_tda` |
| 2 | ~~xtb-gfn-skeleton~~ | `cpu::semiempirical::XtbGfn0`, `reg_h2o_gfn0` |
| 3 | ~~anisotropic-cell~~ | `AnisotropicCellRelax`, `reg_aniso_lj_box` |
| 4 | ~~qmmm-mechanical~~ | `cpu::QmmmMechanical`, `reg_qmmm_h2_ar` |
| 5 | ~~metadynamics-1d~~ | `cpu::bias::Metadynamics1D`, `reg_metad_double_well` |
| 6 | ~~dpd-multicomponent~~ | `PairDpd` multi-type `a_matrix`, `reg_dpd_ab_mix` |
| 7 | ~~tda-excited-forces~~ | `RHF::tda_forces` FD, `reg_h2_tda_force` |
| 8 | ~~ri-mp2-skeleton~~ | public `RHF::mp2` regression `reg_h2o_mp2_sto3g` |

## Closed earlier

| Queue | Items |
| ----- | ----- |
| method ladder 1–8 | RKS-TDA … H-REMD λ |
| breadth-next 1–9 | soft Coulomb … Brownian |
| breadth-first 1–6 | FIRE … Dimer |
