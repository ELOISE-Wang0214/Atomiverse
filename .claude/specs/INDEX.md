# Specs index

Live specs under `.claude/specs/`. Completed work is deleted; history is the
git log. Method census lives in `docs/roadmap/`.

## Live

*(none — CPU breadth-next queue closed 2026-07-23)*

## Closed (2026-07-23) — CPU breadth-next textbook queue

| # | Spec | Shipped as |
| - | ---- | ---------- |
| 1 | ~~soft-core-coulomb~~ | `cpu::PairCoulSoft`, `reg_coul_soft_dimer` |
| 2 | ~~angle-cv~~ | `cpu::bias::AngleCV`, `reg_umbrella_angle` |
| 3 | ~~fep-soft-dimer~~ | `PairLJSoft::dU_dlambda`, `reg_fep_soft_dimer` |
| 4 | ~~joint-atom-cell~~ | `JointIsotropicRelax`, `reg_joint_lj_zero_p` |
| 5 | ~~string-method~~ | `cpu::path::StringMethod`, `reg_string_muller_brown` |
| 6 | ~~irc-path~~ | `cpu::path::Irc`, `reg_irc_muller_brown` |
| 7 | ~~rpmd-estimator~~ | `Pimd::radius_of_gyration_sq` / analytic, `reg_pimd_free_particle_estimators` |
| 8 | ~~settle-water~~ | `cpu::md::FixSettle`, `reg_settle_tip3p` |
| 9 | ~~brownian-dynamics~~ | `cpu::md::Brownian`, `reg_brownian_einstein` |

Prior closed: CPU breadth 1–6 (FIRE…Dimer), backend-mirror, cuda-resurrect,
cuda-pme-accuracy, HF fossils purged.
