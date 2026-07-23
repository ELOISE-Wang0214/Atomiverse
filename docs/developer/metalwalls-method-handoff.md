# MetalWalls method loop — handoff

Last harness pass: **assign wave 1 — 2026-07-23** (after recover + census same day).

## Wave summary

| Slug | Domain | Status | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ------ | ------ | ---------- | ------ | ---------------------------------------------- |
| `const-pot-bo-parallel-plate` | electrochem E1 | [~] assigned | — | metalwalls | — |
| `freeze-electrode-atoms` | md E1 | [~] assigned | — | lammps | — |
| `const-pot-charge-observer` | electrochem E1 | [~] assigned | — | lammps | — |
| `charging-work-energy` | electrochem E1 | [~] assigned | — | metalwalls | — |
| `const-pot-md-nvt-short` | electrochem E1 | [~] assigned | — | metalwalls | — |
| `gaussian-electrode-charge` | electrochem E2 | [~] assigned | — | metalwalls | — |
| `slab-2d-ewald-foothold` | electrochem E2 | [~] assigned | — | lammps | — |
| `external-efield-md` | electrochem E3 | [~] assigned | — | lammps | — |
| `buckingham-pair-lammps` | forcefield E8 | [~] assigned | — | lammps | — |
| `four-site-tip4p-foothold` | md E9 | [~] assigned | — | lammps | — |

- **Blocked this wave:** —
- **Merge-failed this wave:** —
- **Queue:** `docs/developer/metalwalls-method-queue.md` — **10** `[~]` in_progress; pending tail remains
- **Incomplete remaining at assign:** **0** resumes; **10** newly assigned (all `is_resume=false`)
- **Capability map:** `.claude/notes/metalwalls-capability-map.md` (E1–E12)
- **Orphans (not incomplete for this loop):**
  - `regressions/md/{dpd_ab_mix,dpd_groot_warren,fene_wca_chain}/`,
    `regressions/forcefield/yukawa_dimer/` — method-ladder shipped reg bodies,
    no `REFERENCE.md` (hygiene; not on MW queue; do not invent `[~]` rows)
  - local `wf/*` PySCF/post-HF branch tips fully merged into `dev` (stale names)

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin
- `ctest` never launches LAMMPS or MetalWalls

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## Next assign candidates (pending head after wave 1)

1. `const-pot-cg-large-electrode` (E2, oracle=metalwalls; dep BO dense)
2. `electrode-charge-predictor` (E2, oracle=metalwalls)
3. `finite-field-single-electrode` (E3, oracle=metalwalls; dep E1+efield)
4. `pim-induced-dipole-bulk` (E4, oracle=metalwalls)
5. `dipole-electrode-joint` (E5; dep E1+PIM)
6. `piston-electrode-pz` / `thomas-fermi-electrode` (E6/E7)
7. `fumi-tosi-nacl` / `steele-wall-argon` / `born-repulsion-foothold` (E8)
8. `capacitance-charge-fluctuation` / `potential-profile-observer` (E10; dep E1)
9. `surface-tension-coupling` (shared classical; oracle=lammps)
10. `constant-force-pull-ion` (still pending if not taken; oracle=lammps)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.
Default **staff_count=10**.
