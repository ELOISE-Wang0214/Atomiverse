# MetalWalls method loop — handoff

Last harness pass: **recover 2026-07-23** (workflow seeded same day; no prior wave).

## Wave summary

| Slug | Domain | Status | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ------ | ------ | ---------- | ------ | ---------------------------------------------- |
| — | — | — | — | — | — |

- **Blocked this wave:** —
- **Merge-failed this wave:** —
- **Queue:** see `docs/developer/metalwalls-method-queue.md` (seed pending; in_progress empty)
- **Incomplete remaining:** **0** (recover reconfirm 2026-07-23)
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

## Next assign candidates (pending head)

1. `const-pot-bo-parallel-plate` (E1, oracle=metalwalls / ELECTRODE)
2. `freeze-electrode-atoms` (oracle=lammps)
3. `coul-long-pme-nacl-lammps` (oracle=lammps)
4. `nose-hoover-nvt-lammps` (oracle=lammps)
5. `special-bonds-amber-lammps` (oracle=lammps)
6. `lj-cut-lammps-dimer` (oracle=lammps; skip if already golden under another reg)
7. `const-pot-charge-observer`
8. `charging-work-energy`
9. `const-pot-md-nvt-short`
10. `external-efield-md`

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.
Default **staff_count=10**.
