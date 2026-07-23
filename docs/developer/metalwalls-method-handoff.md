# MetalWalls method loop — handoff

Last harness pass: **wave 1 — 2026-07-23** on `dev`
(`chore(harness): metalwalls wf wave 1`).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/metalwalls-method-loop.md`.

## Wave 1 summary

| Slug | Domain | Status | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ------ | ------ | ---------- | ------ | ---------------------------------------------- |
| `const-pot-bo-parallel-plate` | electrochem E1 | **merged → done** | `regressions/md/const-pot-bo-parallel-plate` | lammps | **perf_exempt** |
| `freeze-electrode-atoms` | md E1 | **merged → done** | `regressions/md/freeze-electrode-atoms` | lammps | **3.98** |
| `const-pot-charge-observer` | electrochem E1 | **merged → done** | `regressions/md/const-pot-charge-observer` | lammps | **perf_exempt** |
| `charging-work-energy` | electrochem E1 | **merged → done** | `regressions/md/charging-work-energy` | literature | **perf_exempt** |
| `const-pot-md-nvt-short` | electrochem E1 | **merged → done** | `regressions/md/const-pot-md-nvt-short` | metalwalls | **perf_exempt** |
| `gaussian-electrode-charge` | electrochem E2 | **merged → done** | `regressions/md/gaussian-electrode-charge` | literature | **perf_exempt** |
| `slab-2d-ewald-foothold` | electrochem E2 | **merged → done** | `regressions/md/slab-2d-ewald-foothold` | lammps | **6.49** |
| `external-efield-md` | electrochem E3 | **merged → done** | `regressions/md/external-efield-md` | lammps | **0.60** |
| `buckingham-pair-lammps` | forcefield E8 | **merged → done** | `regressions/forcefield/buckingham-pair-lammps` | lammps | **perf_exempt** |
| `four-site-tip4p-foothold` | md E9 | **merged → done** | `regressions/md/four-site-tip4p-foothold` | lammps | **0.338** |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Incomplete remaining:** **0** (no open `[~]`, no open MD/electrochem specs)
- **Queue after harness:** pending = **14**; in_progress = **0**; done = **10** (wave 1); blocked = **4**
- **Roadmap (CPU ○→✓ on `dev`):** Buckingham pair; constant electric-field MD
  (`FixEfield`); TIP4P-class virtual M-site (pair style foothold — not general
  multi-topology virtual-site dynamics product). E1 const-pot stack ships as
  Components on `cpu::md::MD` (`ConstPotBO`, `FixFreeze`, observers, charging
  work); slab Ewald via `KSpaceEwald::slab_factor`.
- **Specs closed:** all 10 wave-1 slugs (acceptance met on `dev`; deleted under
  `.claude/specs/`).
- **Capability map:** `.claude/notes/metalwalls-capability-map.md` (E1–E12) —
  no lasting rule rewrite this harness.
- **Orphans (not incomplete for this loop):**
  - `regressions/md/{dpd_ab_mix,dpd_groot_warren,fene_wca_chain}/`,
    `regressions/forcefield/yukawa_dimer/` — method-ladder shipped reg bodies,
    no `REFERENCE.md` (hygiene; not on MW queue; do not invent `[~]` rows)
  - local `wf/*` tips for wave-1 slugs fully merged into `dev` (stale branch names OK)

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin (or literature when no binary)
- `ctest` never launches LAMMPS or MetalWalls

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## Next assign candidates (pending head after wave 1)

1. `const-pot-cg-large-electrode` (E2, oracle=metalwalls; dep BO dense ✓)
2. `electrode-charge-predictor` (E2, oracle=metalwalls)
3. `finite-field-single-electrode` (E3, oracle=metalwalls; dep E1+efield ✓)
4. `pim-induced-dipole-bulk` (E4, oracle=metalwalls)
5. `dipole-electrode-joint` (E5; dep E1+PIM)
6. `piston-electrode-pz` / `thomas-fermi-electrode` (E6/E7)
7. `fumi-tosi-nacl` / `steele-wall-argon` / `born-repulsion-foothold` (E8; Buck ✓)
8. `capacitance-charge-fluctuation` / `potential-profile-observer` (E10; dep E1 ✓)
9. `surface-tension-coupling` (shared classical; oracle=lammps)
10. `constant-force-pull-ion` (shared classical; oracle=lammps)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only.
Default **staff_count=10**. Spec-first: write `.claude/specs/<slug>/` **before**
code on each worker branch.
