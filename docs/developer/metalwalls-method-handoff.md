# MetalWalls method loop — handoff

Last harness: **wave 2 — 2026-07-24** on `dev`
(`chore(harness): metalwalls wf wave 2`). Last assign: **wave 3 — 2026-07-24**
(`chore(queue): metalwalls assign wave 3`). Last recover: **2026-07-24** on `dev`
(`chore(queue): metalwalls recover incomplete`).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/metalwalls-method-loop.md`.

## Recover inventory (pre assign wave 3, 2026-07-24)

| Source | Result |
| ------ | ------ |
| A queue `## in_progress` `[~]` | **0** (before assign) |
| B open `.claude/specs/` MD/electrochem | **0** (wave-2 specs deleted; PySCF deferred specs ignored) |
| C unmerged `wf/<slug>` (this loop) | **0** (wave-1/2 tips `ahead_of_dev=0`) |
| D handoff incomplete rows | **none** |
| E half-written `regressions/{md,forcefield}/` wave-3 | **0** |

**Resume set:** empty. **Assign wave 3** promoted all **6** remaining pending
rows → `## in_progress` `[~]` (`is_resume=false`). Staff default 10; only 6
open rows existed. **incomplete_remaining=6** after assign.

## Wave 1 summary (historical — harnessed 2026-07-23)

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

## Wave 2 summary (harnessed 2026-07-24)

| Slug | Domain | Status | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ------ | ------ | ---------- | ------ | ---------------------------------------------- |
| `const-pot-cg-large-electrode` | electrochem E2 | **merged → done** | `regressions/md/const-pot-cg-large-electrode` | metalwalls | **perf_exempt** |
| `electrode-charge-predictor` | electrochem E2 | **merged → done** | `regressions/md/electrode-charge-predictor` | metalwalls | **perf_exempt** |
| `finite-field-single-electrode` | electrochem E3 | **merged → done** | `regressions/md/finite-field-single-electrode` | metalwalls | **perf_exempt** |
| `pim-induced-dipole-bulk` | polarizable E4 | **merged → done** | `regressions/md/pim-induced-dipole-bulk` | metalwalls | **perf_exempt** |
| `piston-electrode-pz` | electrochem E6 | **merged → done** | `regressions/md/piston-electrode-pz` | metalwalls | **perf_exempt** |
| `thomas-fermi-electrode` | electrochem E7 | **merged → done** | `regressions/md/thomas-fermi-electrode` | metalwalls | **perf_exempt** |
| `fumi-tosi-nacl` | forcefield E8 | **merged → done** | `regressions/forcefield/fumi-tosi-nacl` | lammps | **perf_exempt** |
| `steele-wall-argon` | forcefield E8 | **merged → done** | `regressions/forcefield/steele-wall-argon` | lammps | **perf_exempt** |
| `born-repulsion-foothold` | forcefield E8 | **merged → done** | `regressions/forcefield/born-repulsion-foothold` | lammps | **perf_exempt** |
| `capacitance-charge-fluctuation` | electrochem E10 | **merged → done** | `regressions/md/capacitance-charge-fluctuation` | literature | **perf_exempt** |

- **Blocked this wave:** none
- **Merge-failed this wave:** none
- **Incomplete remaining after harness wave 2:** **0**
- **Queue after harness wave 2:** pending = **6**; in_progress = **0**; done = **20** (wave 1 + wave 2); blocked = **4**
- **Unblocked by wave 2:** `dipole-electrode-joint` (PIM dep satisfied)
- **Roadmap (CPU ○→✓ on `dev` after wave 2):** Born–Mayer–Huggins pair
  (`PairBorn`); multi-type Fumi–Tosi Born+Coul path; Steele 10-4-3 wall
  (`FixWallLJ1043`); bulk PIM self-consistent induced dipoles + Thole (TT
  optional). Electrode stack additions: `ConstPotCG`, ASPC charge predictor,
  single-electrode finite-field + `ClassicalPolarization`, `PistonElectrode`,
  Thomas–Fermi hardness on `ConstPotBO`, `CapacitanceChargeFluctuation`
  observer — all Components on `cpu::md::MD` / Evaluate (not new Driver
  categories; charge solve not inside `KSpace*`).
- **Specs closed:** all 10 wave-2 slugs (acceptance met on `dev`; deleted under
  `.claude/specs/`).
- **Capability map:** `.claude/notes/metalwalls-capability-map.md` (E1–E12) —
  path referenced by loop; recreate if missing on worker trees. No lasting
  rule change this harness (Steele oracle note remains: LAMMPS `wall/lj1043`
  not `wall/lj93`).
- **Orphans (not incomplete for this loop):**
  - `regressions/md/{dpd_ab_mix,dpd_groot_warren,fene_wca_chain}/`,
    `regressions/forcefield/yukawa_dimer/` — method-ladder shipped reg bodies,
    no `REFERENCE.md` (hygiene; not on MW queue; do not invent `[~]` rows)
  - local `wf/*` tips for wave-1/2 slugs fully merged into `dev` (stale branch
    names OK)

## Assign wave 3 (2026-07-24) — staffed in_progress

| Slug | Domain | oracle | is_resume | Priority note |
| ---- | ------ | ------ | --------- | ------------- |
| `dipole-electrode-joint` | electrochem E5 | metalwalls | false | **dep met** (PIM ✓); JCP 2022 main hardness |
| `potential-profile-observer` | electrochem E10 | metalwalls | false | Poisson \(\Psi(z)\); literature dual OK |
| `oscillating-efield-md` | electrochem E3 | lammps | false | AC residual after constant efield ✓ |
| `virtual-site-multi-topology` | md E9 | lammps | false | beyond TIP4P M-site foothold |
| `constant-force-pull-ion` | md (shared) | lammps | false | EDL transport twin; not MW-only |
| `surface-tension-coupling` | md (shared) | lammps | false | planar interface ensemble |

- **Assigned:** **6** (all remaining pending; staff default 10)
- **Resumed:** **0**
- **Blocked skipped:** `cuda-const-pot`, `aim-daim-product`, `ml-pinnwall-electrode`, `periodic-hybrid-electrode-dft`
- **Queue after assign wave 3:** pending = **0**; in_progress = **6**; done = **20**; blocked = **4**
- **incomplete_remaining:** **6**
- **Workers:** create `wf/<slug>` under
  `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/grok-worktrees/`;
  **spec-first** (`.claude/specs/<slug>/` + acceptance) before code; MD-only build

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin (or literature when no binary)
- `ctest` never launches LAMMPS or MetalWalls
- Wave-2 oracle tags: CG/predictor/finite-field/PIM/piston/TF → **metalwalls**;
  capacitance → **literature**; Fumi–Tosi / Steele / Born → **lammps**
  (`born/coul/*`, `wall/lj1043`, `pair born`)
- Wave-3 oracle tags: dipole⊕electrode joint + \(\Psi(z)\) → **metalwalls**;
  oscillating efield / virtual sites / force-pull / surface tension → **lammps**

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## Next (wave 3 workers — no further pending)

1. `dipole-electrode-joint` (E5, oracle=metalwalls; **dep met:** PIM ✓)
2. `potential-profile-observer` (E10, oracle=metalwalls / literature Poisson)
3. `oscillating-efield-md` (E3 residual, oracle=lammps)
4. `virtual-site-multi-topology` (E9 residual, oracle=lammps)
5. `constant-force-pull-ion` (shared classical; oracle=lammps)
6. `surface-tension-coupling` (shared classical; oracle=lammps)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only
under `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/grok-worktrees/`
(never bulk under `/home/jicli594`). Default **staff_count=10**. Spec-first:
write `.claude/specs/<slug>/` **before** code on each worker branch.
