# MetalWalls method loop — handoff

Last assign pass: **wave 2 — 2026-07-24** on `dev`
(`chore(queue): metalwalls assign wave 2`). Last harness: **wave 1 — 2026-07-23**.
Last recover: **2026-07-24** on `dev` (`chore(queue): metalwalls recover incomplete`).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/metalwalls-method-loop.md`.

## Recover inventory (2026-07-24)

| Source | Result |
| ------ | ------ |
| A queue `## in_progress` `[~]` | **10** wave-2 staff batch (unchanged; already `[~]`) |
| B open `.claude/specs/` MD/electrochem | **0** (workers must create on `wf/<slug>`) |
| C unmerged `wf/<slug>` (this loop) | **0** (no wave-2 branches; wave-1 tips `ahead_of_dev=0`) |
| D handoff incomplete rows | wave-2 list only (not merged) |
| E half-written `regressions/{md,forcefield}/` wave-2 | **0** |

**Resume set (staff=10, incomplete ≫ pending):** all 10 wave-2 `[~]` rows.
No demotion; no brand-new pending. Spec-first on each `wf/<slug>` under
`/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/grok-worktrees/`.

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
- **Incomplete remaining (pre–wave-2 assign):** **0** (no open `[~]`, no open MD/electrochem specs)
- **Queue after wave-2 assign:** pending = **4**; in_progress = **10**; done = **10** (wave 1); blocked = **4**
- **Wave 2 in_progress `[~]` (is_resume=false):** `const-pot-cg-large-electrode`,
  `electrode-charge-predictor`, `finite-field-single-electrode`,
  `pim-induced-dipole-bulk`, `piston-electrode-pz`, `thomas-fermi-electrode`,
  `fumi-tosi-nacl`, `steele-wall-argon`, `born-repulsion-foothold`,
  `capacitance-charge-fluctuation`
- **Held for deps (still pending):** `dipole-electrode-joint` (needs PIM)
- **Still pending (slots after staff=10):** `potential-profile-observer`,
  `constant-force-pull-ion`, `surface-tension-coupling`
- **Roadmap (CPU ○→✓ on `dev` after wave 1):** Buckingham pair; constant
  electric-field MD (`FixEfield`); TIP4P-class virtual M-site (pair style
  foothold — not general multi-topology virtual-site dynamics product). E1
  const-pot stack ships as Components on `cpu::md::MD` (`ConstPotBO`,
  `FixFreeze`, observers, charging work); slab Ewald via
  `KSpaceEwald::slab_factor`.
- **Specs closed:** all 10 wave-1 slugs (acceptance met on `dev`; deleted under
  `.claude/specs/`). Wave-2 workers **create** `.claude/specs/<slug>/` on
  `wf/<slug>` before code.
- **Capability map:** `.claude/notes/metalwalls-capability-map.md` (E1–E12) —
  path referenced by loop; recreate if missing on worker trees.
- **Orphans (not incomplete for this loop):**
  - `regressions/md/{dpd_ab_mix,dpd_groot_warren,fene_wca_chain}/`,
    `regressions/forcefield/yukawa_dimer/` — method-ladder shipped reg bodies,
    no `REFERENCE.md` (hygiene; not on MW queue; do not invent `[~]` rows)
  - local `wf/*` tips for wave-1 slugs fully merged into `dev` (stale branch names OK)

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin (or literature when no binary)
- `ctest` never launches LAMMPS or MetalWalls
- Wave-2 oracle tags: CG/predictor/finite-field/PIM/piston/TF/capacitance →
  **metalwalls**; Fumi–Tosi / Steele / Born → **lammps** (`born/coul/*`,
  `wall/lj93`, `pair born`)

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## Next after wave 2 (pending only)

1. `dipole-electrode-joint` (E5, oracle=metalwalls; **dep:** PIM)
2. `potential-profile-observer` (E10, oracle=metalwalls / literature Poisson)
3. `constant-force-pull-ion` (shared classical; oracle=lammps)
4. `surface-tension-coupling` (shared classical; oracle=lammps)

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only
under `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/grok-worktrees/`
(never bulk under `/home/jicli594`). Default **staff_count=10**. Spec-first:
write `.claude/specs/<slug>/` **before** code on each worker branch.
