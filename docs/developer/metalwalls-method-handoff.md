# MetalWalls method loop — handoff

Last harness: **wave 3 — 2026-07-24** on `dev`
(`chore(harness): metalwalls wf wave 3`). Prior: harness wave 2
(`chore(harness): metalwalls wf wave 2`); assign wave 3
(`chore(queue): metalwalls assign wave 3`).

**Workflow order (binding):** Recover incomplete first → Census →
Assign (`incomplete ≫ pending`) → Worktrees (**spec-first**) → Merge → Harness →
Tutorial. See `docs/developer/metalwalls-method-loop.md`.

## Recover inventory (post harness wave 3, 2026-07-24)

| Source | Result |
| ------ | ------ |
| A queue `## in_progress` `[~]` | **0** |
| B open `.claude/specs/` MD/electrochem | **0** (wave-3 specs deleted; PySCF deferred specs ignored) |
| C unmerged `wf/<slug>` (this loop) | **0** (wave-1/2/3 tips ancestor of `dev`) |
| D handoff incomplete rows | **none** |
| E half-written `regressions/{md,forcefield}/` | **0** |

**Resume set:** empty. **incomplete_remaining=0**. Classical MW staff queue
exhausted (pending empty; blocked 4 remain out of scope).

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

## Wave 3 summary (harnessed 2026-07-24)

| Slug | Domain | Status | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ------ | ------ | ---------- | ------ | ---------------------------------------------- |
| `constant-force-pull-ion` | md (shared) | **merged → done** | `regressions/md/constant-force-pull-ion` | lammps | **0.12** |
| `oscillating-efield-md` | electrochem E3 | **merged → done** | `regressions/md/oscillating-efield-md` | lammps | **2.12** |
| `potential-profile-observer` | electrochem E10 | **merged → done** | `regressions/md/potential-profile-observer` | literature | **perf_exempt** |
| `surface-tension-coupling` | md (shared) | **merged → done** | `regressions/md/surface-tension-coupling` | lammps | **2.88** |
| `virtual-site-multi-topology` | md E9 | **merged → done** | `regressions/md/virtual-site-multi-topology` | lammps | **0.084** |
| `dipole-electrode-joint` | electrochem E5 | **merged → done** | `regressions/md/dipole-electrode-joint` | metalwalls | **perf_exempt** |

### Incomplete remaining (wave 3)

*(none — all six acceptance-met and merged)*

- **Blocked this wave:** none of the staffed six; global blocked set still
  `cuda-const-pot`, `aim-daim-product`, `ml-pinnwall-electrode`,
  `periodic-hybrid-electrode-dft` (out of CPU MD loop)
- **Merge-failed this wave:** none
- **Incomplete remaining after harness wave 3:** **0**
- **Queue after harness wave 3:** pending = **0**; in_progress = **0**;
  done = **26** (wave 1 + 2 + 3); blocked = **4**
- **Roadmap (CPU ○→✓ on `dev` after wave 3):** constant-surface-tension
  coupling / NPγT (`SurfaceTensionNPT`, couple xy); oscillating E-field
  (`FixOscillatingEfield`); constant-force pull (`FixAddForce`); multi-topology
  virtual sites (`VirtualSites` Weighted3 + OutOfPlane3 TIP4P/TIP5P);
  dipole⊕electrode joint (`DipoleElectrodeJoint`); slab \(\Psi(z)\) observer
  (`PotentialProfile` / planar Poisson). All Components on `cpu::md::MD` /
  Evaluate — not new Driver categories; charge solve not inside `KSpace*`.
- **Specs closed:** all 6 wave-3 slugs (acceptance met on `dev`; deleted under
  `.claude/specs/`).
- **Capability map:** no lasting rule change this harness (oracle priority and
  Component-on-MD architecture unchanged).
- **Orphans (not incomplete for this loop):**
  - `regressions/md/{dpd_ab_mix,dpd_groot_warren,fene_wca_chain}/`,
    `regressions/forcefield/yukawa_dimer/` — method-ladder shipped reg bodies,
    no `REFERENCE.md` (hygiene; not on MW queue; do not invent `[~]` rows)
  - local `wf/*` tips for wave-1/2/3 slugs fully merged into `dev` (stale branch
    names OK)

## Oracle

- Overlap with LAMMPS → **LAMMPS** offline pin (`tools/metalwalls_oracle/` + LAMMPS inputs)
- MW-only → **MetalWalls** offline pin (or literature when no binary)
- `ctest` never launches LAMMPS or MetalWalls
- Wave-3 oracle tags: force-pull / oscillating efield / surface tension /
  virtual sites → **lammps**; \(\Psi(z)\) → **literature**; dipole⊕electrode
  joint → **metalwalls**

## Build (workers)

MD-only flags: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
`ATV_BUILD_REGRESSIONS=ON`, **HF/PW/GPW/CUDA/Torch OFF**.

## Next (after wave 3 — queue empty)

Classical MetalWalls staff queue is **exhausted** on `dev` (26 done, 0 pending,
0 in_progress). Next shift options:

1. **Census** from map E11/E12 or residual MD ○ (constant-velocity SMD, RESPA,
   …) only if product-scheduled — do not invent rows casually.
2. **Tutorial** wave for wave-3 methods if docs lag code.
3. Unblock tracks remain out of loop: CUDA const-pot, AIM/DAIM, ML PINN wall,
   periodic hybrid DFT electrodes.

Parent stays on **`dev`**. Workers use isolation worktrees `wf/<slug>` only
under `/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/grok-worktrees/`
(never bulk under `/home/jicli594`). Default **staff_count=10**. Spec-first:
write `.claude/specs/<slug>/` **before** code on each worker branch.
