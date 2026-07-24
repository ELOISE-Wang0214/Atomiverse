# MetalWalls capability queue

Durable worklist for workflow `metalwalls-method-loop`. Do not delete done rows.

**Recover first:** every shift finishes `## in_progress` `[~]`, open
`.claude/specs/` (MD/electrochemistry scope), and unmerged `wf/<slug>` before
new pending. Census only rewrites **pending**. Assign: **incomplete ≫ pending**.

**Oracle rule:** if LAMMPS implements the same algorithm/feature, pin and fair
10→bench against **LAMMPS**. MetalWalls only for MW-only rows (mark
`oracle=metalwalls` in the row rationale).

**Build rule:** workers use **MD-only** cmake (HF/PW/GPW/CUDA/Torch **OFF**,
MD+KSPACE+TESTS+REGRESSIONS **ON**). See `metalwalls-method-loop.md`.

Census source: `.claude/notes/metalwalls-capability-map.md` slices **E1–E12**,
plus classical MD/FF ○ rows that LAMMPS or MetalWalls can oracle offline.

Census: **2026-07-24** (parent `dev`, post–harness wave 3). Recover re-check
same day: **incomplete_count=0** (A–E empty). All open staff rows from assign
wave 3 are **merged → done**. Open MD/electrochem specs: **none** (wave-3
specs deleted). Unmerged in-scope `wf/<slug>` with `ahead_of_dev>0`: **none**
(tips ancestor of `dev`). Blocked (4) skipped. **pending empty** — no further
classical MW queue rows without a new census inventing work.
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), wave-1 done (10), wave-2 done (10), wave-3 done (6),
  SCF/post-HF, GPW/PW, CUDA product, AIM/DAIM (E11 blocked), ML electrodes
  (E12 blocked), hybrid DFT electrodes, non-electrochem MD ○ (RESPA, FSSH,
  PIMD variants, constant-velocity SMD residual, …)
- Handoff next-list items `coul-long-pme-nacl-lammps` / `nose-hoover-nvt-lammps` /
  `special-bonds-amber-lammps` / `lj-cut-lammps-dimer` **not pending** (already
  covered by excluded twins)

## pending

*(empty after harness wave 3 — 2026-07-24; staffed MW classical queue exhausted)*

### E1 — constant-potential v1 (map §1 P0)

*(wave-1 E1 core → done 2026-07-23)*

### E2 — electrostatic rigor (after E1 foothold)

*(wave-2 CG + predictor → done 2026-07-24)*

### E3 — external field (map §1 P1; MD roadmap ✓ constant + oscillating E-field)

*(wave-2 finite-field single electrode → done 2026-07-24; wave-3 oscillating → done 2026-07-24)*

### E4 / E5 — polarization (map §2; FF roadmap ✓ bulk PIM + Thole + dipole⊕electrode)

*(wave-2 PIM bulk → done 2026-07-24; wave-3 dipole⊕electrode joint → done 2026-07-24)*

### E6 / E7 — electrode mechanics & metallicity

*(wave-2 piston + Thomas–Fermi → done 2026-07-24)*

### E8 — salt / wall potentials (map §3; FF ✓ Buckingham + Born BMH + Fumi–Tosi + Steele footholds)

*(wave-2 Fumi–Tosi + Steele + Born → done 2026-07-24)*

### E9 — molecular water at electrodes (MD roadmap ✓ TIP4P + multi-topology virtual sites)

*(four-site-tip4p-foothold → done wave 1; virtual-site-multi-topology → done wave 3)*

### E10 — ecosystem (after E1 observers)

*(wave-2 capacitance / charge-fluctuation → done 2026-07-24; wave-3 potential profile → done 2026-07-24)*

### Shared classical MD ○ still open (LAMMPS primary; electrochem-adjacent)

*(wave-3 constant-force pull + surface-tension coupling → done 2026-07-24; residual ○ outside this loop: constant-velocity SMD, RESPA, …)*

## in_progress

*(empty after harness wave 3 — 2026-07-24)*

## blocked

- [!] `cuda-const-pot` | blocked — CUDA electrode track is separate; this loop is CPU MD-only
- [!] `aim-daim-product` | blocked until map E11 is explicitly scheduled; not E1
- [!] `ml-pinnwall-electrode` | blocked — map E12 research; external ML interface
- [!] `periodic-hybrid-electrode-dft` | blocked — quantum electrodes → CP2K/PySCF loops, not classical MW loop

## done

### Wave 1 — 2026-07-23 (merged → `dev`)

| slug | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ---------- | ------ | ---------------------------------------------- |
| [x] `const-pot-bo-parallel-plate` | `regressions/md/const-pot-bo-parallel-plate` | lammps (ELECTRODE madelung) | **perf_exempt** (no ELECTRODE binary) |
| [x] `freeze-electrode-atoms` | `regressions/md/freeze-electrode-atoms` | lammps | **3.98** |
| [x] `const-pot-charge-observer` | `regressions/md/const-pot-charge-observer` | lammps | **perf_exempt** (micro-sum vs startup) |
| [x] `charging-work-energy` | `regressions/md/charging-work-energy` | literature (MW JOSS / Scalfi) | **perf_exempt** |
| [x] `const-pot-md-nvt-short` | `regressions/md/const-pot-md-nvt-short` | metalwalls | **perf_exempt** (no MW/ELECTRODE twin) |
| [x] `gaussian-electrode-charge` | `regressions/md/gaussian-electrode-charge` | literature (Reed / MW Gaussian) | **perf_exempt** |
| [x] `slab-2d-ewald-foothold` | `regressions/md/slab-2d-ewald-foothold` | lammps | **6.49** |
| [x] `external-efield-md` | `regressions/md/external-efield-md` | lammps | **0.60** |
| [x] `buckingham-pair-lammps` | `regressions/forcefield/buckingham-pair-lammps` | lammps | **perf_exempt** (dimer startup) |
| [x] `four-site-tip4p-foothold` | `regressions/md/four-site-tip4p-foothold` | lammps | **0.338** |

### Wave 2 — 2026-07-24 (merged → `dev`)

| slug | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ---------- | ------ | ---------------------------------------------- |
| [x] `const-pot-cg-large-electrode` | `regressions/md/const-pot-cg-large-electrode` | metalwalls (dual: dense BO) | **perf_exempt** (no MW binary) |
| [x] `electrode-charge-predictor` | `regressions/md/electrode-charge-predictor` | metalwalls (ASPC) | **perf_exempt** |
| [x] `finite-field-single-electrode` | `regressions/md/finite-field-single-electrode` | metalwalls / literature | **perf_exempt** |
| [x] `pim-induced-dipole-bulk` | `regressions/md/pim-induced-dipole-bulk` | metalwalls | **perf_exempt** |
| [x] `piston-electrode-pz` | `regressions/md/piston-electrode-pz` | metalwalls | **perf_exempt** |
| [x] `thomas-fermi-electrode` | `regressions/md/thomas-fermi-electrode` | metalwalls | **perf_exempt** |
| [x] `fumi-tosi-nacl` | `regressions/forcefield/fumi-tosi-nacl` | lammps | **perf_exempt** (4-ion startup) |
| [x] `steele-wall-argon` | `regressions/forcefield/steele-wall-argon` | lammps (`wall/lj1043`) | **perf_exempt** |
| [x] `born-repulsion-foothold` | `regressions/forcefield/born-repulsion-foothold` | lammps | **perf_exempt** |
| [x] `capacitance-charge-fluctuation` | `regressions/md/capacitance-charge-fluctuation` | literature | **perf_exempt** |

### Wave 3 — 2026-07-24 (merged → `dev`)

| slug | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
| ---- | ---------- | ------ | ---------------------------------------------- |
| [x] `constant-force-pull-ion` | `regressions/md/constant-force-pull-ion` | lammps (`fix addforce`) | **0.12** |
| [x] `oscillating-efield-md` | `regressions/md/oscillating-efield-md` | lammps (`fix efield` variable) | **2.12** |
| [x] `potential-profile-observer` | `regressions/md/potential-profile-observer` | literature (planar Poisson) | **perf_exempt** |
| [x] `surface-tension-coupling` | `regressions/md/surface-tension-coupling` | lammps (`fix npt` couple xy) | **2.88** |
| [x] `virtual-site-multi-topology` | `regressions/md/virtual-site-multi-topology` | lammps (tip4p path) | **0.084** |
| [x] `dipole-electrode-joint` | `regressions/md/dipole-electrode-joint` | metalwalls | **perf_exempt** (no MW binary) |
