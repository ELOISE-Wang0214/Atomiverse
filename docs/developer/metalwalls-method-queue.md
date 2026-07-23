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

Census: **2026-07-23** (parent `dev`). **Harness wave 1 2026-07-23:** all **10**
staff-batch `[~]` merged → `## done` (regressions + pins on `dev`). **Assign
wave 2 2026-07-24:** **10** unblocked pending → `## in_progress` `[~]` (staff=10).
Skipped `dipole-electrode-joint` (dep **PIM** still open). Open MD/electrochem
specs: **none** (workers re-create on `wf/<slug>`). Unmerged in-scope wave-1
`wf/<slug>` tips: HEADs already merged into `dev` (not incomplete).
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), SCF/post-HF, GPW/PW, CUDA product, AIM/DAIM (E11 blocked),
  ML electrodes (E12 blocked), hybrid DFT electrodes
- Handoff next-list items `coul-long-pme-nacl-lammps` / `nose-hoover-nvt-lammps` /
  `special-bonds-amber-lammps` / `lj-cut-lammps-dimer` **not assigned** (already
  covered by excluded twins)

## pending

### E1 — constant-potential v1 (map §1 P0)

*(wave-1 E1 core → done 2026-07-23)*

### E2 — electrostatic rigor (after E1 foothold)

*(wave-2 E2 CG + predictor → in_progress)*

### E3 — external field (map §1 P1; MD roadmap ✓ constant E-field foothold)

*(wave-2 finite-field single electrode → in_progress)*

### E4 / E5 — polarization (map §2; FF roadmap ○ induced dipoles / Thole)

- [ ] `dipole-electrode-joint` | Coupled dipole ⊕ electrode \(Q\) solve | electrochem | IL + graphene walls | — | \(Q\), \(\mu\) | **oracle=metalwalls**; map E5; **dep:** E1 ✓ + PIM (**hold** until `pim-induced-dipole-bulk` lands)

### E6 / E7 — electrode mechanics & metallicity

*(wave-2 piston + Thomas–Fermi → in_progress)*

### E8 — salt / wall potentials (map §3; FF ✓ Buckingham foothold)

*(wave-2 Fumi–Tosi + Steele + Born → in_progress)*

### E9 — molecular water at electrodes (MD roadmap ✓ TIP4P M-site foothold)

*(four-site-tip4p-foothold → done wave 1; general multi-topology virtual sites still open)*

### E10 — ecosystem (after E1 observers)

- [ ] `potential-profile-observer` | Electrostatic potential profile (slab) | electrochem | parallel plate | — | \(\Psi(z)\) curve | **oracle=metalwalls** / literature Poisson; map E10; **dep:** charges + long-range ✓

### Shared classical MD ○ still open (LAMMPS primary; electrochem-adjacent)

- [ ] `constant-force-pull-ion` | Constant-force / steered ion pull (noneq MD ○) | md | ion between walls | — | \(\langle z\rangle\), work proxy | **oracle=lammps** `fix smd` / `fix addforce`; useful EDL transport twin — **not** MW-only
- [ ] `surface-tension-coupling` | Constant surface-tension ensemble ○ | md | LJ slab | — | \(\gamma\) / \(L_z\) | **oracle=lammps** `fix npt` surface tension path; optional for planar interfaces

## in_progress

### Wave 2 staff batch — assigned 2026-07-24 (parent `dev`)

- [~] `const-pot-cg-large-electrode` | CG electrode charge solve | electrochem | larger wall / many electrode sites | — | \(Q\) match dense path | **oracle=metalwalls**; dep: BO dense inverse path ✓
- [~] `electrode-charge-predictor` | Electrode charge predictor (ASPC-like) | electrochem | short MD | — | fewer CG/BO iters | **oracle=metalwalls**; map E2
- [~] `finite-field-single-electrode` | Finite-field + single electrode (3D) | electrochem | single wall + electrolyte | — | \(Q\), polarization | **oracle=metalwalls** Dufils-style; **dep:** E1 ✓ + efield ✓
- [~] `pim-induced-dipole-bulk` | PIM self-consistent induced dipoles (bulk) | polarizable | molten salt / IL toy | — | \(E_{\mathrm{pol}}\), forces | **oracle=metalwalls**; map E4; Thole/TT damping in-scope
- [~] `piston-electrode-pz` | Piston electrode constant \(P_z\) | electrochem | IL + graphene / parallel plate | — | \(\langle P_z\rangle\) | **oracle=metalwalls**; map E6; **dep:** E1 ✓ + electrode forces
- [~] `thomas-fermi-electrode` | Thomas–Fermi finite metallicity | electrochem | gold capacitor toy | — | \(C\), \(Q\) response | **oracle=metalwalls**; map E7
- [~] `fumi-tosi-nacl` | Fumi–Tosi molten salt | forcefield | NaCl melt | — | \(E\) | **oracle=lammps** (`pair born/coul/*`); map E8; **dep:** Buckingham pair ✓
- [~] `steele-wall-argon` | Steele wall potential | forcefield | Ar / graphite wall | — | \(E\), \(F\) | **oracle=lammps** (`wall/lj93` Steele 10-4-3); map E8
- [~] `born-repulsion-foothold` | Born repulsion wall/ion term | forcefield | ion near wall | — | \(E\), \(F\) | **oracle=lammps** (`pair born`); map E8 companion
- [~] `capacitance-charge-fluctuation` | Capacitance / charge-fluctuation analysis | electrochem | parallel-plate EDL | — | \(C\) estimate | **oracle=metalwalls** tools; map E10; **dep:** E1 ✓ + charge observer ✓

## blocked

- [!] `cuda-const-pot` | blocked — CUDA electrode track is separate; this loop is CPU MD-only
- [!] `aim-daim-product` | blocked until map E11 is explicitly scheduled; not E1
- [!] `ml-pinnwall-electrode` | blocked — map E12 research; external ML interface
- [!] `periodic-hybrid-electrode-dft` | blocked — quantum electrodes → CP2K/PySCF loops, not classical MW loop

## done

### Wave 1 — 2026-07-23 (merged → `dev`)

| Slug | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
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
