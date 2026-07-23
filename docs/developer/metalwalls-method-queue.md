# MetalWalls capability queue

Durable worklist for workflow `metalwalls-method-loop`. Do not delete done rows.

**Recover first:** every shift finishes `## in_progress` `[~]`, open
`.claude/specs/` (MD/electrochemistry scope), and unmerged `wf/<slug>` before
new pending. Census only rewrites **pending**. Assign: **incomplete ≫ pending**.

**Oracle rule:** if LAMMPS implements the same algorithm/feature, pin and fair
bench against **LAMMPS**. MetalWalls only for MW-only rows (mark
`oracle=metalwalls` in the row rationale).

**Build rule:** workers use **MD-only** cmake (HF/PW/GPW/CUDA/Torch **OFF**,
MD+KSPACE+TESTS+REGRESSIONS **ON**). See `metalwalls-method-loop.md`.

Census source: `.claude/notes/metalwalls-capability-map.md` slices **E1–E12**,
plus classical MD/FF ○ rows that LAMMPS or MetalWalls can oracle offline.
(Map path referenced by loop; **missing on this tree** — pending rebuilt from
prior census rows, handoff, `tools/metalwalls_oracle/README.md`, and
`docs/roadmap/{molecular-dynamics,force-fields}.md` ○.)

Census: **2026-07-24** (parent `dev`, post–harness wave 2). Assign **wave 3**
promoted the last **6** pending rows → `## in_progress` `[~]` (staff default 10;
only 6 remained). **incomplete_count=6**. Open MD/electrochem specs: **none**
yet (workers write `.claude/specs/<slug>/` on `wf/<slug>` before code). Unmerged
in-scope `wf/<slug>` with `ahead_of_dev>0` for wave-3 slugs: **none** at assign
time. Blocked (4) skipped. No brand-new census rows invented.
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), wave-1 done (10), wave-2 done (10), wave-3 in_progress (6),
  SCF/post-HF, GPW/PW, CUDA product, AIM/DAIM (E11 blocked), ML electrodes
  (E12 blocked), hybrid DFT electrodes, non-electrochem MD ○ (RESPA, FSSH,
  PIMD variants, …)
- Handoff next-list items `coul-long-pme-nacl-lammps` / `nose-hoover-nvt-lammps` /
  `special-bonds-amber-lammps` / `lj-cut-lammps-dimer` **not pending** (already
  covered by excluded twins)

## pending

*(empty after assign wave 3 — 2026-07-24; all remaining open rows staffed)*

### E1 — constant-potential v1 (map §1 P0)

*(wave-1 E1 core → done 2026-07-23)*

### E2 — electrostatic rigor (after E1 foothold)

*(wave-2 CG + predictor → done 2026-07-24)*

### E3 — external field (map §1 P1; MD roadmap ✓ constant E-field foothold)

*(wave-2 finite-field single electrode → done 2026-07-24; wave-3 oscillating residual → in_progress)*

### E4 / E5 — polarization (map §2; FF roadmap ✓ bulk PIM + Thole foothold)

*(wave-2 PIM bulk → done 2026-07-24; wave-3 dipole⊕electrode joint → in_progress)*

### E6 / E7 — electrode mechanics & metallicity

*(wave-2 piston + Thomas–Fermi → done 2026-07-24)*

### E8 — salt / wall potentials (map §3; FF ✓ Buckingham + Born BMH + Fumi–Tosi + Steele footholds)

*(wave-2 Fumi–Tosi + Steele + Born → done 2026-07-24)*

### E9 — molecular water at electrodes (MD roadmap ✓ TIP4P M-site foothold)

*(four-site-tip4p-foothold → done wave 1; wave-3 multi-topology virtual sites → in_progress)*

### E10 — ecosystem (after E1 observers)

*(wave-2 capacitance / charge-fluctuation → done 2026-07-24; wave-3 potential profile → in_progress)*

### Shared classical MD ○ still open (LAMMPS primary; electrochem-adjacent)

*(wave-3 constant-force pull + surface-tension coupling → in_progress)*

## in_progress

### Wave 3 — assign 2026-07-24 (staff=6; all prior pending; is_resume=false)

- [~] `dipole-electrode-joint` | Coupled dipole ⊕ electrode \(Q\) solve | electrochem | IL + graphene walls | — | \(Q\), \(\mu\) | **oracle=metalwalls**; map E5; **dep met:** E1 ✓ + PIM ✓ (`pim-induced-dipole-bulk`)
- [~] `potential-profile-observer` | Electrostatic potential profile (slab) | electrochem | parallel plate | — | \(\Psi(z)\) curve | **oracle=metalwalls** / literature Poisson; map E10; charges + long-range ✓
- [~] `oscillating-efield-md` | Oscillating / AC electric field MD | electrochem | ions / LJ+Coul fluid | — | \(F(t)=qE(t)\), \(\Delta E\) | **oracle=lammps** `fix efield` variable / equal-style; map E3 residual
- [~] `virtual-site-multi-topology` | General multi-topology virtual sites (beyond TIP4P M-site pair) | md | TIP4P/TIP5P textbook box | — | \(E\), geometry / force partition | **oracle=lammps** tip4p / virtual-site path; map E9 residual
- [~] `constant-force-pull-ion` | Constant-force / steered ion pull (noneq MD ○) | md | ion between walls | — | \(\langle z\rangle\), work proxy | **oracle=lammps** `fix smd` / `fix addforce`; EDL transport twin
- [~] `surface-tension-coupling` | Constant surface-tension ensemble ○ | md | LJ slab | — | \(\gamma\) / \(L_z\) | **oracle=lammps** `fix npt` surface tension path

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

### Wave 2 — 2026-07-24 (merged → `dev`)

| Slug | Regression | oracle | ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) |
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
