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

Census: **2026-07-24 wave 4 gap-census** (parent `dev` @ post–harness wave 3).
Recover re-check: **incomplete_count=0** (A–E empty). Open MD/electrochem
specs: **none**. Unmerged in-scope `wf/<slug>` with `ahead_of_dev>0`: **none**.
Blocked (4) preserved. **Foothold ≠ closed:** waves 1–3 shipped 26 textbook /
observer / Component footholds; production-depth and residual map ○ rows reopened
below (NeighborList salt/Gaussian pairs, MaZe, D-field, electrode BO forces,
charge tools, CV-SMD, aqueous const-pot product, QEq / XL-PIM / TI).
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), wave-1 done (10), wave-2 done (10), wave-3 done (6),
  SCF/post-HF, GPW/PW, CUDA product, AIM/DAIM (E11 blocked), ML electrodes
  (E12 blocked), hybrid DFT electrodes, **full 2D-PBC product stack** (slab
  Yeh–Berkowitz ✓ is the foothold; true 2D Ewald out of scope), **PLUMED**
  bridge, cross-language product ports, non-electrochem MD residual (RESPA,
  FSSH, PIMD variants, …)
- Handoff next-list items `coul-long-pme-nacl-lammps` / `nose-hoover-nvt-lammps` /
  `special-bonds-amber-lammps` / `lj-cut-lammps-dimer` **not pending** (already
  covered by excluded twins)

## pending

### E4 / E5 — polarization depth (staffed later if wave-4 capacity remains)

- [ ] `extended-lagrangian-dipoles` | Extended-Lagrangian / mass-zero dipole dynamics | polarizable | bulk polarizable toy | — | \(\mu(t)\), energy | **oracle=metalwalls**; map §2 XL dipoles ○ (SCF PIM ✓ is foothold)
- [ ] `qeq-fluctuating-charge` | QEq / fluctuating-charge (≠ const-pot electrode) | polarizable | small organic / ionic toy | — | \(q_i\), \(E\) | **oracle=lammps** `fix qeq` / charge equilibration; map §2 FQ/QEq ○

### E10 — TI residual (after charge-tool observers)

- [ ] `thermo-integration-capacitance` | Soft-core / λ TI workflow for electrode capacitance | electrochem | parallel-plate | — | \(\Delta A\) / \(C\) proxy | **oracle=literature**; map E10 TI workflow; soft-core pairs already ✓

## in_progress

**Wave-4 assign (2026-07-24)** — production depth first (NeighborList salt/Gaussian,
MaZe, electrode BO forces, CV-SMD, aqueous const-pot+TIP4P, D-field, charge tools).
Staff_count=10; incomplete was 0. Spec-first on `wf/<slug>` before impl.

### E1 / E2 — const-pot depth

- [~] `maze-charge-dynamics` | MaZe mass-zero electrode charge dynamics | electrochem | parallel-plate / molten-salt toy | — | \(Q(t)\), energy drift proxy | **oracle=metalwalls**; map §1 MaZe ★ after BO/CG ✓
- [~] `electrode-atom-forces` | Forces on electrode atoms from BO \(Q\) (term breakdown) | electrochem | movable wall toy | — | \(F\) on electrode sites | **oracle=metalwalls**; map §1 P1; needed for true movable const-pot walls (piston uses fluid \(F_z\) only)
- [~] `global-electroneutrality-bo` | Global electroneutrality \(\sum Q_i=-\sum q_{\mathrm{elec}}\) on `ConstPotBO` path | electrochem | charged electrolyte + walls | — | neutrality residual | **oracle=metalwalls**; map §1 ~ unify with legacy `ConstPotentialBO`
- [~] `gaussian-coul-neighborlist` | NeighborList-backed `PairGaussianCoul` (production depth) | electrochem | multi-site Gaussian wall | — | \(E\), \(F\) vs O(N²) twin | **oracle=literature** Reed/MW Gaussian; map §1 Gaussian ~ O(N²) foothold only

### E3 — external / finite field residual

- [~] `displacement-dfield-md` | Displacement \(\mathbf{D}\)-field finite-field MD | electrochem | planar capacitor | — | \(\langle P\rangle\), \(Q\) response | **oracle=metalwalls**; map §1 P2 Displacement D-field ○

### E8 — salt / wall production depth

- [~] `born-buck-neighborlist-bulk` | NeighborList `PairBorn`/`PairBuck` + bulk Fumi–Tosi molten-salt pin | forcefield | NaCl melt box | — | \(E\), \(P\), short NVE drift | **oracle=lammps** `born`/`buck` + bulk; map §3 Fumi–Tosi/Born ~; FF roadmap note O(N²)

### E9 — aqueous electrodes product

- [~] `const-pot-tip4p-edl` | Const-pot parallel-plate + TIP4P / virtual-site water EDL | electrochem | planar walls + TIP4P | — | \(Q\), \(\Psi(z)\), short NVT | **oracle=metalwalls**; map §9 TIP4P+const-pot ○; deps: E1 ✓ + vsites ✓

### E10 — observe / analysis residual

- [~] `charge-charge-autocorr` | Electrode charge–charge autocorrelation \(C_{QQ}(t)\) | electrochem | parallel-plate EDL | — | \(C_{QQ}(t)\) / \(\tau\) | **oracle=literature** / metalwalls tools; map §1/§8 P2
- [~] `elastance-condition-observer` | Empty capacitance / elastance \(\mathrm{cond}(A)\) debug observer | electrochem | empty parallel plate | — | \(\mathrm{cond}(A)\), empty \(C\) | **oracle=literature**; map §8 P1 empty capacitance ○

### Shared classical MD ○ (LAMMPS primary; electrochem-adjacent)

- [~] `constant-velocity-smd` | Constant-velocity pulling / SMD residual | md | ion / tagged atom | — | \(z(t)\), work proxy | **oracle=lammps** `fix smd` / `fix move`; MD roadmap ○ CV pull; after constant-force ✓

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
