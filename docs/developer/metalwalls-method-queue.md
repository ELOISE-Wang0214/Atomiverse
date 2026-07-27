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

Census: **2026-07-27 gap-census (wave-5 refresh)** (parent `dev`; waves 1–4
done=36; `## in_progress` empty).
**Tree audit (this census):** wave-4 land verified on `dev` (`96e16823`) — all
10 Components/regressions present: `ConstPotMaZe`,
`ConstPotBO::write_electrode_forces` / `::global_electroneutrality`,
`PairGaussianCoul` + `PairBorn` NeighborList bulk paths, `FixDfield`,
`FixMove`, `ChargeChargeAutocorr`, `ElastanceObserver`,
`reg_const-pot-tip4p-edl`. Map §1/§3/§8/§9 marks flipped accordingly
(2026-07-27 re-audit). Residual pending re-verified open against tree: no QEq,
no Voronoi, no XL dipoles, per-electrode \(\eta\) still a single shared
`ConstPotBO::eta`, no halo/MPI handling in `ConstPot*`, no PIM-NPT / IL-EDL /
nanoporous / TI regressions. Open MD/electrochem specs: **none**. Unmerged
in-scope `wf/<slug>` ahead_of_dev: **none**. **Preserved:** in_progress (0);
blocked (4); done (36).
**Foothold ≠ closed (drives new rows):** `ConstPotMaZe` is a CG+ASPC adiabatic
foothold — Coretti symplectic mass-zero constrained dynamics deferred
in-header; `FixDfield` is a vacuum \(E\equiv D\) foothold — no self-consistent
\(E = D - 4\pi P\) update; MW molten-salt+Al+MaZe example scenario now
reachable (MaZe + Fumi–Tosi NL footholds) but unbuilt.
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), waves 1–4 done (36), SCF/post-HF, GPW/PW, CUDA
  product, AIM/DAIM (E11 blocked), ML electrodes (E12 blocked), hybrid DFT
  electrodes, **full 2D-PBC product stack** (slab Yeh–Berkowitz ✓ is the
  foothold; true 2D Ewald out of scope), **PLUMED** bridge, cross-language
  product ports, non-electrochem MD residual (RESPA, FSSH, PIMD variants, …),
  XFT / linear-molecule / multipole-FF (P3 or off-spine)
- Handoff next-list items `coul-long-pme-nacl-lammps` / `nose-hoover-nvt-lammps` /
  `special-bonds-amber-lammps` / `lj-cut-lammps-dimer` **not pending** (already
  covered by excluded twins)

## pending

### E9 / scenarios — production EDL (beyond wave-4 aqueous staff)

- [ ] `nanoporous-carbon-const-pot` | Nanoporous carbon geometry + const-pot (product) | electrochem | slit / pore carbon + electrolyte | — | \(Q\), occupancy proxy | **oracle=metalwalls**; map §9 nanoporous ○ (after E1 production depth)
- [ ] `maze-molten-salt-al` | Molten-salt + Al wall + MaZe MW example scenario | electrochem | Fumi–Tosi NL salt + Al const-pot walls | — | \(Q\), \(g(r)\) proxy, short NVT | **oracle=metalwalls**; map §9 molten salt ~ (MaZe + Fumi–Tosi NL footholds landed wave-4; product box unbuilt — new 2026-07-27)

### E10 — TI residual (after charge-tool observers)

- [ ] `thermo-integration-capacitance` | Soft-core / λ TI workflow for electrode capacitance | electrochem | parallel-plate | — | \(\Delta A\) / \(C\) proxy | **oracle=literature**; map E10 TI workflow; soft-core pairs already ✓

## in_progress

*(empty — wave 5 merged and harnessed 2026-07-27)*

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

### Wave 4 — 2026-07-24 (manual land → `dev` tree)

| slug | Regression | oracle | ratio |
| ---- | ---------- | ------ | ----- |
| [x] `maze-charge-dynamics` | `regressions/md/maze-charge-dynamics` | metalwalls (CG+ASPC foothold) | perf_exempt |
| [x] `electrode-atom-forces` | `regressions/md/electrode-atom-forces` | literature (HF on elyt) | perf_exempt |
| [x] `global-electroneutrality-bo` | `regressions/md/global-electroneutrality-bo` | metalwalls/Scalfi | perf_exempt |
| [x] `gaussian-coul-neighborlist` | `regressions/md/gaussian-coul-neighborlist` | literature | perf_exempt |
| [x] `displacement-dfield-md` | `regressions/md/displacement-dfield-md` | metalwalls D-field foothold | perf_exempt |
| [x] `born-buck-neighborlist-bulk` | `regressions/md/born-buck-neighborlist-bulk` | lammps born | perf_exempt |
| [x] `const-pot-tip4p-edl` | `regressions/md/const-pot-tip4p-edl` | metalwalls | perf_exempt |
| [x] `charge-charge-autocorr` | `regressions/md/charge-charge-autocorr` | literature | perf_exempt |
| [x] `elastance-condition-observer` | `regressions/md/elastance-condition-observer` | literature | perf_exempt |
| [x] `constant-velocity-smd` | `regressions/md/constant-velocity-smd` | lammps fix move | perf_exempt |

### Wave 5 — 2026-07-27 (merged → `dev`)

| slug | Regression | oracle | ratio |
| ---- | ---------- | ------ | ----- |
| [x] `per-electrode-eta` | `regressions/md/per-electrode-eta` | literature (Reed 2007 DOI 10.1063/1.2464084; Scalfi 2020 DOI 10.1063/5.0028232; MW JOSS 2020) | perf_exempt |
| [x] `const-charge-electrode` | `regressions/md/const-charge-electrode` | lammps (ELECTRODE `conq` twin) | perf_exempt |
| [x] `const-pot-charge-halo-mpi` | `regressions/md/const-pot-charge-halo-mpi` | literature / self-consistency | perf_exempt |
| [x] `maze-symplectic-charge-dynamics` | `regressions/md/maze-symplectic-charge-dynamics` | metalwalls (Coretti 2020 literature pin) | perf_exempt |
| [x] `constant-d-polarization-feedback` | `regressions/md/constant-d-polarization-feedback` | metalwalls (Zhang–Sprik constant-D) | perf_exempt |
| [x] `extended-lagrangian-dipoles` | `regressions/md/extended-lagrangian-dipoles` | metalwalls (iEL/SCF literature pins — Albaugh/Niklasson; MW MaZe route cited) | perf_exempt (no MW binary; no LAMMPS iEL/SCF twin) |
| [x] `qeq-fluctuating-charge` | `regressions/md/qeq-fluctuating-charge` | lammps (`fix qeq/point`) | **0.45** |
| [x] `pim-npt-polarizable-il` | `regressions/md/pim-npt-polarizable-il` | metalwalls (offline numpy Thole-PIM pin; MW JOSS 2020, Thole 1981) | perf_exempt |
| [x] `voronoi-volume-tf` | `regressions/md/voronoi-volume-tf` | metalwalls | perf_exempt |
| [x] `const-pot-il-edl-box` | `regressions/md/const-pot-il-edl-box` | metalwalls | perf_exempt |
