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

Census: **2026-07-23** (parent `dev`). **incomplete_count=0** for this loop.
- `## in_progress` `[~]`: empty
- live `.claude/specs/` (MD/electrochem/const-pot/forcefield): **none**
  (`INDEX.md` Live empty)
- unmerged `wf/<slug>` in this loop’s scope: **none** (local `wf/*` tips are
  PySCF/post-HF ancestors of `dev`, out of scope; no metalwalls/const-pot/
  electrode/efield branches local or on `origin`)
- handoff incomplete rows: **0**
- half-written `regressions/{md,forcefield,longrange}/` for *this* queue: **none**
- **Excluded from pending:** already ✓ classical twins with in-tree pins
  (`lj_argon`, `npt_pme_nacl` Madelung, `gaff_peo_fragment` specials, NH/NPT
  integrators), SCF/post-HF, GPW/PW, CUDA product, AIM/DAIM (E11 blocked),
  ML electrodes (E12 blocked), hybrid DFT electrodes

## pending

### E1 — constant-potential v1 (map §1 P0) — assign first

- [ ] `const-pot-bo-parallel-plate` | Constant-potential BO solve (parallel plate) | electrochem | NaCl slab + two planar walls | — | electrode \(Q\), \(\Delta\Psi\), neutrality | **oracle=metalwalls** (LAMMPS ELECTRODE only if package available offline); map E1; Component on `md::MD`; private \(A,Q\); write `kAtom_charge`; freeze walls
- [ ] `freeze-electrode-atoms` | Freeze / immobile electrode atoms | md | any textbook | — | \(x\) fixed, \(v=0\) | **oracle=lammps** `fix setforce 0` / `fix freeze`; Component on integrate path; E1 dep for wall atoms
- [ ] `const-pot-charge-observer` | Electrode total charge observer | electrochem | same foothold | — | \(Q_{\mathrm{tot}}(t)\) per electrode | **oracle=lammps** electrode charge dump if available; else **oracle=metalwalls**; dep: charge solve
- [ ] `charging-work-energy` | Charging work \(-\sum\Psi_i Q_i\) in energy | electrochem | parallel plate | — | energy channel / meta | **oracle=metalwalls** / Scalfi literature; dep: const-pot
- [ ] `const-pot-md-nvt-short` | Const-pot + MD NVT short traj | electrochem | NaCl + parallel walls | — | stable \(E\), \(Q\) sign, no NaN | **oracle=metalwalls** (or lammps electrode MD); dep: BO + freeze + CoulLong/PME

### E2 — electrostatic rigor (after E1 foothold)

- [ ] `gaussian-electrode-charge` | Gaussian charge kernels on electrode atoms | electrochem | parallel plate | — | \(E\), \(Q\) vs point-charge \(A\) | **oracle=metalwalls**; map E2; document point-charge approx if deferred
- [ ] `const-pot-cg-large-electrode` | CG electrode charge solve | electrochem | larger wall / many electrode sites | — | \(Q\) match dense path | **oracle=metalwalls**; dep: BO dense inverse path
- [ ] `electrode-charge-predictor` | Electrode charge predictor (ASPC-like) | electrochem | short MD | — | fewer CG/BO iters | **oracle=metalwalls**; map E2
- [ ] `slab-2d-ewald-foothold` | 2D / slab Ewald (or `kspace_modify slab`) | electrochem | slab NaCl / parallel plate | — | \(E\) vs 3D+slab correction | **oracle=lammps** `kspace_modify slab` when available; else **oracle=metalwalls** 2D; map E2

### E3 — external field (map §1 P1; MD roadmap ○ constant E-field)

- [ ] `external-efield-md` | Applied uniform E-field on MD charges | electrochem | ions / LJ+Coul fluid | — | \(F=qE\), \(\Delta E\) | **oracle=lammps** `fix efield` (LAMMPS primary); map E3
- [ ] `finite-field-single-electrode` | Finite-field + single electrode (3D) | electrochem | single wall + electrolyte | — | \(Q\), polarization | **oracle=metalwalls** Dufils-style; **dep:** E1 + efield

### E4 / E5 — polarization (map §2; FF roadmap ○ induced dipoles / Thole)

- [ ] `pim-induced-dipole-bulk` | PIM self-consistent induced dipoles (bulk) | polarizable | molten salt / IL toy | — | \(E_{\mathrm{pol}}\), forces | **oracle=metalwalls**; map E4; Thole/TT damping in-scope
- [ ] `dipole-electrode-joint` | Coupled dipole ⊕ electrode \(Q\) solve | electrochem | IL + graphene walls | — | \(Q\), \(\mu\) | **oracle=metalwalls**; map E5; **dep:** E1 + PIM

### E6 / E7 — electrode mechanics & metallicity

- [ ] `piston-electrode-pz` | Piston electrode constant \(P_z\) | electrochem | IL + graphene / parallel plate | — | \(\langle P_z\rangle\) | **oracle=metalwalls**; map E6; **dep:** E1 + electrode forces
- [ ] `thomas-fermi-electrode` | Thomas–Fermi finite metallicity | electrochem | gold capacitor toy | — | \(C\), \(Q\) response | **oracle=metalwalls**; map E7

### E8 — salt / wall potentials (map §3; FF ○ Buckingham-class)

- [ ] `buckingham-pair-lammps` | Buckingham pair (Born–Mayer foothold) | forcefield | ion dimer / bulk toy | — | \(E\), \(F\) | **oracle=lammps** `pair_style buck`; textbook precursor to Fumi–Tosi
- [ ] `fumi-tosi-nacl` | Fumi–Tosi molten salt | forcefield | NaCl melt | — | \(E\) | **oracle=metalwalls** / literature; map E8
- [ ] `steele-wall-argon` | Steele wall potential | forcefield | Ar / graphite wall | — | \(E\), \(F\) | **oracle=metalwalls** (Steele) or literature; map E8
- [ ] `born-repulsion-foothold` | Born repulsion wall/ion term | forcefield | ion near wall | — | \(E\), \(F\) | **oracle=metalwalls**; map E8 companion

### E9 — molecular water at electrodes (MD roadmap ○ virtual sites)

- [ ] `four-site-tip4p-foothold` | Four-site / virtual-site water (TIP4P-class) | md | TIP4P dimer or box | — | \(E\), geometry | **oracle=lammps** `tip4p` / virtual sites; map E9; aqueous electrode dep

### E10 — ecosystem (after E1 observers)

- [ ] `capacitance-charge-fluctuation` | Capacitance / charge-fluctuation analysis | electrochem | parallel-plate EDL | — | \(C\) estimate | **oracle=metalwalls** tools; map E10; **dep:** E1 + charge observer
- [ ] `potential-profile-observer` | Electrostatic potential profile (slab) | electrochem | parallel plate | — | \(\Psi(z)\) curve | **oracle=metalwalls** / literature Poisson; map E10; **dep:** charges + long-range

### Shared classical MD ○ still open (LAMMPS primary; electrochem-adjacent)

- [ ] `constant-force-pull-ion` | Constant-force / steered ion pull (noneq MD ○) | md | ion between walls | — | \(\langle z\rangle\), work proxy | **oracle=lammps** `fix smd` / `fix addforce`; useful EDL transport twin — **not** MW-only
- [ ] `surface-tension-coupling` | Constant surface-tension ensemble ○ | md | LJ slab | — | \(\gamma\) / \(L_z\) | **oracle=lammps** `fix npt` surface tension path; optional for planar interfaces

## in_progress

(empty)

## blocked

- [!] `cuda-const-pot` | blocked — CUDA electrode track is separate; this loop is CPU MD-only
- [!] `aim-daim-product` | blocked until map E11 is explicitly scheduled; not E1
- [!] `ml-pinnwall-electrode` | blocked — map E12 research; external ML interface
- [!] `periodic-hybrid-electrode-dft` | blocked — quantum electrodes → CP2K/PySCF loops, not classical MW loop

## done

(seed later from harness)
