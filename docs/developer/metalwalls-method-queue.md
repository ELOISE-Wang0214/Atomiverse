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

Census: 2026-07-23 (seed). Parent `dev`. **No incomplete yet** — seed pending
from map E1 footholds + shared MD twins.

**Recover reconfirm 2026-07-23:** incomplete_count=**0** for this loop.
- `## in_progress` `[~]`: empty
- live `.claude/specs/` (MD/electrochem/const-pot/forcefield): none
- unmerged `wf/<slug>` in this loop’s scope: none (local `wf/*` tips are
  PySCF/post-HF ancestors of `dev`, out of scope; no metalwalls/const-pot
  branches local or on `origin`)
- handoff incomplete rows: 0
- half-written `regressions/md|forcefield/<slug>/` for *this* queue: none
  (REFERENCE-less but shipped method-ladder cases noted under orphans only)

## pending

### E1 — constant-potential v1 (map §1 P0)

- [ ] `const-pot-bo-parallel-plate` | Constant-potential BO solve (parallel plate) | electrochem | NaCl slab + two walls | — | electrode Q / ΔΨ | **oracle=metalwalls** (or LAMMPS ELECTRODE if package available offline); map E1; freeze electrodes; electroneutrality; write `kAtom_charge`
- [ ] `const-pot-charge-observer` | Electrode total charge observer | electrochem | same foothold | — | Q_tot(t) | **oracle=lammps** if fix electrode charge dump; else metalwalls; dep: charge solve
- [ ] `freeze-electrode-atoms` | Freeze / immobile electrode atoms | md | any | — | x fixed | **oracle=lammps** `fix setforce 0` / `fix freeze`; Component on integrate path
- [ ] `charging-work-energy` | Charging work −Σ Ψ_i Q_i in energy | electrochem | parallel plate | — | E channel | **oracle=metalwalls** literature/Scalfi; dep: const-pot
- [ ] `const-pot-md-nvt-short` | Const-pot + MD NVT short traj | electrochem | NaCl + walls | — | stable E, Q sign | **oracle=metalwalls** or lammps electrode MD; dep: BO + freeze + CoulLong/PME

### Shared classical (LAMMPS primary)

- [ ] `lj-cut-lammps-dimer` | Pair LJ cut energy/force | forcefield | Ar dimer | — | E, F | **oracle=lammps** `pair_style lj/cut`; regression pin; may already exist — re-pin under this loop only if missing REFERENCE
- [ ] `coul-long-pme-nacl-lammps` | CoulLong + PME NaCl | forcefield | rock-salt | — | E | **oracle=lammps** `lj/cut/coul/long` + `kspace pppm/pme`; dep: KSPACE
- [ ] `nose-hoover-nvt-lammps` | Nose–Hoover NVT thermostat | md | LJ fluid | — | ⟨T⟩ | **oracle=lammps** `fix nvt`; statistical tol
- [ ] `special-bonds-amber-lammps` | special_bonds 1-4 scale | forcefield | PE fragment | — | E | **oracle=lammps** special_bonds; GAFF path

### E2 / E3 footholds (after E1)

- [ ] `const-pot-cg-large-electrode` | CG electrode charge solve | electrochem | larger wall | — | Q | **oracle=metalwalls**; dep: BO dense inverse path
- [ ] `external-efield-md` | Applied E-field on MD charges | electrochem | ions in box | — | force on q | **oracle=lammps** `fix efield` when overlap; else metalwalls external_field
- [ ] `slab-2d-ewald-foothold` | 2D / slab Ewald foothold | electrochem | slab | — | E | **oracle=lammps** `kspace_modify slab` if available; else metalwalls 2D; map E2

### Later slices (do not assign until deps green)

- [ ] `pim-induced-dipole-bulk` | PIM induced dipoles bulk | polarizable | molten salt / IL | — | E_pol | **oracle=metalwalls**; map E4
- [ ] `dipole-electrode-joint` | Dipole ⊕ electrode joint solve | electrochem | IL + graphene | — | Q, μ | **oracle=metalwalls**; map E5; **dep:** E1 + PIM
- [ ] `thomas-fermi-electrode` | Thomas–Fermi metallicity | electrochem | gold capacitor toy | — | C | **oracle=metalwalls**; map E7
- [ ] `piston-electrode-pz` | Piston electrode constant Pz | electrochem | IL + graphene | — | Pz | **oracle=metalwalls**; map E6
- [ ] `steele-wall-argon` | Steele wall potential | forcefield | Ar / graphite | — | E, F | **oracle=metalwalls** (Steele) or literature; map E8
- [ ] `fumi-tosi-nacl` | Fumi–Tosi molten salt | forcefield | NaCl | — | E | **oracle=metalwalls** / literature; map E8
- [ ] `four-site-tip4p-foothold` | Four-site / virtual site water | md | TIP4P | — | E | **oracle=lammps** tip4p; map E9

## in_progress

(empty at seed)

## blocked

- [!] `cuda-const-pot` | blocked — CUDA electrode track is separate; this loop is CPU MD-only
- [!] `aim-daim-product` | blocked until map E11 is explicitly scheduled; not E1
- [!] `ml-pinnwall-electrode` | blocked — map E12 research; external ML interface
- [!] `periodic-hybrid-electrode-dft` | blocked — quantum electrodes → CP2K/PySCF loops, not classical MW loop

## done

(seed later from harness)
