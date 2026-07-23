<!-- mol:note:topic:metalwalls-capability-map -->
## [2026-07-23] MetalWalls → Atomiverse capability map

**Why.** MetalWalls is a full electrochemical classical MD code, not “solve
\(AQ=b\) once.” Agents must not under-scope (claim MW = constant-potential
kernel only) or over-scope (rewrite ATV’s MD stack to clone Fortran
`system`). This page is the canonical inventory and delivery slicing.

**Rule (binding for electrochemistry / constant-potential work):**

1. **Map first.** Any plan or spec that cites MetalWalls must state which
   **slice** (E1–E12 below) it covers and which map rows it intentionally
   leaves out.
2. **Architecture.** Constant-potential electrodes are a **Component** on
   existing `cpu::md::MD` (and later CUDA MD if needed) — **not** a new
   Driver category, not a public `*Config` bag, not logic welded into
   `KSpace*::compute`. Parameters = ctor args / members; \(A\), \(Q\), CG
   workspace = private; charges write `frame.span(kAtom_charge)`.
3. **E1 is the minimum landable product**, not the full MW feature set.
4. **Do not reimplement** ATV strengths already marked ATV+ (PME/CUDA,
   GAFF bonded stack, SCF/PW/GPW, MC/path methods) “for parity.”

**Evidence.** MW: GitLab `ampere2/metalwalls` branch `release` (src +
`example/LIST`), JOSS 2020 (10.21105/joss.02373), JCP 2022 polarizable
interfaces. ATV: tree as of 2026-07-23 (`docs/roadmap/force-fields.md`,
`cpu/{dynamics,potentials,longrange,bias,observe}/`).

**Supersedes:** conversation-only summaries that treated MW as “just
constant-potential BO.”

---

### Legend

| Mark | Meaning |
|:---:|---|
| ✓ | Present, production-usable |
| ~ | Partial / different form / not the same physics |
| ○ | Absent |
| ★ | MetalWalls specialty (electrochemical depth) |
| ATV+ | Atomiverse stronger; MW weak or absent |

**Priority** (if targeting electrochemical MD): **P0** now · **P1**
production parity · **P2** enhancement · **P3** ecosystem · **—** non-goal
or already covered better in ATV.

---

### 0. Positioning

| Dimension | MetalWalls | Atomiverse |
|---|---|---|
| Role | Electrochemical-interface classical MD | Dual-backend MD + electronic structure |
| Depth | Constant potential / PIM / AIM / TF / piston / 2D | Driver architecture, PME, force catalogue, CUDA, SCF/PW/GPW |
| Stack | Fortran + MPI + optional Python/PLUMED | C++17 + optional CUDA + molrs |
| Shape | Monolithic `system` + specialty modules | Driver / Component / Frame |

Parity ≠ cloning the program: **add electrochemical depth** on ATV’s MD
spine; keep ATV advantages.

---

### 1. Electrodes and constant potential ★

| Capability | MW | ATV | ATV home / gap | Pri |
|---|:---:|:---:|---|:---:|
| Constant-potential electrodes (fluctuating \(Q\)) | ✓★ | ○ | New Component → `kAtom_charge` | **P0** |
| Gaussian charge on electrode atoms | ✓★ | ○ | MD absent; `longrange/gaussian_extent.h` is SCF/CFMM only | **P0–P1** |
| Constant-charge electrodes (control) | ✓ | ~ | Fixed `charge` works; no electrode semantics | P1 |
| BO matrix inverse for \(Q\) | ✓ | ○ | Private workspace + factorize/invert | **P0** |
| CG solve for \(Q\) | ✓ | ○ | Large / movable electrodes | P1 |
| MaZe mass-zero charge dynamics | ✓★ | ○ | New algorithm | P2 |
| Electrode charge predictor | ✓ | ○ | ASPC-like; SCF ASPC is reference only | P1 |
| Electrode electroneutrality \(\sum Q_i=0\) | ✓ | ○ | Augmented linear system | **P0** |
| Global electroneutrality (incl. electrolyte) | ✓ | ○ | Constraint variant | P1 |
| Multi-electrode, per-electrode \(V,\eta\) | ✓ | ○ | Select by type/tag | **P0** |
| Forces on electrode atoms (term breakdown) | ✓ | ○ | Needed for piston; often zero if frozen | P1 |
| Freeze / immobile ions | ✓ (`mobile`) | ○ | Need `FreezeAtoms` or integrate mask | **P0** |
| Piston electrode + constant \(P_z\) | ✓★ | ○ | Rigid piston ≠ generic barostat | P2 |
| Thomas–Fermi finite metallicity | ✓★ | ○ | \(l_{\mathrm{TF}}\) response | P2 |
| Voronoi volume (TF) | ✓ | ○ | TF companion | P2 |
| Applied E-field (MD) | ✓★ | ○ | SCF has external field; MD does not | P1 |
| Displacement D-field | ✓★ | ○ | Finite-field electrochemistry | P2 |
| Finite-field + single electrode (3D) | ✓★ | ○ | Dufils-style setup | P2 |
| 2D PBC electrochemical cell | ✓★ | ○ | Ewald/min-image are 3D today | P1 |
| Charging work \(-\sum\Psi_i Q_i\) in energy | ✓ | ○ | `kEnergy` / meta channel | **P0** |
| Capacitance / charge-fluctuation analysis | ✓ tools | ○ | Observer + tools | P1 |
| Potential-profile post-processing | ✓ tools | ○ | Analysis | P1 |
| Charge–charge autocorrelation | ✓ tools | ○ | Analysis | P2 |
| ML heterogeneous electrodes (PiNNwall, …) | ecosystem ✓ | ○ | External interface | P3 |

**MD phase / halo (implementation constraint):** `kAtom_charge` already
borders (`kColComm|kColBorder|kColWriteBack`). Updating \(Q\) at `pre_force`
after `HaloSend` leaves **stale ghost charges** under MPI. Prefer solve
before halo, or `HaloRefresh` including charge after solve.

---

### 2. Polarizable electrolyte and floating charge ★

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| PIM induced dipoles (SCF) | ✓★ | ○ | Roadmap: self-consistent induced dipoles ○ | P1 |
| Extended-Lagrangian dipoles | ~ MaZe-related | ○ | Roadmap ○ | P2 |
| **Coupled** dipole ⊕ electrode \(Q\) solve | ✓★ | ○ | `dipoles_and_electrodes` — main hardness | **P1** |
| Thole / Tang–Toennies damping | ✓ | ○ | Roadmap Thole ○ | P1 |
| Per-species polarizability | ✓ | ○ | Frame column or component table | P1 |
| AIM aspherical ions | ✓★ | ○ | | P3 |
| DAIM deformable ions | ✓★ | ○ | | P3 |
| Multipole force field | ✓ | ~ | `observe/multipole` ≠ multipole FF | P2 |
| Fluctuating-charge / QEq | ~ electrode only | ○ | Roadmap QEq/FQ ○; ≠ const-potential electrode | P2 |
| Drude / AMOEBA | ○ | ○ | Both roadmap-level | — |
| QM/MM polarizable | ○ | ~ | `qmmm_polarizable` — ATV path, not MW | ATV+ |

---

### 3. Classical nonbonded

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| LJ 12-6 | ✓ | ✓⚡ CPU/CUDA | ATV stronger | — |
| LJ soft-core | ✓ | ✓ `PairLJSoft` | | — |
| Coulomb cut | ✓ | ✓ `PairCoulCut` | | — |
| Coulomb soft-core | ✓ | ✓ `PairCoulSoft` | | — |
| Ewald real + reciprocal | ✓ | ✓ `CoulLong` + `KSpaceEwald` | | — |
| PME / SPME | ○/~ | ✓ CPU/CUDA | **ATV+** | ATV+ |
| PPPM | ○ | ✓ | | ATV+ |
| 2D Ewald / 2D PBC | ✓★ | ○ | Electrochem default in MW | P1 |
| \(k=0\) slab correction | ✓ | ○ | With 2D | P1 |
| Fumi–Tosi | ✓★ | ○ | Molten salts | P2 |
| XFT | ✓ | ○ | | P3 |
| Born repulsion | ✓ | ○ | | P2 |
| Steele wall | ✓★ | ○ | Graphite wall | P2 |
| Morse / Yukawa / WCA / DPD | ○ | ✓ | | ATV+ |
| EAM / Stillinger–Weber | ○ | ✓ | | ATV+ |
| `special_bonds` scales | ~ | ✓ | GAFF path mature | ATV+ |

---

### 4. Bonded / molecular models

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| Harmonic bond / angle | ✓ | ✓ CPU/CUDA | | — |
| Dihedrals | ✓ | ✓ multi/harmonic + CHARMM | ATV broader | ATV+ |
| Impropers | ✓ | ✓ harmonic | | — |
| FENE | ○ | ✓ | | ATV+ |
| Four-site water (TIP4P-class) | ✓★ | ○ | Virtual site | P1 (aqueous electrodes) |
| Linear molecule model | ✓ | ○ | | P3 |
| Molecule-type table | ✓ | ~ | ATV per-atom topology columns | ~ |
| Virtual sites | ✓ four-site | ○ | | P1 |

---

### 5. Constraints / integration / ensembles

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| Velocity Verlet-class | ✓ | ✓ CPU/CUDA | | — |
| RATTLE | ✓ | ✓ `FixRattle` | | — |
| SETTLE | ○ | ✓ | | ATV+ |
| Nose–Hoover chain NVT | ✓ | ✓ `NoseHoover` | | — |
| Langevin / CSVR | ○ | ✓ | | ATV+ |
| NPT (MTK-class) | ✓ | ✓ `NoseHooverNPT` etc. | | — |
| NPH / PR NPT | ○ | ✓ | | ATV+ |
| Piston NPT (electrodes) | ✓★ | ○ | | P2 |
| Velocity init / scale | ✓ | ✓ | | — |
| Deform / Brownian / NVE limit | ○ | ✓ | | ATV+ |

---

### 6. Long-range and communication

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| 3D Ewald | ✓ | ✓ | | — |
| MPI MD + domain decomp | ✓ | ✓ | | — |
| Neighbor lists | ✓ | ✓ | | — |
| Charge halo correctness after \(Q\) update | ✓ | ✓ schema; **timing** must be designed | See §1 | **P0** |
| CUDA MD + PME | ○ | ✓ | **ATV+** | ATV+ |
| FMM open-boundary | ○ | ~ operators / SCF | Not MD `KSpace` | — |

---

### 7. Driver / workflow

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| Dedicated MD loop | ✓ monolith | ✓ `md::MD` nestable | | ATV+ |
| File input language | ✓ runtime/data | ~ library API | | — |
| Python driving | ✓ f2py | ~ nanobind wheel | Different shape | P3 |
| Thermodynamic integration workflow | ✓★ + soft-core | ~ soft-core only | | P2 |
| PLUMED | ✓ | ○ | ATV has in-tree bias/metadynamics | ~ |
| Collective vars / metadynamics | via PLUMED | ✓ `cpu/bias/*` | | ATV+ |
| Continuable runs | ✓ | ✓ | | — |

---

### 8. Observe / I/O / post-process

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| Thermo / trajectories | ✓ | ✓ + molrs | | — |
| RDF / MSD / VACF | ~ | ✓ | | ATV+ |
| Virial / pressure tensor | ✓ | ✓ | | — |
| Capacitance scripts | ✓★ | ○ | | P1 |
| Empty capacitance / matrix condition number | ✓ tools | ○ | Debug \(A\) | P1 |
| Electrode charge time series | ✓ | ○ | Observer | **P0** |
| Widom / WHAM | ○ | ✓ | | ATV+ |
| Electronic observers (Mulliken, …) | ○ | ✓ | | ATV+ |

---

### 9. Scenario coverage (MW `example/`)

| Scenario | MW | ATV today | Gap |
|---|:---:|---|---|
| Bulk non-polarizable IL | ✓ | ✓ LJ+Coul+PME | — |
| Bulk polarizable IL / PIM-NPT | ✓ | ○ | PIM |
| Planar graphene + const-pot + IL | ✓ | ○ | E1 + preferably 2D |
| Nanoporous carbon + const-pot | ✓ | ○ | E1 + geometry |
| Blue energy / CDC desalination | ✓ | ○ | E1 |
| Finite-field single electrode | ✓ | ○ | E-field + electrode |
| TF gold capacitor | ✓ | ○ | Thomas–Fermi |
| Molten salt + Al + MaZe | ✓ | ○ | Const-pot + MaZe + salt FF |
| Piston constant \(P_z\) | ✓ | ○ | Piston |
| PLUMED constrained ion | ✓ | ~ ATV bias | |
| TIP4P water | ✓ | ○ no four-site | Virtual sites |
| AIM tests | ✓ | ○ | AIM |
| AIMD / Gaussian SCF / PW / GPW | ○ | ✓ | **ATV+** |
| CUDA large LJ/PME | ○ | ✓ | **ATV+** |
| NEB / REMD / PIMD / MC | ○ | ✓ | **ATV+** |

---

### 10. Delivery slices (use these in specs)

| Slice | Map coverage | Goal |
|---|---|---|
| **E1 const-pot v1** | §1 P0: BO \(Q\), multi-electrode \(V\), neutrality, freeze, charging work, Observer, bind `md::MD` | Parallel-plate EDL runs |
| **E2 electrostatic rigor** | Gaussian kernels, 2D Ewald, predictor, CG | Numerical near MW small cases |
| **E3 external field** | E-field, finite-field single electrode | Dufils-class setups |
| **E4 polarization** | PIM + damping | Bulk polarizable IL |
| **E5 coupling** | Dipole ⊕ electrode joint solve | JCP 2022 main capability |
| **E6 electrode mechanics** | Piston \(P_z\), electrode forces | Piston examples |
| **E7 metallicity** | Thomas–Fermi | Imperfect metal |
| **E8 salt / walls** | Fumi–Tosi, Steele, Born | LiCl–Al, graphite wall |
| **E9 molecular water** | Four-site / virtual sites | TIP4P at electrodes |
| **E10 ecosystem** | Capacitance tools, TI workflow, optional PLUMED bridge | Paper workflows |
| **E11 AIM/DAIM** | §2 AIM | Deep molten-salt users only |
| **E12 ML electrodes** | PiNNwall-class | Research frontier |

**E1 only** is the “simple” kernel discussed in chat — roughly the P0 rows of
§1, **not** the full map.

---

### 11. Forbidden when implementing from this map

- New `ElectrochemistryDriver` / category base  
- Public `ElectrodeConfig` / `ConstPotentialConfig` bags  
- Absorbing charge solve into `KSpaceEwald` / `KSpacePME`  
- Claiming MetalWalls equivalence from point-charge \(A\) without documenting
  the approximation  
- Ignoring MPI ghost charge after \(Q\) update  
- Full Fortran `system` port or dual MD architectures

---

### Delivery workflow

- Loop: `.grok/workflows/metalwalls-method-loop.rhai` (`/metalwalls-method-loop`)
- Design / queue / handoff: `docs/developer/metalwalls-method-*.md`
- Oracle: `tools/metalwalls_oracle/` — **LAMMPS first** when feature overlaps;
  MetalWalls only for MW-only physics
- Build: **MD-only** (`ATV_BUILD_MD`+`KSPACE`+tests/regressions; HF/PW/GPW/CUDA/Torch OFF)
- Default parallel workers: **staff_count=10**

### References (stable)

- Repo: https://gitlab.com/ampere2/metalwalls (`release`)  
- JOSS: Marin-Laflèche et al., JOSS 5, 2373 (2020)  
- JCP 2022: polarizable electrolyte + metallic electrodes  
- ATV roadmap polarizable/FQ rows: `docs/roadmap/force-fields.md`  

<!-- mol:note:topic:metalwalls-capability-map -->
