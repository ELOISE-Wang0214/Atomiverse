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
5. **Foothold ≠ closed.** A textbook regression + Component does **not**
   flip a row to ✓ if production depth is still missing (e.g. O(N²)
   pairs without NeighborList, Yeh–Berkowitz slab vs full 2D-PBC product).
   Use **~** until production depth lands; open a deeper slug.

**Evidence.** MW: GitLab `ampere2/metalwalls` branch `release` (src +
`example/LIST`), JOSS 2020 (10.21105/joss.02373), JCP 2022 polarizable
interfaces. ATV: tree as of **2026-07-24 census wave 5** on `dev`
(waves 1–3 shipped 26 classical electrochem/MD/FF slugs; wave-4 staff
unstarted — no new Components/regressions; see
`docs/developer/metalwalls-method-queue.md`).

**Supersedes:** conversation-only summaries that treated MW as “just
constant-potential BO." Map ATV marks re-audited **2026-07-24** (wave-5
gap-census): marks unchanged vs wave 4 (no landings); pending refill added
per-electrode \(\eta\), const-charge control, MPI charge-halo, Voronoi TF,
PIM-NPT IL, IL EDL box, nanoporous product.

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
| Constant-potential electrodes (fluctuating \(Q\)) | ✓★ | ✓ | `ConstPotBO` / `ConstPotCG` → `kAtom_charge` | **P0** |
| Gaussian charge on electrode atoms | ✓★ | ~ | `PairGaussianCoul` + Gaussian \(A\) in BO; **O(N²)** pair path | **P0–P1** |
| Constant-charge electrodes (control) | ✓ | ~ | Fixed charge + `FixFreeze`; no separate “const-Q electrode” type | P1 |
| BO matrix inverse for \(Q\) | ✓ | ✓ | Dense Cholesky path on `ConstPotBO` | **P0** |
| CG solve for \(Q\) | ✓ | ✓ | `ConstPotCG` / `ConstPotBO::use_cg` | P1 |
| MaZe mass-zero charge dynamics | ✓★ | ○ | No `MaZe` charge dynamics Component yet | P2 |
| Electrode charge predictor | ✓ | ✓ | `ChargeASPC` / `predictor_order` on BO | P1 |
| Electrode electroneutrality \(\sum Q_i=0\) | ✓ | ✓ | Augmented system on BO/CG | **P0** |
| Global electroneutrality (incl. electrolyte) | ✓ | ~ | `ConstPotentialBO` form; unify on `ConstPotBO` path | P1 |
| Multi-electrode, per-electrode \(V,\eta\) | ✓ | ~ | Multi-type \(V\) ✓; shared \(\eta\) in v1 | **P0** |
| Forces on electrode atoms (term breakdown) | ✓ | ○ | Piston uses fluid \(F_z\); no full BO \(Q\)-force breakdown | P1 |
| Freeze / immobile ions | ✓ (`mobile`) | ✓ | `FixFreeze` | **P0** |
| Piston electrode + constant \(P_z\) | ✓★ | ✓ | `PistonElectrode` | P2 |
| Thomas–Fermi finite metallicity | ✓★ | ✓ | `ConstPotBO::l_tf`, `site_density` | P2 |
| Voronoi volume (TF) | ✓ | ○ | TF companion still absent | P2 |
| Applied E-field (MD) | ✓★ | ✓ | `FixEfield` + `FixOscillatingEfield` | P1 |
| Displacement D-field | ✓★ | ○ | Finite-field \(\mathbf{D}\) not shipped | P2 |
| Finite-field + single electrode (3D) | ✓★ | ✓ | Dufils path (`ConstPotBO` + \(E_z\)) | P2 |
| 2D PBC electrochemical cell | ✓★ | ~ | Yeh–Berkowitz **slab** only; full 2D product **out of loop scope** | P1 |
| Charging work \(-\sum\Psi_i Q_i\) in energy | ✓ | ✓ | `ChargingWork` | **P0** |
| Capacitance / charge-fluctuation analysis | ✓ tools | ✓ | `CapacitanceChargeFluctuation` | P1 |
| Potential-profile post-processing | ✓ tools | ✓ | `PotentialProfile` \(\Psi(z)\) | P1 |
| Charge–charge autocorrelation | ✓ tools | ○ | No charge–charge \(C_{QQ}(t)\) observer | P2 |
| Empty capacitance / matrix condition | ✓ tools | ○ | No \(\mathrm{cond}(A)\) / empty-\(C\) debug observer | P1 |
| ML heterogeneous electrodes (PiNNwall, …) | ecosystem ✓ | ○ | External interface — **blocked E12** | P3 |

**MD phase / halo (implementation constraint):** `kAtom_charge` already
borders (`kColComm|kColBorder|kColWriteBack`). Updating \(Q\) at `pre_force`
after `HaloSend` leaves **stale ghost charges** under MPI. Prefer solve
before halo, or `HaloRefresh` including charge after solve.

---

### 2. Polarizable electrolyte and floating charge ★

| Capability | MW | ATV | Notes | Pri |
|---|:---:|:---:|---|:---:|
| PIM induced dipoles (SCF) | ✓★ | ✓ | `PimInducedDipole` bulk foothold | P1 |
| Extended-Lagrangian dipoles | ~ MaZe-related | ○ | Roadmap ○ — deeper after SCF PIM | P2 |
| **Coupled** dipole ⊕ electrode \(Q\) solve | ✓★ | ✓ | `DipoleElectrodeJoint` | **P1** |
| Thole / Tang–Toennies damping | ✓ | ✓ | Thole + optional `tt_b` on PIM | P1 |
| Per-species polarizability | ✓ | ✓ | `alpha_by_type` table on PIM | P1 |
| AIM aspherical ions | ✓★ | ○ | **blocked E11** until scheduled | P3 |
| DAIM deformable ions | ✓★ | ○ | **blocked E11** | P3 |
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
| 2D Ewald / 2D PBC | ✓★ | ~ | Yeh–Berkowitz slab foothold; full 2D **out of scope** this loop | P1 |
| \(k=0\) slab correction | ✓ | ✓ | `KSpaceEwald::slab_factor` | P1 |
| Fumi–Tosi | ✓★ | ~ | `PairBorn` multi-type + Coul; **O(N²)** textbook pin | P2 |
| XFT | ✓ | ○ | | P3 |
| Born repulsion | ✓ | ~ | `PairBorn` dimer foothold; needs NeighborList bulk | P2 |
| Steele wall | ✓★ | ✓ | `FixWallLJ1043` (`wall/lj1043`) | P2 |
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
| Four-site water (TIP4P-class) | ✓★ | ✓ | `PairLJCutTip4pCut` + `VirtualSites` Weighted3 | P1 (aqueous electrodes) |
| Linear molecule model | ✓ | ○ | | P3 |
| Molecule-type table | ✓ | ~ | ATV per-atom topology columns | ~ |
| Virtual sites | ✓ four-site | ✓ | Weighted3 + OutOfPlane3 (TIP4P/TIP5P) | P1 |

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
| Piston NPT (electrodes) | ✓★ | ✓ | `PistonElectrode` (rigid \(P_z\); ≠ Langevin piston) | P2 |
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
| Capacitance scripts | ✓★ | ✓ | `CapacitanceChargeFluctuation` | P1 |
| Empty capacitance / matrix condition number | ✓ tools | ○ | Debug \(A\) / empty \(C\) still open | P1 |
| Electrode charge time series | ✓ | ✓ | `ElectrodeCharge` observer | **P0** |
| Widom / WHAM | ○ | ✓ | | ATV+ |
| Electronic observers (Mulliken, …) | ○ | ✓ | | ATV+ |

---

### 9. Scenario coverage (MW `example/`)

| Scenario | MW | ATV today | Gap |
|---|:---:|---|---|
| Bulk non-polarizable IL | ✓ | ✓ LJ+Coul+PME | — |
| Bulk polarizable IL / PIM-NPT | ✓ | ~ | PIM bulk ✓; NPT product twin open |
| Planar graphene + const-pot + IL | ✓ | ~ | E1 Components ✓; production EDL box open |
| Nanoporous carbon + const-pot | ✓ | ○ | Geometry + E1 (product) |
| Blue energy / CDC desalination | ✓ | ○ | Product workflow |
| Finite-field single electrode | ✓ | ✓ | Dufils foothold |
| TF gold capacitor | ✓ | ✓ | `thomas-fermi-electrode` |
| Molten salt + Al + MaZe | ✓ | ○ | MaZe + salt NL depth |
| Piston constant \(P_z\) | ✓ | ✓ | `piston-electrode-pz` |
| PLUMED constrained ion | ✓ | ~ ATV bias | PLUMED bridge **out of scope** |
| TIP4P water | ✓ | ✓ | four-site + multi-topology vsites |
| TIP4P + const-pot electrodes | ✓ | ○ | Product aqueous EDL still open |
| AIM tests | ✓ | ○ | **blocked E11** |
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
