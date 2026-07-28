# SE / tight-binding method queue

Durable worklist for workflow `tb-method-loop`. Do not delete done rows.
Census rewrites **pending** from `docs/roadmap/semiempirical.md` ○ rows that
**OpenMOPAC / DFTB+ / xtb / literature** can oracle offline; preserve
done / blocked / in_progress.

**Build rule:** SE-only cmake (`SEMIEMPIRICAL` + `EVALUATE` + TESTS/REGRESSIONS
**ON**; HF/PW/GPW/MD/MC/MINIMIZER/VIB/PATH/KSPACE/CUDA/TORCH **OFF**).
See `tb-method-loop.md`.

**Oracle rule (family-locked):**

| Family | Primary pin package |
| ------ | ------------------- |
| NDDO (MNDO/AM1/PM3/RM1/PM6/PM7) | **OpenMOPAC** — never DFTB+ |
| DFTB0 / DFTB2 / DFTB3 / spin / periodic DFTB | **DFTB+** — never OpenMOPAC |
| GFN-xTB | **xtb** — not a DFTB3 successor; schedule late |
| Unified API design (no energy pin) | SCINE Sparrow (architecture) |
| Driver embed / MD / periodic patterns | CP2K (architecture only) |

`REFERENCE.md` must state `oracle: mopac|dftb+|xtb|literature` and
`oracle_version`.

**Assign priority (binding ladder):** DFTB0 → DFTB2 depth → DFTB3 → spin DFTB
→ NDDO PM ladder → **only then** GFN-xTB. Do not staff GFN2 while DFTB0/3 and
analytic DFTB forces are still open unless operator `theme=` overrides.
**Recover incomplete `[~]` before new pending** (`incomplete ≫ pending`).

## GapCensus — 2026-07-27 wave-3 (parent `dev`)

**Branch:** `dev` @ `cec4f6ae` (SE product surface unchanged since harness wave-2
`ddc1ccc1`; HEAD is concurrent esmm/driver land).

### Sources re-diffed

| Source | Finding |
| ------ | ------- |
| `docs/roadmap/semiempirical.md` ○ | PM3/6/7/OM3; closed-form NDDO multipole F; Hessian; DFTB0/3/D3/stress/TD/noncollinear/+U/SOC; GFN1/2/IPEA/spin/PBC/analytic G/H; GFN-FF; sTDA |
| `include/atomiverse/cpu/semiempirical/` | Public: `AM1`/`RM1` (+evaluators), `DFTB2`, `XtbGfn0`. **No** `DFTB0`/`DFTB3`/`PM*`/GFN1/2 types on `dev` |
| `SemiEmpiricalMethod` | **AM1\|RM1 only** — PM tables still need enum + param rows |
| `basis_data/` | mio-style sketch H/C/N/O + mio-1-1 H–H pin + 3ob-3-1 CHNO |
| `NDDOSolver::compute_forces` | **fixed-P central FD** of E (product path ✓; closed-form multipole dI/dR still ○) |
| `dftb2.cpp` forces on `dev` | Restricted MioSketch: **closed-form** HF+Pulay+γ+rep; open-shell Mio + all ThreeOb: **re-SCC FD** |
| Analytic F under Γ-PBC | Energy uses lattice image sums; `analytic_forces_restricted` pair walk is **min_image only** → deepen `dftb2-periodic-forces` |
| Spin path | Collinear `multiplicity` on MioSketch; ThreeOb rejects open shell |
| `xtb_gfn0` | H/C/O toy — not production GFN0 / xtb |
| Regressions on `dev` | `am1_*`, `rm1_water_energy`, `h2_dftb2`, `h2o_dftb2_mio`, `dftb2-3ob-organic`, `dftb2-analytic-forces`, `dftb2-dftbplus-energy-pin`, `dftb2-periodic-gamma`, `spin-dftb2-o2`, `h2o_gfn0` |
| REFERENCE hygiene | Wave-1/2 DFTB regs have `oracle:`; `h2_dftb2` **no** REFERENCE; `h2o_gfn0` **no** REFERENCE; AM1/RM1 lack hard mopac E pin + often lack `oracle:` line |
| `tools/tb_oracle/` on `dev` | `ch4_3ob_dftb2`, `h2_mio_scc`, `dftb2-periodic-gamma`, force notes; **no** `inputs/mopac/` or `inputs/xtb/` trees |
| OpenMOPAC surface | MNDO / AM1 / PM3 / RM1 / PM6 / PM7 / MNDO-d — **not** OM3, **not** DFTB |
| DFTB+ surface | non-SCC, SCC-DFTB2, DFTB3, spin, PBC/k, forces/stress, TD-DFTB, +U, SOC, dispersion, … |
| xtb surface | GFN0/1/2, IPEA, analytic G; GFN-FF classical (scope-out); sTDA needs GFN GS |
| Live SE specs on `dev` | **none** |
| TB isolation worktrees | **none** (Clearup clean; only unrelated `driver-api-02` tree) |
| Unmerged TB tips | **6** with commits ahead of `dev` (3 older ladder tips behind=41; 3 force/stress tips behind=25) |

### Tree vs roadmap (wave-3)

| Capability | Roadmap | Tree on `dev` | Verdict |
| ---------- | ------- | ------------- | ------- |
| AM1 / RM1 energy + ΔHf | ✓ | ship; no hard mopac E pin | deepen → mopac pins |
| NDDO nuclear forces (product) | ✓ | fixed-P central FD | usable; **not** closed-form multipole |
| SCC-DFTB2 H/C/N/O E | ✓ | MioSketch + ThreeOb + mio-1-1 H–H pin | energy product path |
| DFTB analytic gradient | ✓ (restricted) | closed-form on **restricted MioSketch only** | deepen ThreeOb + open-shell F (tips) |
| Spin-polarized DFTB | ✓ (collinear) | MioSketch dual αβ + Mulliken m | foothold; open-shell F tip; noncollinear blocked |
| Periodic DFTB | ✓ (Γ energy) | lattice-sum H⁰/S/γ/rep on MioSketch | **landed**; k-mesh + lattice-consistent F still open |
| DFTB analytic stress | ○ | absent on `dev` | tip `wf/dftb2-analytic-stress` |
| Non-SCC DFTB0 | ○ | tip only | **incomplete** |
| DFTB3 | ○ | tip only | **incomplete** |
| Dispersion D3 | ○ | tip only | **incomplete** |
| TD-DFTB | ○ | absent | **gap** |
| GFN0-xTB | ✓ toy | sketch H/C/O | toy → Phase C |
| GFN1/2 / IPEA / spin / PBC GFN | ○ | absent | **gaps** (late) |
| PM3 / PM6 / PM7 | ○ | no types | Phase B **oracle=mopac** |
| OM3 | ○ | no type | **blocked** |
| GFN-FF / sTDA / SE Hessian (VIB) | ○ | — | **blocked** (scope or SE-only cmake) |

### Foothold ≠ closed (wave-3)

- Three force/stress pending slugs now have unmerged tips → **promote to `[~]`**
  (`dftb2-analytic-forces-threeob`, `dftb2-spin-analytic-forces`,
  `dftb2-analytic-stress`); do **not** leave them as pending while tips exist.
- Γ-PBC energy ships with `with_forces=false` pin; analytic F pair geometry is
  min_image — lattice-image forces under PBC remain a product gap → **new**
  `dftb2-periodic-forces`.
- Older ladder tips (`dftb0-h2`, `dftb3-water-mio`, `dftb2-dispersion-d3`) still
  need **Recover/rebase** (behind=41) before merge.
- NDDO hard OpenMOPAC E pins and closed-form multipole F still open (Phase B).

### gaps_remaining (wave-3 census snapshot)

**18** pending + **6** incomplete (`[~]`) + **9** blocked (scope-out /
forbidden included) → **nonzero; continue loop**. `new_slugs_added=1`
(`dftb2-periodic-forces`). Roadmap ○/✓ **unchanged** this wave (no SE product
merge onto `dev` since wave-2). **Recover incomplete ≫ staff new pending.**

## GapCensus — 2026-07-27 wave-2 (parent `dev`) — historical

**Branch:** `dev` @ `77d80486` at census time (pre later lands; see wave-3).

### Sources re-diffed

| Source | Finding |
| ------ | ------- |
| `docs/roadmap/semiempirical.md` ○ | PM3/6/7/OM3; closed-form NDDO F; Hessian; DFTB0/3/D3/PBC/stress/TD/noncollinear/+U/SOC; GFN1/2/IPEA/spin/PBC/analytic G/H; GFN-FF; sTDA |
| `include/atomiverse/cpu/semiempirical/` | Public types: `AM1`/`RM1` (+ evaluators), `DFTB2`, `XtbGfn0`. **No** `DFTB0`/`DFTB3`/`PM3`/`PM6`/`PM7`/GFN1/2 on `dev` |
| `SemiEmpiricalMethod` | Enum is **AM1\|RM1 only** — PM tables need enum + param rows |
| `basis_data/` | `dftb_mio_tables.h` mio-**style** sketch H/C/N/O; `dftb_3ob_tables.h` 3ob-3-1 organic CHNO |
| `NDDOSolver::compute_forces` | **fixed-density central FD** of E (not multipole dI/dR closed-form) |
| `dftb2.cpp` forces | Restricted MioSketch: **closed-form** HF+Pulay+γ+rep; open-shell Mio: **re-SCC FD**; ThreeOb: **re-SCC FD** |
| Spin path | Collinear `multiplicity` on MioSketch only; ThreeOb rejects open shell |
| `xtb_gfn0` | H/C/O toy onsite+hop+rep — not production GFN0 / xtb |
| Regressions on `dev` | `am1_*`, `rm1_water_energy`, `h2_dftb2`, `h2o_dftb2_mio`, `dftb2-3ob-organic`, `dftb2-analytic-forces`, `spin-dftb2-o2`, `h2o_gfn0` |
| REFERENCE hygiene | Wave-1 regs have `oracle:`; `h2_dftb2` **no** REFERENCE; `h2o_gfn0` **no** REFERENCE; AM1/RM1 lack hard mopac E pin |
| `tools/tb_oracle/` | Family README + `inputs/dftbplus/ch4_3ob_dftb2` + force note on `dev`; mopac/xtb trees empty; more twins only on `wf/*` tips |
| OpenMOPAC surface | MNDO / AM1 / PM3 / RM1 / PM6 / PM7 / MNDO-d — **not** OM3, **not** DFTB |
| DFTB+ surface | non-SCC, SCC-DFTB2, DFTB3, spin, PBC/k, forces/stress, TD-DFTB, +U, SOC, dispersion, … |
| xtb surface | GFN0/1/2, IPEA, analytic G; GFN-FF classical (scope-out); sTDA needs GFN GS |
| Live SE specs on `dev` | **none** (wave-1 + wave-2 harness specs deleted) |
| TB isolation worktrees | **none** (Clearup done); metalwalls `wf/*` worktrees unrelated |
| Unmerged TB tips | **3** still ahead of `dev` by 1–2 commits, **behind ~39** after wave-2 land (need rebase before merge) |

### Tree vs roadmap (wave-2)

| Capability | Roadmap | Tree on `dev` | Verdict |
| ---------- | ------- | ------------- | ------- |
| AM1 / RM1 energy + ΔHf | ✓ | ship; no hard mopac E pin | deepen → mopac pins |
| NDDO nuclear forces (product) | ✓ | fixed-P central FD | usable; **not** closed-form multipole |
| SCC-DFTB2 H/C/N/O E | ✓ | MioSketch + ThreeOb + **mio-1-1 H–H pin** | energy product path; H–H absolute pin **done** wave-2 |
| DFTB analytic gradient | ✓ (restricted) | closed-form on **restricted MioSketch only** | deepen ThreeOb + open-shell F |
| Spin-polarized DFTB | ✓ (collinear) | MioSketch dual αβ + Mulliken m | foothold; open-shell F still FD; noncollinear blocked |
| Production 3ob organic SK | — (prose) | `DftbSkSet::ThreeOb` + CH4 pin | landed wave-1 (DFTB2 SK, not DFTB3) |
| Non-SCC DFTB0 | ✓ | `dftb0.h`; reads the transcribed mio-1-1 H–H tables | **landed** 2026-07-28 (`757155d0`+`2c1502b6`); pin = the DFTB+ 24.1 value, licensed by SCC vanishing for homonuclear H₂ |
| DFTB3 | ✓ | `dftb3.h` + `internal/dftb3_scc.h` | **landed** 2026-07-28 (`68df18a8`..`148a17ed`); third-order potential fixed to ∂E³/∂Δq (Euler 4.4e-16), SCC now converges (Anderson mixing) |
| Dispersion D3 | ✓ | D3(BJ) applied in `DFTB2::launch` with its own analytic gradient | **landed** 2026-07-28 (`53ff52a0`); E_disp vs s-dftd3 1.5e-12 eV; absent from `kStress` (gap) |
| Periodic DFTB | ✓ (Γ) | lattice-sum H⁰/S/γ/rep on MioSketch | **landed** wave-2; k-mesh → `dftb2-kpoints` |
| DFTB analytic stress | ✓ | closed-form virial, restricted MioSketch | **landed** 2026-07-28 (`24d43574`); analytic vs strain-FD 1.1e-10 eV/Å³ |
| TD-DFTB | ○ | absent | **gap** |
| GFN0-xTB | ✓ toy | sketch H/C/O | toy → Phase C deepen |
| GFN1/2 / IPEA / spin / PBC GFN | ○ | absent | **gaps** (late) |
| PM3 / PM6 / PM7 | ○ | no types | Phase B **oracle=mopac** |
| OM3 | ○ | no type | **blocked** |
| GFN-FF / sTDA / SE Hessian (VIB) | ○ | — | **blocked** (scope or SE-only cmake) |

### Foothold ≠ closed (new deepen slugs this wave)

Wave-1 lands closed restricted Mio forces + collinear spin **energy**/mulliken
gates + 3ob energy pin — **not** product force completeness:

- ThreeOb nuclear forces remain re-SCC FD → `dftb2-analytic-forces-threeob`
- Open-shell nuclear forces remain re-SCC FD → `dftb2-spin-analytic-forces`
- DFTB3 energy tip exists but forces product still open → `dftb3-analytic-forces`
- Γ-point tip exists; DFTB+ k-mesh surface still open → `dftb2-kpoints`
- RM1 ships without OpenMOPAC hard E pin → `rm1-water-mopac-energy`

### gaps_remaining (wave-2 census snapshot; post harness wave-2 land)

**20** pending + **3** incomplete (`[~]`) + **9** blocked (includes scope-out /
forbidden-pattern rows) → **nonzero; continue loop**. Never claim queue
exhausted. Wave-2 harness closed `dftb2-dftbplus-energy-pin` +
`dftb2-periodic-gamma` → done. Next: **Recover** the three remaining `[~]`
tips (rebase onto `dev`) before staffing new Phase A pending.

## GapCensus — 2026-07-27 wave-1 (parent `dev`) — historical

**Branch:** `dev` @ `d80fd953` at census time (pre wave-1 land). Post-harness
HEAD is later; see `## done` / handoff.

### Sources re-diffed

| Source | Finding |
| ------ | ------- |
| `docs/roadmap/semiempirical.md` ○ | NDDO PM6/PM7/OM3/Hessian; DFTB ladder beyond energy foothold; full GFN product |
| `include/atomiverse/cpu/semiempirical/` | Types: `AM1`/`RM1` + evaluators, `DFTB2`, `XtbGfn0` toy; **no** `DFTB0`/`DFTB3`/`PM6`/`PM7`/GFN2 (at census) |
| `basis_data/dftb_mio_tables.h` | mio-**style** sketch (not DFTB.org mio-1-1 bit-identical); H/C/N/O only |
| `src/.../nddo_solver.cpp` `compute_forces` | **fixed-density central FD** of E (not closed-form multipole dI/dR) |
| `src/.../dftb2.cpp` `with_forces` | was re-SCC central FD; **now** closed-form restricted (post wave-1) |
| `src/.../xtb_gfn0.cpp` | H/C/O onsite+hop+rep **sketch** — not production GFN0/xtb |
| Regressions | `am1_*`, `rm1_water_energy`, `h2_dftb2`, `h2o_dftb2_mio`, + wave-1: `dftb2-3ob-organic`, `dftb2-analytic-forces`, `spin-dftb2-o2` |
| `tools/tb_oracle/` | README family matrix; offline input twins growing on unmerged tips |
| OpenMOPAC surface | MNDO / AM1 / PM3 / RM1 / PM6 / PM7 / MNDO-d — **not** OM3, **not** DFTB |
| DFTB+ surface | non-SCC, SCC-DFTB2, DFTB3, spin, PBC/k, forces/stress, TD-DFTB, +U, SOC, … |
| xtb surface | GFN0/1/2, IPEA, analytic G; GFN-FF classical (scope-out); sTDA needs GFN GS |

### Tree vs roadmap (post harness wave-1 land on `dev`)

| Capability | Roadmap | Tree | Verdict |
| ---------- | ------- | ---- | ------- |
| AM1 / RM1 energy + ΔHf | ✓ | ship; REFERENCEs lack hard mopac E pin | energy foothold; deepen → mopac pin + closed-form F |
| NDDO nuclear forces | ✓ product path | fixed-P central FD in `NDDOSolver::compute_forces` | usable F; **not** closed-form multipole → deepen slug |
| SCC-DFTB2 H/C/N/O E | ✓ | `DFTB2` + sketch SK + **3ob-3-1** static tables | energy foothold; offline DFTB+ pin hygiene still open for mio sketch |
| DFTB analytic gradient | ✓ (restricted) | closed-form HF+Pulay+γ+rep on MioSketch closed shell | **landed** wave-1; open-shell still FD; stress still ○ |
| Spin-polarized DFTB | ✓ (collinear) | `multiplicity` / dual αβ / Mulliken m; `reg_spin-dftb2-o2` | **landed** wave-1 foothold; noncollinear blocked |
| Production 3ob organic SK | — | `DftbSkSet::ThreeOb` + `dftb_3ob_tables.h` + CH4 DFTB+ pin | **landed** wave-1 (DFTB2 SK fidelity, not DFTB3 product) |
| Non-SCC DFTB0 | ✓ | `cpu::semiempirical::DFTB0` on `dev` | **landed** 2026-07-28 |
| DFTB3 / D3 / PBC / stress / TD | ○ | absent on `dev` (tips on some `wf/*`) | **gaps** |
| GFN0-xTB | ✓ toy | sketch H/C/O | toy only → Phase C |
| GFN1 / GFN2 / IPEA / spin / PBC GFN | ○ | absent | **gaps** (late) |
| PM6 / PM7 / PM3 | ○ | no types | Phase B **oracle=mopac** |
| OM3 | ○ | no type | **blocked** — not OpenMOPAC; no primary pin |
| GFN-FF / sTDA / SE Hessian (VIB) | ○ | — | **blocked** (scope or SE-only cmake) |

### Oracle / pin hygiene (process debt — not product rows)

- REFERENCEs: wave-1 lands have `oracle:` lines; older footholds may still lack hard pins.
- `reg_h2_dftb2`, `reg_h2o_gfn0` hygiene / offline pins still deepen via pending slugs.
- SK mio tables are **not** interchangeable with DFTB.org mio-1-1 for pin identity
  (3ob path pins 3ob-3-1 under pure DFTB2 SCC).

### gaps_remaining (post harness wave-1)

**15** pending + **5** incomplete (`[~]` tips) + **9** blocked (nonzero → continue loop).
Census snapshot was wave-1 pre-land; **foothold ≠ closed** for remaining ladder.
**Never claim queue exhausted** — next GapCensus owns discovery.

## Recover — 2026-07-27 (harness wave-1)

| Source | Result |
| ------ | ------ |
| Queue greens merged | `dftb2-3ob-organic`, `dftb2-analytic-forces`, `spin-dftb2-o2` → **done** |
| Unmerged TB `wf/*` with commits ahead of `dev` | **5** → keep **`[~]`** (see `## in_progress`) |
| Isolation worktrees (TB) | none (Clearup already ran) |
| Live SE/TB specs after close | none (three wave-1 specs deleted) |

**Next:** Recover/Merge incomplete tips or GapCensus+Assign remaining pending ladder.

## Recover / harness — 2026-07-27 (harness wave-2)

| Source | Result |
| ------ | ------ |
| Queue greens merged | `dftb2-dftbplus-energy-pin`, `dftb2-periodic-gamma` → **done** |
| Unmerged TB `wf/*` with commits ahead of `dev` | **3** → keep **`[~]`** (`dftb0-h2`, `dftb3-water-mio`, `dftb2-dispersion-d3`) |
| Isolation worktrees (TB) | none (Clearup already ran) |
| Live SE/TB specs after close | none (two wave-2 specs deleted) |
| Pending | **20** (unchanged this harness; GapCensus owns refill) |
| Blocked | **9** |

**Next:** Recover remaining incomplete (rebase onto `dev`) before new pending staff;
GapCensus owns discovery — **queue not exhausted**.

## Recover note — 2026-07-27 (GapCensus wave-3)

| Source | Result |
| ------ | ------ |
| SE product on `dev` | Unchanged since harness wave-2 (no new SE merge) |
| Pending rewrite | Dropped 3 staffed force/stress slugs → `## in_progress`; added `dftb2-periodic-forces` |
| Incomplete TB `wf/*` | **6** (was 3): +threeob F, +spin F, +stress |
| Isolation worktrees (TB) | none |
| Live SE/TB specs on `dev` | none |
| Roadmap ○/✓ | prose deepen only (PBC forces); no ✓ flips |

**Next:** **Recover/rebase** the six `[~]` tips onto `dev` before staffing
remaining Phase A pending; mopac twins still empty for Phase B.

## Assign — 2026-07-27 wave 3 (parent `dev` @ `223dcd4c`)

**Priority:** incomplete `[~]` ≫ open SE specs (none) ≫ pending Phase A.
**Staff = 8** (cap). Ladder: DFTB0 / DFTB2 depth / DFTB3 energy before NDDO/GFN.

| # | Slug | Kind | Action |
| - | ---- | ---- | ------ |
| 1 | `dftb0-h2` | resume | Rebase `wf/dftb0-h2` (ahead=1 behind=45) → dual gate |
| 2 | `dftb3-water-mio` | resume | Rebase `wf/dftb3-water-mio` (ahead=1 behind=45) → dual gate |
| 3 | `dftb2-dispersion-d3` | resume | Rebase `wf/dftb2-dispersion-d3` (ahead=2 behind=45) → dual gate |
| 4 | `dftb2-analytic-forces-threeob` | resume | Rebase `wf/dftb2-analytic-forces-threeob` (ahead=1 behind=29) → dual gate |
| 5 | `dftb2-spin-analytic-forces` | resume | Rebase `wf/dftb2-spin-analytic-forces` (ahead=1 behind=29) → dual gate |
| 6 | `dftb2-analytic-stress` | resume | Rebase `wf/dftb2-analytic-stress` (ahead=1 behind=29) → dual gate |
| 7 | `dftb2-periodic-forces` | **new** | Promote from pending; lattice-image F under Γ-PBC |
| 8 | `dftb2-kpoints` | **new** | Promote from pending; k-mesh beyond Γ |

**Skipped this wave:** `dftb3-analytic-forces` (needs `dftb3-water-mio` land);
`td-dftb-h2` (ground-state depth first); all Phase B NDDO / Phase C GFN
(ladder); blocked (9).

**After staff:** incomplete `[~]` = **8**; pending = **16**; blocked = **9**.
Workers: isolation worktree on nobackup; branch `wf/<slug>`; **spec-first**
`.claude/specs/<slug>/`; rebase resumes onto current `dev` before impl;
SE-only cmake; offline DFTB+ pins only. **No push.**

## pending

### Phase A — DFTB ladder (oracle = DFTB+)

Staff these before NDDO product depth and before GFN-xTB.
**Prefer Recover of `## in_progress` over new pending.**

- [ ] `dftb3-analytic-forces` | Analytic DFTB3 nuclear forces after energy foothold | dftb3 | H2O off-eq | F | **oracle=dftb+**; after `dftb3-water-mio` lands
- [ ] `td-dftb-h2` | Linear-response TD-DFTB foothold | td-dftb | H2/H2O | excitation [eV] | **oracle=dftb+|literature**; after ground-state SCC solid

### Phase B — NDDO product depth (oracle = OpenMOPAC)

AM1/RM1 energy footholds ship. Forces are fixed-density FD — deepen with
OpenMOPAC pins and closed-form integrals. Do **not** use DFTB+ for these rows.

- [ ] `am1-water-mopac-energy` | Hard OpenMOPAC total energy / ΔHf pin for H2O AM1 | nddo/am1 | H2O | E [eV] or ΔHf | **oracle=mopac**; existing reg is finite/self-FD only; seed `tools/tb_oracle/inputs/mopac/`
- [ ] `rm1-water-mopac-energy` | Hard OpenMOPAC total energy / ΔHf pin for H2O RM1 | nddo/rm1 | H2O | E [eV] or ΔHf | **oracle=mopac**; current reg is method-discriminating only
- [ ] `nddo-closed-form-forces` | Closed-form multipole / integral-derivative NDDO nuclear forces (replace fixed-P FD) | nddo | H2O off-eq | F vs mopac + FD | **oracle=mopac**
- [ ] `pm3-organic-energy` | PM3 (or MNDO) on NDDO ladder for table pipeline | nddo | small organic | E | **oracle=mopac**; bridge MNDO→PM6; extend `SemiEmpiricalMethod`
- [ ] `pm6-water-energy` | PM6 energy vs OpenMOPAC | nddo/pm6 | H2O | E | **oracle=mopac**
- [ ] `pm6-water-forces` | PM6 forces vs OpenMOPAC / FD | nddo/pm6 | H2O off-eq | F | **oracle=mopac**
- [ ] `pm7-water-energy` | PM7 energy vs OpenMOPAC | nddo/pm7 | H2O | E | **oracle=mopac**; after PM6 foothold

### Phase C — GFN-xTB (oracle = xtb) — late only

xTB is **not** DFTB3++. Staff only when Phase A ground-state DFTB is green
or operator forces `theme=gfn|xtb`.

- [ ] `gfn0-beyond-toy` | Lift GFN0 past H/C/O sketch tables toward production GFN0 (+ REFERENCE/`oracle:` line; `h2o_gfn0` still lacks REFERENCE) | xtb | multi-element | E | **oracle=xtb|literature**
- [ ] `gfn1-xtb-organic` | GFN1-xTB energy | xtb | small organic | E [eV] | **oracle=xtb**
- [ ] `gfn2-xtb-water` | Production GFN2-xTB energy | xtb | H2O | E [eV] | **oracle=xtb**
- [ ] `gfn2-xtb-forces` | GFN2-xTB analytic forces | xtb | H2O off-eq | F | **oracle=xtb**
- [ ] `ipea-xtb-small` | IPEA-xTB energy foothold | xtb | small | E | **oracle=xtb**; after GFN2 energy
- [ ] `periodic-gfn-xtb` | Periodic GFN-xTB Γ / PBC energy | xtb | small cell | E | **oracle=xtb**
- [ ] `spin-gfn-xtb` | Spin-polarized GFN-xTB | xtb | open shell | E / spin | **oracle=xtb**; after GFN2

## in_progress

Staffed assign wave 3 (2026-07-27). Unmerged TB tips with commits ahead of
`dev` (Clearup removed worktrees; keep `[~]` until merge or drop). **Do not
demote to pending.** Resumes must **rebase** onto current `dev` before merge
(older tips behind≈45; force/stress tips behind≈29). New staff have no tip yet.

- [~] `dftb0-h2` | Non-SCC DFTB / DFTB0 energy + force skeleton | dftb0 | H2 | E, F | **oracle=dftb+**; tip `wf/dftb0-h2` @ `bcff8c4a` ahead=1 behind=45; **resume assign-w3**
- [~] `dftb3-water-mio` | SCC-DFTB3 energy (third-order) | dftb3 | H2O | E [eV] | **oracle=dftb+**; tip `wf/dftb3-water-mio` @ `f26ba40b` ahead=1 behind=45; **resume assign-w3**
- [~] `dftb2-dispersion-d3` | Dispersion-corrected DFTB (D3(BJ) water dimer) | dftb2 | water dimer | E_disp / E_tot | **oracle=dftb+|literature**; tip `wf/dftb2-dispersion-d3` @ `edabfa64` ahead=2 behind=45; **resume assign-w3**
- [~] `dftb2-analytic-forces-threeob` | Closed-form nuclear F on `DftbSkSet::ThreeOb` (replace re-SCC FD) | dftb2 | CH4/H2O | F vs FD | **oracle=dftb+|literature**; tip `wf/dftb2-analytic-forces-threeob` @ `211dc386` ahead=1 behind=29; **resume assign-w3**
- [~] `dftb2-spin-analytic-forces` | Closed-form collinear open-shell nuclear F (replace re-SCC FD) | dftb2/spin | O2 off-eq | F vs FD | **oracle=dftb+|literature**; tip `wf/dftb2-spin-analytic-forces` @ `1e8cac93` ahead=1 behind=29; **resume assign-w3**
- [~] `dftb2-analytic-stress` | Analytic cell stress | dftb2 | mol or Γ cell | stress vs FD | **oracle=dftb+|literature**; tip `wf/dftb2-analytic-stress` @ `d2da3fe9` ahead=1 behind=29; **resume assign-w3**
- [~] `dftb2-periodic-forces` | Lattice-image-consistent nuclear F under Γ-PBC (match H⁰/S/γ energy lattice sums; not min_image-only) | dftb2 | C-chain / dense cell | F vs FD | **oracle=dftb+|literature**; **new assign-w3** (no tip yet)
- [~] `dftb2-kpoints` | Periodic DFTB k-mesh beyond Γ | dftb2 | small cell | E (and optional F) | **oracle=dftb+**; **new assign-w3** (no tip yet)

## blocked

- [!] `noncollinear-spin-dftb` | blocked — needs collinear spin-DFTB first (**collinear product foothold exists; noncollinear still blocked**)
- [!] `dftb-plus-u` | blocked — needs SCC base + Hubbard-U design note
- [!] `spin-orbit-dftb` | blocked — needs specialized SOC / SK design; not ladder-head
- [!] `om3-nddo` | blocked — OM3 not in OpenMOPAC primary surface; no family-locked pin (Thiel OMx out of loop oracles)
- [!] `se-analytic-hessian` | blocked — needs `ATV_BUILD_VIB` / `cpu::vib::Hessian`; SE-only cmake keeps VIB **OFF**
- [!] `gfn-xtb-analytic-hessian` | blocked — same VIB gate; xtb semi-numeric Hessian is product-later
- [!] `gfn-ff` | blocked — classical GFN-FF, **scope-out** of electronic SE/TB loop (never pending)
- [!] `stda-xtb` | blocked — needs production GFN ground state first
- [!] `mopac-as-dftb-oracle` | **forbidden pattern** — never pin DFTB to OpenMOPAC (never pending)

## done

_(footholds + wave-1/2 product depth on `dev`)_

- [x] `dftb2-dftbplus-energy-pin` | Offline DFTB+ 24.1 H₂ total-energy pin on mio-1-1 H–H static SK (`dftb_mio_hh_sk.h`) | `reg_dftb2-dftbplus-energy-pin` | **oracle=dftb+**
- [x] `dftb2-periodic-gamma` | Γ-point periodic SCC-DFTB2 lattice-sum energy (H⁰/S, γ, repulsive) | `reg_dftb2-periodic-gamma` | **oracle=literature** (+ family dftb+); k-mesh → `dftb2-kpoints`
- [x] `dftb2-3ob-organic` | Static 3ob-3-1 organic CHNO SK (`DftbSkSet::ThreeOb`) + CH4 DFTB+ 24.1 pin | `reg_dftb2-3ob-organic` | **oracle=dftb+**
- [x] `dftb2-analytic-forces` | Closed-form DFTB2 nuclear gradient (HF+Pulay+SCC γ+repulsive); replace re-SCC FD on restricted MioSketch | `reg_dftb2-analytic-forces` | **oracle=dftb+** (formula / module boundary)
- [x] `spin-dftb2-o2` | Collinear unrestricted spin-DFTB2; O₂ triplet Mulliken m | `reg_spin-dftb2-o2` | **oracle=literature** (+ family dftb+)
- [x] `h2_dftb2` / `h2o_dftb2_mio` | SCC-DFTB2 H/C/N/O static SK sketch + forces path | absolute H–H pin closed via `dftb2-dftbplus-energy-pin`
- [x] `h2o_gfn0` | GFN0-xTB **toy** | not product GFN0/1/2
- [x] `am1_*` / `rm1_water_energy` | NDDO AM1/RM1 energy + fixed-P FD forces | further NDDO → OpenMOPAC (PM6/PM7, mopac E pin, closed-form F)
