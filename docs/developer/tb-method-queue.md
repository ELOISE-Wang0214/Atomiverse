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
| Non-SCC DFTB0 | ○ | no type on `dev` (tip on `wf/dftb0-h2`) | **gap** / incomplete tip |
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

## pending

### Phase A — DFTB ladder (oracle = DFTB+)

Staff these before NDDO product depth and before GFN-xTB.

- [ ] `dftb2-analytic-stress` | Analytic cell stress | dftb2 | mol or Γ cell | stress vs FD | **oracle=dftb+|literature**
- [ ] `td-dftb-h2` | Linear-response TD-DFTB foothold | td-dftb | H2/H2O | excitation [eV] | **oracle=dftb+|literature**; after ground-state SCC solid

### Phase B — NDDO product depth (oracle = OpenMOPAC)

AM1/RM1 energy footholds ship. Forces are fixed-density FD — deepen with
OpenMOPAC pins and closed-form integrals. Do **not** use DFTB+ for these rows.

- [ ] `am1-water-mopac-energy` | Hard OpenMOPAC total energy / ΔHf pin for H2O AM1 | nddo/am1 | H2O | E [eV] or ΔHf | **oracle=mopac**; existing reg is finite/self-FD only
- [ ] `nddo-closed-form-forces` | Closed-form multipole / integral-derivative NDDO nuclear forces (replace fixed-P FD) | nddo | H2O off-eq | F vs mopac + FD | **oracle=mopac**
- [ ] `pm3-organic-energy` | PM3 (or MNDO) on NDDO ladder for table pipeline | nddo | small organic | E | **oracle=mopac**; bridge MNDO→PM6
- [ ] `pm6-water-energy` | PM6 energy vs OpenMOPAC | nddo/pm6 | H2O | E | **oracle=mopac**
- [ ] `pm6-water-forces` | PM6 forces vs OpenMOPAC / FD | nddo/pm6 | H2O off-eq | F | **oracle=mopac**
- [ ] `pm7-water-energy` | PM7 energy vs OpenMOPAC | nddo/pm7 | H2O | E | **oracle=mopac**; after PM6 foothold

### Phase C — GFN-xTB (oracle = xtb) — late only

xTB is **not** DFTB3++. Staff only when Phase A ground-state DFTB is green
or operator forces `theme=gfn|xtb`.

- [ ] `gfn0-beyond-toy` | Lift GFN0 past H/C/O sketch tables toward production GFN0 | xtb | multi-element | E | **oracle=xtb|literature**
- [ ] `gfn1-xtb-organic` | GFN1-xTB energy | xtb | small organic | E [eV] | **oracle=xtb**
- [ ] `gfn2-xtb-water` | Production GFN2-xTB energy | xtb | H2O | E [eV] | **oracle=xtb**
- [ ] `gfn2-xtb-forces` | GFN2-xTB analytic forces | xtb | H2O off-eq | F | **oracle=xtb**
- [ ] `ipea-xtb-small` | IPEA-xTB energy foothold | xtb | small | E | **oracle=xtb**; after GFN2 energy
- [ ] `periodic-gfn-xtb` | Periodic GFN-xTB Γ / PBC energy | xtb | small cell | E | **oracle=xtb**
- [ ] `spin-gfn-xtb` | Spin-polarized GFN-xTB | xtb | open shell | E / spin | **oracle=xtb**; after GFN2

## in_progress

Unmerged TB tips with commits ahead of `dev` (Clearup removed worktrees; keep
`[~]` until merge or drop decision). **Do not demote to pending.**

- [~] `dftb0-h2` | Non-SCC DFTB / DFTB0 energy + force skeleton | dftb0 | H2 | E, F | **oracle=dftb+**; tip `wf/dftb0-h2`
- [~] `dftb2-dftbplus-energy-pin` | Offline DFTB+ energy pin on ATV DFTB2 path | dftb2 | H2 / H2O | E [eV] | **oracle=dftb+**; tip `wf/dftb2-dftbplus-energy-pin`
- [~] `dftb3-water-mio` | SCC-DFTB3 energy (third-order) | dftb3 | H2O | E [eV] | **oracle=dftb+**; tip `wf/dftb3-water-mio`
- [~] `dftb2-dispersion-d3` | Dispersion-corrected DFTB | dftb2 | water dimer | E_disp / E_tot | **oracle=dftb+|literature**; tip `wf/dftb2-dispersion-d3`
- [~] `dftb2-periodic-gamma` | Periodic DFTB Γ-point energy | dftb2 | small cell | E | **oracle=dftb+**; tip `wf/dftb2-periodic-gamma`

## blocked

- [!] `noncollinear-spin-dftb` | blocked — needs collinear spin-DFTB first (**unblocked enough for collinear product; noncollinear still blocked**)
- [!] `dftb-plus-u` | blocked — needs SCC base + Hubbard-U design note
- [!] `spin-orbit-dftb` | blocked — needs specialized SOC / SK design; not ladder-head
- [!] `om3-nddo` | blocked — OM3 not in OpenMOPAC primary surface; no family-locked pin (Thiel OMx out of loop oracles)
- [!] `se-analytic-hessian` | blocked — needs `ATV_BUILD_VIB` / `cpu::vib::Hessian`; SE-only cmake keeps VIB **OFF**
- [!] `gfn-xtb-analytic-hessian` | blocked — same VIB gate; xtb semi-numeric Hessian is product-later
- [!] `gfn-ff` | blocked — classical GFN-FF, **scope-out** of electronic SE/TB loop
- [!] `stda-xtb` | blocked — needs production GFN ground state first
- [!] `mopac-as-dftb-oracle` | **forbidden pattern** — never pin DFTB to OpenMOPAC

## done

_(footholds + wave-1 product depth on `dev`)_

- [x] `dftb2-3ob-organic` | Static 3ob-3-1 organic CHNO SK (`DftbSkSet::ThreeOb`) + CH4 DFTB+ 24.1 pin | `reg_dftb2-3ob-organic` | **oracle=dftb+**
- [x] `dftb2-analytic-forces` | Closed-form DFTB2 nuclear gradient (HF+Pulay+SCC γ+repulsive); replace re-SCC FD on restricted MioSketch | `reg_dftb2-analytic-forces` | **oracle=dftb+** (formula / module boundary)
- [x] `spin-dftb2-o2` | Collinear unrestricted spin-DFTB2; O₂ triplet Mulliken m | `reg_spin-dftb2-o2` | **oracle=literature** (+ family dftb+)
- [x] `h2_dftb2` / `h2o_dftb2_mio` | SCC-DFTB2 H/C/N/O static SK sketch + forces path | deepen via `dftb2-dftbplus-energy-pin` (**dftb+**)
- [x] `h2o_gfn0` | GFN0-xTB **toy** | not product GFN0/1/2
- [x] `am1_*` / `rm1_water_energy` | NDDO AM1/RM1 energy + fixed-P FD forces | further NDDO → OpenMOPAC (PM6/PM7, mopac E pin, closed-form F)
