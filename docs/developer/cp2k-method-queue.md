# CP2K capability queue

Durable worklist for workflow `cp2k-method-loop`. Do not delete done rows.

**Recover first:** every shift starts by finishing `## in_progress` `[~]`, open
`.claude/specs/` (GPW/PW/periodic), and unmerged `wf/<slug>` branches. Census
only rewrites **pending**. Assign priority: **incomplete ≫ pending**. Spec is
source of truth for workers.

Census rewrites **pending** from `docs/roadmap/dft.md` §10–12 (periodic / GPW /
PW / response / infra) and related MD ○ footholds that **CP2K can oracle
offline** (Quickstep GPW, plane-wave DFT, forces/stress, OT, D3, BOMD, phonons).
**done** / **blocked** / **in_progress** history is preserved (Census must not
erase incomplete).

Out of scope for this loop (use `pyscf-method-queue.md`): molecular Gaussian
post-HF, pure PySCF continuum/charge methods without CP2K twin, classical MD/FF,
CUDA-Kokkos product tracks.

Census: 2026-07-24 (parent `dev`). Source: `docs/roadmap/dft.md` §10–12 +
§3/§8 GGA/dispersion ○, `docs/roadmap/index.md` GPW/PW caveats,
`docs/roadmap/molecular-dynamics.md` §AIMD ○, `docs/roadmap/optimization.md`
variable-cell productization.
**Harness wave 1 2026-07-23:** 10 staff-batch slugs merged → `## done`.
**Assign wave 2 2026-07-23:** 6 incomplete resumes + 4 dep-safe pending →
`## in_progress` (staff 10; includes `gpw-kpoints-si`).
**Harness wave 2 2026-07-24:** 9 merged slugs → `## done` (stress, BOMD, PW
analytic stress, FD phonons, DFT+U, bands/DOS, k-UKS, PBEsol, D4). Incomplete
remaining: **`gpw-kpoints-si`** only (not in harness wave-2 staff list; no
merge/reg).
**Assign wave 3 2026-07-24:** 1 incomplete resume (`gpw-kpoints-si`) + 9
dep-safe pending → `## in_progress` (staff 10).
**Recover 2026-07-24:** incomplete_count=**10** (all wave-3 `[~]`); open
CP2K/GPW/PW specs=**0**; unmerged CP2K `wf/*` ahead of `dev`=**0**.
**Census 2026-07-24:** rewrote **pending only** from remaining roadmap ○ with
offline CP2K twins; **excluded** every slug already in `## in_progress` /
`## done` / `## blocked` (no demote of incomplete; no promote into staff).
Excluded also molecular post-HF, pure PySCF, classical MD/FF, CUDA-Kokkos,
GAPW/PAW/periodic hybrid/meta-GGA. Assign must finish incomplete first.
**Census 2026-07-25:** no CP2K/GPW/PW merges into `dev` since wave-2 harness
(only metalwalls waves + molecular scf-review chain — out of scope); all 11
prior pending re-verified as true gaps on `dev` (GPW `XcKind` still
PADE|PBE|PBEsol, GPW stress still strain-FD, no TS/MBD/vdW-DF, no
DFPT/dispersions, no periodic MP2/RPA/GW); wave-3 `[~]` ×10 preserved (no
regressions/specs/`wf/*` for any — still incomplete). Added 2 slugs from
remaining §3/§8 ○ with offline CP2K twins: `gpw-rpbe-water`, `gpw-vdw-df-ice`.
Open specs `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp` are **PySCF-queue**
deferred items (per `.claude/specs/INDEX.md`) — not this loop's. Roadmap
already matches `dev` — no ○/✓ flips.
**Assign wave (claude) 2026-07-25:** staffed 6 of the 10 wave-3 `[~]` resumes
(no pending promoted; all picks already in_progress): `gpw-kpoints-si`
(oldest resume, wave-2 leftover), `gpw-variable-cell-ice`, `dfpt-phonon-si`,
`xl-bomd-water`, `berry-polarization`, `periodic-mp2-si`. Held `[~]` for next
wave: `wannier-centers` / `cpmd-water` (surface overlap with
`berry-polarization` / `xl-bomd-water` picks — k-path and AIMD propagator
files), `mpi-gpw-scaling-4rank` / `mpi-pw-scaling-4rank` (infra; contend with
`gpw-kpoints-si` / `periodic-mp2-si` solver trees). No pick depends on
another pick. All picks are write-first: no open spec / no `wf/*` / no
regression exists for any.
**Harness wave 1 (claude) 2026-07-25:** 5 merged → done (k-points GPW,
variable-cell, DFPT Γ, Berry polarization, periodic MP2 foothold);
xl-bomd-water WIP-checkpointed; 71 merged wf/* pruned.
**Census 2026-07-25 (post harness wave 1):** no CP2K/GPW/PW merges since
`859d6e28` (HEAD); all 13 prior pending re-verified as true gaps on `dev`
(GPW `XcKind` still PADE|PBE|PBEsol; `gpw_stress` analytic fixed-P; no
TS/MBD/vdW-DF; DFPT q = 0 only; no Wannier/CPMD/RPA/GW; no PW variable-cell
or Fe-smearing regressions). Wave-3 `[~]` ×5 preserved (xl-bomd-water WIP
`61d235f4` untested; others no spec/`wf/*`/regression). Added 2 slugs for
declared-but-unqueued follow-ups: `periodic-ri-mp2-si` (foothold ≠ capability
closed — §10 Periodic MP2 stays ○ until RI/k-point lands) and
`gpw-variable-cell-molopt` (the "MOLOPT re-twin queued" claim in the
`gpw-variable-cell-ice` done row now has a real queue row). §3 BP86/B97/OLYP/
HCTH/SOGGA remain molecular-path (PySCF-loop oracle) — not added. Roadmap
already matches `dev` (harness `859d6e28` synced Berry/DFPT/k-points/
variable-cell; MP2 ○ kept) — no ○/✓ flips.
**Assign wave (claude) 2026-07-25 (wave 2):** staffed 6 dep-safe picks — all
5 wave-3 `[~]` resumes (`xl-bomd-water` = WIP `61d235f4` on
`wf/xl-bomd-water`, resume build+test first; `wannier-centers` / `cpmd-water`
now free of the prior surface holds — `berry-polarization` merged, xl-bomd
impl checkpointed; `mpi-gpw-scaling-4rank` / `mpi-pw-scaling-4rank` — the
contending `gpw-kpoints-si` / `periodic-mp2-si` solver trees merged) + 1
pending promoted: `gpw-analytic-stress` (dep `gpw-stress-water` ✓).
Skipped `second-gen-cpmd-water` (dep `cpmd-water` is itself a pick this
wave). Open specs `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp` stay
PySCF-loop deferred — not staffed. Note: `cpmd-water` and `xl-bomd-water`
may both touch AIMD propagator files — harness resolves at merge.

**Bootstrap + GapCensus 2026-07-25 (parent `dev` @ `5bdae27d`):** recover
clean — incomplete=`0`, open CP2K/GPW/PW specs=`0` (only PySCF deferred
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`), unmerged `wf/*`=`0`, isolation
worktrees=`0`, leftover `build*`=`none`. All 14 prior pending re-verified as
true gaps on tree: GPW `XcKind` still `{PADE,PBE,PBEsol}` (no BLYP/revPBE/RPBE);
no TS/MBD/vdW-DF symbols; DFPT remains Γ-only (`phonon-dispersion-si` open);
no bulk-Fe + spin+smearing reg; no PW CELL_OPT product reg; Periodic MP2
capability still ○ (Γ foothold `periodic-mp2-si` ✓ — RI/k-point + RPA/GW
pending); second-gen CP after CPMD ✓. All 30 done-slug regression dirs present
under `regressions/{gpw,pw,response,dynamics,scf}/`. Roadmap matches `dev`
(no ○/✓ flips). **Pending unchanged (14).** Ready for Assign wave — prefer
surface-disjoint staff (one XcKind GGA slug/wave; TS before MBD/vdW-DF;
`second-gen-cpmd-water` free now that `cpmd-water` ✓).
**Recover 2026-07-25 (parent `dev` @ `f6838676`):** reconfirmed clean —
incomplete_count=**0**; `## in_progress` empty (no `[~]`); open CP2K/GPW/PW
specs=**0** (PySCF deferred only); unmerged `wf/*`=**0**; isolation
worktrees=**0**; leftover builds=**none**; all 30 done-slug regression dirs
present; no half-regs for any of the 14 pending. Do **not** invent methods —
next is GapCensus/Assign from **pending**.

**GapCensus wave 1 2026-07-25 (parent `dev` @ `f6838676`):** recover clean
(incomplete=`0`, CP2K specs=`0`, `wf/*`=`0`, builds=`none`). Tree re-diff vs
`docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: all **14** prior pending still true
gaps (no half-regs; `XcKind` still `{PADE,PBE,PBEsol}`; no TS/MBD/vdW-DF;
DFPT Γ-only; Periodic MP2 capability still ○; no Fe-smearing / PW CELL_OPT
regs). **Added 4** tree-declared follow-ons with offline CP2K twins:
`gpw-uks-pbe`, `gpw-kpoints-uks`, `gpw-kpoint-forces` (from `gpw_solver.h`
scope: "UKS PBE-GGA, k-point UKS, and analytic k-point forces are follow-on"),
`periodic-wannier-si` (roadmap §11 MLWF residual after molecular Boys
foothold). **Not added:** BP86/B97/OLYP/HCTH (molecular PySCF-oracle path);
XDM/VV10 (hold behind TS/MBD/vdW-DF order); GAPW/PAW/periodic hybrid/meta-GGA
(blocked); molecular post-HF / classical MW / CUDA ports (other loops).
Roadmap: **no ○/✓ table flips** (rows already match tree); fixed stale
`docs/roadmap/index.md` GPW blurb (still said "no Gaussian k-points /
strain-FD stress" after k-points + analytic stress landed). **Pending = 18.**

**Harness wave 1 2026-07-25 (post GapCensus staff; parent `dev` @ `83a100df`):**
9 staffed slugs merged → `## done` (`gpw-blyp-water`, `gpw-uks-pbe`,
`pw-uks-fe-smearing`, `pw-variable-cell-si`, `gpw-ts-ice`,
`phonon-dispersion-si`, `periodic-wannier-si`, `second-gen-cpmd-water`,
`periodic-ri-mp2-si`). Specs deleted on close. **Pending residual = 9**
(revPBE/RPBE, k-UKS/k-forces, MBD/vdW-DF, MOLOPT re-twin, RPA/GW).
`## in_progress` empty. Roadmap flips: TS ✓, Periodic MP2 ✓ (Γ RI product;
multi-k residual), MLWF ✓ (molecular Boys + periodic PW MV), 2nd-gen CP ✓,
phonon dispersion productized (IFC Fourier; finite-q DFPT still ○), GPW BLYP
+ UKS-PBE + bulk Fe smear + PW CELL_OPT prose. **Queue not exhausted** —
GapCensus owns next discovery.

**GapCensus wave 2 2026-07-25 (parent `dev` @ `507bb5f8`):** recover —
incomplete_count was **0** in queue, but unmerged `wf/gpw-kpoint-forces`
(ahead 1: `d15a475d` feat analytic k-point forces + reg + oracle input;
merge-base `9ff2f9a8` pre-harness) → promote to **`## in_progress`**
(incomplete ≫ pending; do not leave half-landed `wf/*` as pending). Open
CP2K/GPW/PW specs on `dev` = **0** (PySCF deferred only); isolation
worktrees = **0**; leftover builds = **none**. Tree re-diff vs
`docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending (minus
k-forces) still true gaps — `XcKind` still `{PADE,PBE,PBEsol,BLYP}` (no
revPBE/RPBE); no MBD/vdW-DF symbols; DFPT remains **q = 0** only; multi-k
MP2 still throws; `gpw_stress` throws multi-k + UKS; force assembly is
closed-shell Γ RKS (UKS forces unpinned). **Added 4** tree-declared
follow-ons with offline CP2K twins: `dfpt-finite-q-si`,
`periodic-mp2-kpoints`, `gpw-kpoint-stress`, `gpw-uks-forces`. **Not
added:** BP86/B97/OLYP/HCTH (molecular PySCF path); XDM/VV10 (hold behind
MBD/vdW-DF order); GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular
post-HF / classical MW / CUDA ports (other loops). Roadmap: **no ○/✓ table
flips** (nothing new merged since harness); stale §10 intro + DFPT prose
refreshed. **Pending = 12**; **in_progress = 1**.

**Harness wave 2 2026-07-25 (parent `dev` @ `f5739506`):** 10 staffed
slugs merged → `## done` (`gpw-kpoint-forces`, `gpw-revpbe-si`,
`gpw-kpoints-uks`, `gpw-uks-forces`, `gpw-kpoint-stress`, `gpw-mbd-ice`,
`periodic-mp2-kpoints`, `periodic-rpa-si`, `dfpt-finite-q-si`,
`gpw-variable-cell-molopt`). Specs deleted on close. **Pending residual =
3** (`gpw-rpbe-water`, `gpw-vdw-df-ice`, `periodic-gw-si`).
`## in_progress` empty. Roadmap flips: revPBE ✓, free-atom MBD ✓, Periodic
RPA ✓ (Γ dRPA@KS product), finite-q DFPT productized; multi-k MP2 residual
closed; GPW k-UKS / k-forces / k-stress / Γ-UKS forces / MOLOPT CELL_OPT
prose. **Queue not exhausted** — GapCensus owns next discovery.

## pending

- [ ] `gpw-rpbe-water` | GPW RPBE GGA | gpw | water | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 RPBE ○; CP2K XC_FUNCTIONAL RPBE; shares `XcKind` plumbing — prefer after `gpw-revpbe-si` ✓ (merge-conflict avoidance)
- [ ] `gpw-vdw-df-ice` | Nonlocal vdW-DF correlation | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp / e_tot Ha | roadmap §8 vdW-DF ○; CP2K VDW_POTENTIAL NON_LOCAL (DRSLL); **after** TS ✓ / MBD ✓ order
- [ ] `periodic-gw-si` | Periodic GW | gpw / response | Si (tiny cell) | MOLOPT+GTH | QP gaps | roadmap §10 ○ Periodic GW; CP2K GW; **dep:** solid GPW + bands foothold ✓ + RPA ✓; not molecular G0W0 (PySCF loop)

## in_progress

_(empty — harness wave 2 closed all staffed `[~]`)_

## done (harness wave 2 2026-07-25)

- [x] `gpw-kpoint-forces` | GPW analytic forces at Monkhorst–Pack k-points | gpw | `regressions/gpw/gpw-kpoint-forces` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-revpbe-si` | GPW revPBE GGA (`XcKind::RevPBE`) | gpw | `regressions/gpw/gpw-revpbe-si` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-kpoints-uks` | GPW collinear multi-k UKS SCF | gpw | `regressions/gpw/gpw-kpoints-uks` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-uks-forces` | Γ-UKS analytic ionic forces (O₂ PBE) | gpw | `regressions/gpw/gpw-uks-forces` | ratio **~2.52** | 2026-07-25
- [x] `gpw-kpoint-stress` | GPW analytic cell stress at k-points | gpw | `regressions/gpw/gpw-kpoint-stress` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-mbd-ice` | Free-atom MBD@rsSCS dispersion | gpw | `regressions/gpw/gpw-mbd-ice` | ratio **perf_exempt** (libMBD E_disp) | 2026-07-25
- [x] `periodic-mp2-kpoints` | Multi-k RI-GPW MP2 residual | gpw | `regressions/gpw/periodic-mp2-kpoints` | ratio **perf_exempt** | 2026-07-25
- [x] `periodic-rpa-si` | Γ-point GPW direct RPA (dRPA@PBE) | gpw | `regressions/gpw/periodic-rpa-si` | ratio **perf_exempt** | 2026-07-25
- [x] `dfpt-finite-q-si` | Finite-q DFPT phonons at X (PW) | response | `regressions/response/dfpt-finite-q-si` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-variable-cell-molopt` | GPW CELL_OPT MOLOPT re-twin vs CP2K | gpw | `regressions/gpw/gpw-variable-cell-molopt` | ratio **~7.5** | 2026-07-25

## done (harness wave 1 2026-07-25)

- [x] `gpw-blyp-water` | GPW BLYP GGA (`XcKind::BLYP`) | gpw | `regressions/gpw/gpw-blyp-water` | ratio **perf_exempt** (no CP2K wall; dual-gate \|E_BLYP−E_PBE\|) | 2026-07-25
- [x] `gpw-uks-pbe` | Γ-UKS PBE-GGA O₂ | gpw | `regressions/gpw/gpw-uks-pbe` | ratio **0.85** | 2026-07-25
- [x] `pw-uks-fe-smearing` | bulk Fe UKS + metallic smearing | pw | `regressions/pw/pw-uks-fe-smearing` | ratio **perf_exempt** (no fair NC-UPF twin) | 2026-07-25
- [x] `pw-variable-cell-si` | PW isotropic CELL_OPT (Nielsen–Martin stress) | pw | `regressions/pw/pw-variable-cell-si` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-ts-ice` | Tkatchenko–Scheffler free-atom TS dispersion | gpw | `regressions/gpw/gpw-ts-ice` | ratio **perf_exempt** (E_disp pin; libMBD) | 2026-07-25
- [x] `phonon-dispersion-si` | Phonon dispersion ω(q) IFC Fourier (SW Si) | response | `regressions/response/phonon-dispersion-si` | ratio **perf_exempt** | 2026-07-25
- [x] `periodic-wannier-si` | Periodic multi-k MLWF (Marzari–Vanderbilt PW) | response | `regressions/pw/periodic-wannier-si` | ratio **perf_exempt** | 2026-07-25
- [x] `second-gen-cpmd-water` | Second-generation CP (ASPC + partial SCF) | dynamics | `regressions/dynamics/second-gen-cpmd-water` | ratio **perf_exempt** | 2026-07-25
- [x] `periodic-ri-mp2-si` | Γ RI-GPW MP2 product (`mp2_ri`; multi-k residual) | gpw | `regressions/gpw/periodic-ri-mp2-si` | ratio **~0.09** | 2026-07-25

## done (wave parallel land 2026-07-25)

- [x] `mpi-pw-scaling-4rank` | MPI PW k-mesh RR gate (4 ranks) | pw | `regressions/pw/mpi-pw-scaling-4rank` | ratio **perf_exempt** (speedup 2.54×; no CP2K) | 2026-07-25
- [x] `gpw-analytic-stress` | Analytic GPW fixed-P+Pulay stress | gpw | `regressions/gpw/gpw-analytic-stress` | dual-gate vs FD ~3e-10 Ha/Bohr³ | 2026-07-25
- [x] `wannier-centers` | Foster–Boys MLWF foothold (molecular RHF) | scf | `regressions/scf/wannier-centers` | ratio **perf_exempt** | 2026-07-25
- [x] `mpi-gpw-scaling-4rank` | MPI GPW k-point RR industrial gate | gpw | `regressions/gpw/mpi-gpw-scaling-4rank` | ratio **perf_exempt** (energy rank-identical; wall≈1 until grid DD) | 2026-07-25
- [x] `cpmd-water` | Car–Parrinello MD foothold (`Cpmd` + `GPWEvaluator::cpmd`) | dynamics | `regressions/dynamics/cpmd-water` | ratio **perf_exempt** (E_cons drift 2e-3) | 2026-07-25

## blocked

- [!] `periodic-hybrid-gpw` | blocked — no periodic exact exchange (same as periodic HF); HSE/PBE0 periodic ○
- [!] `gapw-dft` | blocked — GAPW not in tree; long-horizon
- [!] `ultrasoft-paw` | blocked — UPF NC + GTH only; PAW/USPP ○
- [!] `cuda-gpw` | blocked — GPW CPU-only; CUDA plane-wave / Kokkos are separate tracks (not this loop)
- [!] `periodic-meta-gga` | blocked — no kinetic-energy density (τ) channel in GPW/PW grid engines

## done

- [x] `gpw-water-pbe-mt` | Γ-RKS GPW energy (MT water) | gpw | `regressions/gpw/gpw-water-pbe-mt` | ratio **perf_exempt** (no CP2K binary at pin host) | 2026-07-23 wave 1
- [x] `gpw-forces-water` | GPW analytic forces | gpw | `regressions/gpw/gpw-forces-water` | ratio **3.56** | 2026-07-23 wave 1
- [x] `gpw-uks-o2-gamma` | Γ-point UKS GPW | gpw | `regressions/gpw/gpw-uks-o2-gamma` | ratio **0.98** | 2026-07-23 wave 1
- [x] `gpw-pbe-d3bj-ice` | GPW PBE-D3(BJ) periodic | gpw | `regressions/gpw/gpw-pbe-d3bj-ice` | ratio **perf_exempt** (E_disp pin; full QS wall unfair) | 2026-07-23 wave 1
- [x] `pw-lda-si` | PW LDA Si k-mesh energy | pw | `regressions/pw/pw-lda-si` | ratio **perf_exempt** (ABINIT-class energy band; no fair CP2K PW twin) | 2026-07-23 wave 1
- [x] `pw-pbe-gga-si` | PW PBE (true GGA) | pw | `regressions/pw/pw-pbe-gga-si` | ratio **perf_exempt** (symmetry+GGA-activation pins; no QE/CP2K binary) | 2026-07-23 wave 1
- [x] `pw-smearing-al` | metallic smearing PW | pw | `regressions/pw/pw-smearing-al` | ratio **0.008** | 2026-07-23 wave 1
- [x] `pw-spin-polarized-o2` | PW spin-polarized KS | pw | `regressions/pw/pw-spin-polarized-o2` | ratio **perf_exempt** (mag=2 + frozen e_tot pin) | 2026-07-23 wave 1
- [x] `pw-bomd-h` | PW BOMD foothold | dynamics | `regressions/pw/pw-bomd-h` | ratio **perf_exempt** (NVE drift envelope; no CP2K twin wall) | 2026-07-23 wave 1
- [x] `ot-minimizer-water` | OT vs DIIS large cell | gpw | `regressions/gpw/ot-minimizer-water` | ratio **~9.0** (fair OT-only, OMP=4) | 2026-07-23 wave 1
- [x] `gpw-stress-water` | GPW cell stress (strain-FD) | gpw | `regressions/gpw/gpw-stress-water` | ratio **4.75** | 2026-07-24 wave 2
- [x] `gpw-bomd-water` | GPW BOMD short traj | dynamics | `regressions/gpw/gpw-bomd-water` | ratio **perf_exempt** (NVE drift envelope; no CP2K AIMD twin wall) | 2026-07-24 wave 2
- [x] `pw-analytic-stress-si` | PW Nielsen–Martin analytic stress | pw | `regressions/pw/pw-analytic-stress-si` | ratio **perf_exempt** (no fair NC-UPF PW twin wall) | 2026-07-24 wave 2
- [x] `fd-phonon-supercell` | FD phonons (supercell Γ) | response | `regressions/pw/fd-phonon-supercell` | ratio **perf_exempt** (SW 64-atom Γ; no CP2K VIBRATIONAL_ANALYSIS wall) | 2026-07-24 wave 2
- [x] `gpw-dftu-nio` | Periodic DFT+U (Dudarev) | gpw | `regressions/gpw/gpw-dftu-nio` | ratio **0.21** | 2026-07-24 wave 2
- [x] `bands-dos-si` | Band structure / DOS | pw | `regressions/pw/bands-dos-si` | ratio **perf_exempt** (literature Γ→X window; no fair GTH twin wall) | 2026-07-24 wave 2
- [x] `pw-uks-kmesh-fe` | k-point UKS plane-wave | pw | `regressions/pw/pw-uks-kmesh-fe` | ratio **perf_exempt** (O₂ multi-k mag=2; no fair NC-UPF twin) | 2026-07-24 wave 2
- [x] `gpw-pbesol-si` | GPW PBEsol GGA | gpw | `regressions/gpw/gpw-pbesol-si` | ratio **perf_exempt** (dual-gate PBE vs PBEsol; no CP2K binary) | 2026-07-24 wave 2
- [x] `gpw-d4-water` | GPW DFT-D4 periodic | gpw | `regressions/gpw/gpw-d4-water` | ratio **perf_exempt** (dftd4 3.7.0 E_disp pin; CP2K Arrhenius build lacks DFTD4) | 2026-07-24 wave 2
- [x] `gpw-kpoints-si` | Monkhorst–Pack k-point GPW RKS (phased lattice sums, KB nonlocal, per-k DIIS; diamond-Si vs CP2K Quickstep) | gpw | `regressions/gpw/gpw-kpoints-si` | ratio **perf_exempt** (no CP2K binary at pin host) | 2026-07-25 claude wave 1
- [x] `gpw-variable-cell-ice` | GPW variable-cell relax (evaluator re-target, stress axes mask, `StressConverged`; a\* = 3.2111 Å vs CP2K CELL_OPT 2026.1) | gpw | `regressions/gpw/gpw-variable-cell-ice` | walls ATV 1.9 s / CP2K 3.9 s (GTH chain; MOLOPT re-twin queued) | 2026-07-25 claude wave 1
- [x] `dfpt-phonon-si` | Γ-point DFPT phonons for PW (Sternheimer PCG, Ewald Hessian, LDA kernel; energy-FD twin 0.07 %) | response | `regressions/response/dfpt-phonon-si` | ratio **perf_exempt** (≈24 s wall; no CP2K DFPT wall recorded) | 2026-07-25 claude wave 1
- [x] `berry-polarization` | KSV Berry-phase polarization on PW (+ Hermitian Davidson/g-vector fixes; symmetry-exact Si pins) | response | `regressions/pw/berry-polarization` | ratio **perf_exempt** (no CP2K binary at pin host; ATV ~18 s) | 2026-07-25 claude wave 1
- [x] `periodic-mp2-si` | Γ-point canonical GPW MP2 **foothold** (`GpwDft::mp2`, Poisson-per-pair (ia\|jb); sc-Si2) — capability row stays ○, RI/k-point + RPA/GW remain | gpw | `regressions/gpw/periodic-mp2-si` | ratio **perf_exempt** (no CP2K binary at pin host) | 2026-07-25 claude wave 1
- [x] `xl-bomd-water` | Niklasson dissipative XL-BOMD on GPW RKS (`XlBomd` K∈[3,9], `GpwDft::run(guess_P)`, `GPWEvaluator::xl_bomd`); H₂O NVE drift envelope | dynamics | `regressions/dynamics/xl-bomd-water` | ratio **perf_exempt** (no CP2K binary; drift 1.2e-5 / 33 s) | 2026-07-25
