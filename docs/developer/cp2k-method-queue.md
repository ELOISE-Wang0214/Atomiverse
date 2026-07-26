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

**GapCensus wave 3 2026-07-25 (parent `dev` @ `def040bd`):** recover
clean — incomplete_count=**0**; `## in_progress` empty; open CP2K/GPW/PW
specs=**0** (PySCF deferred only: `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`);
unmerged `wf/*`=**0**; isolation worktrees=**0**; leftover builds=**none**.
Tree re-diff vs `docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD +
`tools/cp2k_oracle` + `include/atomiverse/cpu/{gpw,planewave}`: prior
pending **3** still true gaps — `XcKind` still
`{PADE,PBE,PBEsol,RevPBE,BLYP}` (no RPBE); no vdW-DF / NON_LOCAL symbols;
no periodic GW. **Added 5** tree-declared / roadmap residuals with offline
CP2K twins: `gpw-kpoint-uks-forces` (forces throw “k-point UKS not
implemented”; header follow-on after Γ-UKS + k-RKS forces ✓),
`gpw-uks-stress` (analytic stress throws UKS at Γ and multi-k),
`periodic-rpa-kpoints` (multi-k RPA residual — `direct_rpa` throws nk>1;
Γ dRPA foothold ≠ multi-k capability), `gpw-mbd-hirshfeld` (roadmap §8
Hirshfeld volumes / density-scaled MBD after free-atom MBD ✓),
`dfpt-incommensurate-si` (finite-q product is commensurate supercell;
incommensurate Sternheimer still header follow-on). **Not added:**
BP86/B97/OLYP/HCTH (molecular PySCF path); XDM/VV10/vdW-DF2 (hold behind
vdW-DF order); GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular
post-HF / classical MW / CUDA ports (other loops). Roadmap: **no ○/✓
table flips** (nothing new merged); open-follow-on + forces/stress/RPA
prose refreshed; stale index Cell “GPW strain-FD” → analytic. **Pending
= 8**; **in_progress = 0**.

**Harness wave 3 2026-07-25 (parent `dev` @ `0a642122`):** 8 staffed
slugs merged → `## done` (`gpw-rpbe-water`, `gpw-kpoint-uks-forces`,
`gpw-uks-stress`, `gpw-vdw-df-ice`, `gpw-mbd-hirshfeld`,
`periodic-rpa-kpoints`, `periodic-gw-si`, `dfpt-incommensurate-si`).
Specs deleted on close. **Pending residual = 0** (staff batch emptied
the census pending list — **not** queue exhausted; GapCensus owns next
discovery). `## in_progress` empty. Roadmap flips: RPBE ✓, vdW-DF ✓
(DRSLL), Periodic GW ✓ (Γ spectral G0W0@PBE; multi-k residual),
Hirshfeld-scaled MBD/TS prose; multi-k dRPA residual closed; Γ-UKS
stress productized (multi-k UKS stress residual); multi-k UKS forces
productized; incommensurate DFPT productized. **Blocked = 5** unchanged.

**GapCensus wave 4 2026-07-26 (parent `dev` @ `decf917f`):** recover
clean — incomplete_count=**0**; `## in_progress` empty; open CP2K/GPW/PW
specs=**0** (PySCF deferred only: `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`);
unmerged `wf/*`=**0**; isolation worktrees=**0**; leftover builds=**none**.
Tree re-diff vs `docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD +
`tools/cp2k_oracle` + `include/atomiverse/cpu/{gpw,planewave}`: prior
pending residual **0** (harness emptied list — not capability-closed).
**True residuals re-opened (9)** with offline CP2K / literature twins:
`gpw-kpoint-uks-stress` (throws multi-k UKS; after Γ-UKS + k-RKS stress ✓),
`periodic-gw-kpoints` (Γ G0W0 foothold ≠ multi-k; `g0w0` throws nk>1),
`gpw-vdw-df2-ice` / `gpw-rvv10-ice` (§8 vdW-DF2 / rVV10 ○ after DRSLL ✓),
`gpw-mbd-periodic` (free-atom MBD ✓; periodic Ewald throw),
`pw-spin-gga` (PW `step_spin` σ=0 LSDA-limit; true collinear spin-GGA
deferred), `gpw-uks-blyp` (UKS spin-GGA wrappers only PBE/PBEsol; RKS
BLYP ✓), `periodic-rpa-finite-q` (k-diagonal multi-k dRPA ✓; full
q-resolved dielectric follow-on), `gpw-dftu-stress` (PLUS_U energy ✓;
stress throws DFT+U). **Not added:** BP86/B97/OLYP/HCTH (molecular
PySCF path); XDM (hold behind rVV10/vdW-DF2 product order — weaker stock
CP2K twin); shear/off-diagonal stress; OT-UKS; CPMD RATTLE; GPW metallic
smearing; UMP2/UGW/URPA; GAPW/PAW/periodic hybrid/meta-GGA (blocked);
molecular post-HF / classical MW / CUDA ports (other loops). Roadmap:
**no ○/✓ table flips** (nothing new merged since harness wave 3); §10/§8
open-residual prose + index Cell/GPW residual line refreshed to match
queue. **Pending = 9**; **in_progress = 0**; **blocked = 5**.

**Harness wave 4 2026-07-26 (parent `dev` @ `85eaa2d9`):** 9 staffed
slugs merged → `## done` (`gpw-kpoint-uks-stress`, `periodic-gw-kpoints`,
`gpw-vdw-df2-ice`, `gpw-rvv10-ice`, `gpw-mbd-periodic`, `pw-spin-gga`,
`gpw-uks-blyp`, `periodic-rpa-finite-q`, `gpw-dftu-stress`). Specs
deleted on close. **Pending residual = 0** (staff batch emptied the
census pending list — **not** queue exhausted; GapCensus owns next
discovery). `## in_progress` empty. Roadmap flips: vdW-DF2 ✓ (LMKLL),
rVV10 ✓; multi-k G0W0 residual closed; periodic Ewald MBD productized;
PW true collinear spin-GGA productized; Γ-UKS BLYP; multi-k UKS stress;
finite-q inter-k dRPA; DFT+U analytic stress. **Blocked = 5** unchanged.

**GapCensus wave 5 2026-07-26 (parent `dev` @ `4437d80f`):** recover
clean — incomplete_count=**0**; `## in_progress` empty; open CP2K/GPW/PW
specs=**0** (PySCF deferred only: `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`);
unmerged `wf/*`=**0**; isolation worktrees=**0**; leftover builds=**none**.
Tree re-diff vs `docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD +
`tools/cp2k_oracle` + `include/atomiverse/cpu/{gpw,planewave}`: prior
pending residual **0** (harness emptied list — not capability-closed).
**True residuals re-opened (9)** with offline CP2K / literature twins:
`gpw-dftu-kpoints` (`build` rejects k-mesh + `plus_u`; Γ DFT+U energy +
stress ✓ ≠ multi-k PLUS_U), `periodic-gw-finite-q` (multi-k G0W0 is
**k-diagonal** only; full inter-\(k\) / finite-\(q\) \(W\) header
follow-on after `periodic-gw-kpoints` ✓ / `periodic-rpa-finite-q` ✓),
`gpw-uks-revpbe` (UKS XcKind still PADE|PBE|PBEsol|BLYP; RKS RevPBE ✓),
`gpw-xdm-ice` (§8 XDM ○ — order hold behind rVV10/vdW-DF2 now cleared),
`gpw-vv10-ice` (stock VV10 ○, distinct from product rVV10),
`cpmd-rattle-water` (CPMD SHAKE-free foothold; RATTLE follow-on),
`gpw-metallic-smearing` (GPW fixed occupations only; PW metallic smear ✓),
`periodic-ump2-si` (`mp2` throws closed-shell; UMP2 follow-on),
`ot-uks-o2` (OT + UKS throws). **Not added:** BP86/B97/OLYP/HCTH
(molecular PySCF path); `gpw-uks-rpbe` (thin after revPBE UKS); shear /
off-diagonal stress; optB88-vdW (needs exchange + nonlocal combo);
UGW/URPA (hold behind UMP2); Hirshfeld×k-points; DFPT-GGA (DFPT LSDA
only); XL-BOMD UKS/k; nonlocal vdW SCF potential / forces; GAPW/PAW/
periodic hybrid/meta-GGA (blocked); molecular post-HF / classical MW /
CUDA ports (other loops). Roadmap: **no ○/✓ table flips** (nothing new
merged since harness wave 4); §10 open-residual prose + index Cell
stale “multi-k UKS residual” refreshed. **Pending = 9**; **in_progress
= 0**; **blocked = 5**.

**Harness wave 5 2026-07-26 (parent `dev` @ `2d43d1c2`):** 6 staffed
slugs merged → `## done` (`gpw-dftu-kpoints`, `periodic-gw-finite-q`,
`gpw-uks-revpbe`, `cpmd-rattle-water`, `periodic-ump2-si`, `ot-uks-o2`).
Specs deleted on close. **Pending residual = 1** (`gpw-xdm-ice`);
**in_progress = 2** unmerged staff `[~]` (`gpw-metallic-smearing`
ahead on `wf/gpw-metallic-smearing`; `gpw-vv10-ice` docs/block note on
`wf/gpw-vv10-ice` — CP2K stock VV10 oracle gap). `## incomplete` kept
as `[~]` (do not demote). Roadmap residual prose: multi-k DFT+U,
finite-q G0W0 \(W\), Γ-UKS revPBE, CPMD electronic RATTLE, UMP2,
OT-UKS productized; **no new ○→✓ table flips** (rows already ✓ or
still-open XDM/VV10). **Blocked = 5** unchanged. **Queue not
exhausted** — GapCensus owns next discovery.

**GapCensus wave 6 2026-07-26 (parent `dev` @ `5847a35a`):** recover —
incomplete_count=**2** (`[~] gpw-metallic-smearing` unmerged
`wf/gpw-metallic-smearing` @ `31b00478` feat FD/MP multi-k RKS + reg;
`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @ `0197163f` docs/block —
CP2K 2026.1 `TYPE VV10` enum reject, only DRSLL|LMKLL|RVV10; do not
demote). Open CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred only:
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**. Tree re-diff vs `docs/roadmap/dft.md` §3/§8/
§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending **1**
(`gpw-xdm-ice`) still true (no XDM symbols; §8 ○). **Added 8**
tree-declared follow-ons with offline CP2K / literature twins:
`gpw-uks-rpbe` (Γ/k UKS XcKind still PADE|PBE|PBEsol|RevPBE|BLYP — no
RPBE spin path / no `xc_grid_*_rpbe_spin`; RKS RPBE ✓),
`gpw-dftu-kpoint-stress` (multi-k DFT+U SCF ✓; `gpw_stress` k-path
throws `plus_u`), `periodic-ump2-kpoints` (Γ UMP2 ✓; multi-k UMP2
throws), `periodic-ugw-si` / `periodic-urpa-si` (closed-shell G0W0/dRPA
✓; UKS throws “UGW/URPA follow-on” after UMP2 foothold),
`gpw-hirshfeld-kpoints` (Hirshfeld TS/MBD Γ RKS ✓; k-points + UKS
throw), `dfpt-gga-si` (PW DFPT LSDA-only gate; GGA kernel residual),
`xl-bomd-uks` (XL-BOMD closed-shell RKS Γ only; UKS/`use_ot` throw).
**Not added:** BP86/B97/OLYP/HCTH (molecular PySCF path); shear /
off-diagonal stress (orthorhombic design); optB88-vdW (exchange +
nonlocal combo after XDM/rVV10 product depth); OT multi-k / OT+DFT+U
(build rejects); second-gen CP UKS; nonlocal vdW SCF potential /
forces (post-SCF E_disp product); multi-k UKS stress GGA expand
(RevPBE/BLYP — thin after PADE/PBE/PBEsol k-UKS stress ✓);
GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular post-HF /
classical MW / CUDA ports (other loops). Roadmap: **no ○/✓ table
flips** (nothing new merged since harness wave 5); §3/§8/§10–12 open-
residual prose + index GPW residual line refreshed for wave-6 queue.
**Pending = 9**; **in_progress = 2**; **blocked = 5**.

**Harness wave 6 2026-07-26 (parent `dev` @ `0869833a`):** 9 staffed
slugs merged → `## done` (`gpw-metallic-smearing`, `gpw-xdm-ice`,
`gpw-uks-rpbe`, `gpw-dftu-kpoint-stress`, `periodic-ump2-kpoints`,
`periodic-ugw-si`, `gpw-hirshfeld-kpoints`, `dfpt-gga-si`, `xl-bomd-uks`).
Specs deleted on close. **Pending residual = 1** (`periodic-urpa-si`);
**in_progress = 1** unmerged staff `[~]` (`gpw-vv10-ice` on
`wf/gpw-vv10-ice` — CP2K stock VV10 oracle gap; do not demote).
Roadmap: **XDM ✓** table flip; residual prose closed GPW metallic smear,
Γ-UKS RPBE, multi-k DFT+U stress, multi-k UMP2, UGW, multi-k Hirshfeld,
DFPT-GGA, XL-BOMD UKS; stock VV10 / URPA stay open. **Blocked = 5**
unchanged. **Queue not exhausted** — GapCensus owns next discovery.

**GapCensus wave 7 2026-07-26 (parent `dev` @ `9bdeddf1`):** recover —
incomplete_count=**1** (`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @
`9be3808a` docs/block — CP2K 2026.1 has no `TYPE VV10`; do not demote).
Open CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred only:
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**. Tree re-diff vs `docs/roadmap/dft.md` §3/§8/
§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending **1**
(`periodic-urpa-si`) still true (`direct_rpa` throws UKS — “URPA is a
follow-on”; closed-shell multi-k dRPA + Γ UGW ✓). **Added 8**
tree-declared follow-ons with offline CP2K / literature twins:
`periodic-ugw-kpoints` (Γ UGW ✓; `g0w0` throws multi-k UKS — “multi-k
UGW is a follow-on”), `gpw-uks-smearing` (multi-k RKS FD/MP ✓; build
rejects UKS+smear), `gpw-hirshfeld-uks` (Hirshfeld TS/MBD Γ+multi-k RKS
✓; UKS throws), `cpmd-uks` (CPMD + RATTLE RKS Γ ✓; evaluator rejects
UKS), `second-gen-cpmd-uks` (SGCP ASPC+partial SCF RKS ✓; UKS throw),
`ot-kpoints` (OT Γ RKS+UKS ✓; k-mesh rejects OT), `ot-dftu` (OT and
DFT+U each product; OT+PLUS_U throws on UKS OT path), `dfpt-spin-si`
(DFPT Γ/finite-q/incommensurate/GGA ✓; collinear spin unsupported).
**Not added:** BP86/B97/OLYP/HCTH (molecular PySCF path); stock VV10
(kept `[~]` incomplete — no demote; no CP2K twin); optB88-vdW (exchange
+ nonlocal combo); shear/off-diagonal stress (orthorhombic design);
XL-BOMD multi-k / OT-XL (hold behind OT-k + XL depth); nonlocal vdW SCF
potential / forces (post-SCF E_disp product depth); DFPT metallic smear
(hold behind spin DFPT); GAPW/PAW/periodic hybrid/meta-GGA (blocked);
molecular post-HF / classical MW / CUDA ports (other loops). Roadmap:
**no ○/✓ table flips** (nothing new merged since harness wave 6); §8/§10
open-residual prose + index GPW residual line refreshed for wave-7
queue. **Pending = 9**; **in_progress = 1**; **blocked = 5**.

**Harness wave 7 2026-07-26 (parent `dev` @ `8cb47664`):** 9 staffed
slugs merged → `## done` (`periodic-urpa-si`, `periodic-ugw-kpoints`,
`gpw-uks-smearing`, `gpw-hirshfeld-uks`, `cpmd-uks`, `second-gen-cpmd-uks`,
`ot-kpoints`, `ot-dftu`, `dfpt-spin-si`). Specs deleted on close.
**Pending residual = 0** (staff batch emptied the census pending list —
**not** queue exhausted; GapCensus owns next discovery).
**in_progress = 1** unmerged staff `[~]` (`gpw-vv10-ice` on
`wf/gpw-vv10-ice` — CP2K stock VV10 oracle gap; do not demote).
Roadmap: **no new ○→✓ table flips** (Periodic RPA/GW/DFPT/OT/smearing/
CPMD/SGCP rows already ✓); residual prose closed URPA@UKS, multi-k UGW,
GPW UKS smear, UKS Hirshfeld, CPMD/SGCP UKS, OT+k, OT+DFT+U, collinear-
spin DFPT; stock VV10 stays open. **Blocked = 5** unchanged. **Queue not
exhausted** — GapCensus owns next discovery.

**GapCensus wave 8 2026-07-26 (parent `dev` @ `9c57a290`):** recover —
incomplete_count=**1** (`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @
`d15cd804` docs/block — CP2K 2026.1 has no `TYPE VV10`; do not demote).
Open CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred only:
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**. Tree re-diff vs `docs/roadmap/dft.md` §3/§8/
§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending residual **0**
(harness emptied list — not capability-closed). **True residuals
re-opened (9)** with offline CP2K / literature twins:
`periodic-urpa-kpoints` (`direct_rpa` throws multi-k UKS — “multi-k URPA
is a follow-on”; Γ URPA + multi-k closed-shell dRPA ✓),
`ot-smearing` (OT and metallic smear each product; `build` throws
OT+smear — fixed integer OT occupations),
`ot-uks-kpoints` (Γ-UKS OT ✓ + multi-k RKS OT ✓; multi-k UKS OT reject),
`ot-dftu-kpoints` (Γ OT+DFT+U ✓ + multi-k OT / multi-k DFT+U DIIS ✓;
multi-k OT+PLUS_U not wired),
`gpw-d4-stress` (D4 energy ✓ `reg_gpw-d4-water`; analytic
`gpw_stress` throws D4 on Γ and multi-k),
`gpw-vdw-scf-potential` (vdW-DF/DF2/rVV10 post-SCF \(E_{\mathrm{disp}}\)
✓; nonlocal potential back-propagation / forces follow-on per
`vdw_df.h`),
`dfpt-metallic-si` (DFPT Γ/finite-q/incommensurate/GGA/spin ✓; metallic
smearing unsupported on Sternheimer path — hold cleared after spin DFPT),
`dfpt-spin-incommensurate` (collinear spin DFPT ✓ at Γ/commensurate;
incommensurate-\(q\) + spin throws),
`gpw-optb88-vdw` (§8 optB88-vdW ○ — exchange + nonlocal combo after
vdW-DF/rVV10 product depth). **Not added:** BP86/B97/OLYP/HCTH
(molecular PySCF path); stock VV10 (kept `[~]` incomplete — no demote;
no CP2K twin); shear/off-diagonal stress (orthorhombic design);
D4/MBD analytic forces (hold behind D4 stress / SCF-vdW order);
XL-BOMD multi-k / OT-XL / SGCP-OT corrector (hold behind OT+smear /
OT-k-UKS depth this wave); UKS multigrid product collocation (infra
multigrid ✓ single-grid product; O(mHa) residual vs CP2K NGRIDS);
spin-polarized / metallic DOS foothold; UGW@UKS-RPBE XC throw (thin
after UKS RPBE SCF); GAPW/PAW/periodic hybrid/meta-GGA (blocked);
molecular post-HF / classical MW / CUDA ports (other loops). Roadmap:
**no ○/✓ table flips** (nothing new merged since harness wave 7); §8/
§10–12 open-residual prose + index GPW residual line refreshed for
wave-8 queue. **Pending = 9**; **in_progress = 1**; **blocked = 5**.

**Harness wave 8 2026-07-26 (parent `dev` @ `a6752916`):** 9 staffed
slugs merged → `## done` (`periodic-urpa-kpoints`, `ot-smearing`,
`ot-uks-kpoints`, `ot-dftu-kpoints`, `gpw-d4-stress`,
`gpw-vdw-scf-potential`, `dfpt-metallic-si`, `dfpt-spin-incommensurate`,
`gpw-optb88-vdw`). Specs deleted on close. **Pending residual = 0**
(staff batch emptied the census pending list — **not** queue exhausted;
GapCensus owns next discovery). **in_progress = 1** unmerged staff
`[~]` (`gpw-vv10-ice` on `wf/gpw-vv10-ice` — CP2K stock VV10 oracle gap;
do not demote). Roadmap: **optB88-vdW ✓**; residual prose closed multi-k
URPA, OT+smear, multi-k UKS OT, multi-k OT+DFT+U, D4 analytic stress,
nonlocal vdW SCF potential/forces, metallic DFPT, collinear-spin
incommensurate DFPT; stock VV10 stays open. **Blocked = 5** unchanged.
**Queue not exhausted** — GapCensus owns next discovery.

**GapCensus wave 9 2026-07-26 (parent `dev` @ `07dd7bf3`):** recover —
incomplete_count=**1** (`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @
`b8a2e6ca` docs/block — CP2K 2026.1 has no `TYPE VV10`; do not demote).
Open CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred only:
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**. Tree re-diff vs `docs/roadmap/dft.md` §3/§8/
§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending residual **0**
(harness emptied list — not capability-closed). **True residuals
re-opened (9)** with offline CP2K / literature twins:
`ot-uks-smearing` (RKS OT+smear ✓ + multi-k UKS OT ✓; `build` /
k-UKS OT path throws UKS OT+smear — integer occupations only),
`ot-uks-dftu` (Γ OT+DFT+U ✓ + multi-k RKS OT+DFT+U ✓ + multi-k UKS OT ✓;
multi-k UKS OT + PLUS_U throws “not wired”),
`ot-dftu-smearing` (OT+smear ✓ + multi-k OT+DFT+U ✓; multi-k OT path
throws smear+PLUS_U together),
`gpw-d4-forces` (D4 energy ✓ + analytic stress ✓; Γ forces skip D4 —
energy-only policy; `dftd4_energy` has no force twin),
`gpw-mbd-forces` (MBD energy ✓ + Γ/periodic stress ✓; forces energy-only
v1 per `mbd_dispersion.h`),
`gpw-uks-multigrid` (RKS multigrid product collocation ✓ when `ngrids>1`;
UKS throws “multigrid product collocation is not wired”),
`dfpt-spin-smearing` (spin DFPT ✓ + metallic DFPT ✓ separately; collinear
spin + smear rejected on Γ Sternheimer path),
`dfpt-finite-q-smearing` (finite-q DFPT ✓ fixed-occ + Γ metallic DFPT ✓;
finite-q + smear throws — “use Γ restricted smear”),
`second-gen-cpmd-ot` (SGCP ASPC+partial SCF RKS/UKS ✓; OT corrector
throws — diagonalization SCF required). **Not added:** BP86/B97/OLYP/
HCTH (molecular PySCF path); stock VV10 (kept `[~]` incomplete — no
demote; no CP2K twin); shear/off-diagonal stress (orthorhombic design);
multi-k UKS stress XcKind expand RevPBE/RPBE/BLYP (thin after
PADE/PBE/PBEsol k-UKS stress ✓); UGW@UKS-RPBE XC throw (thin after UKS
RPBE SCF); spin/metallic bands-DOS footholds (hold behind DFPT
spin+smear depth); XL-BOMD multi-k / OT-XL (hold behind OT combo depth
this wave); multi-k dispersion forces on k-path (hold behind D4/MBD
force product); GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular
post-HF / classical MW / CUDA ports (other loops). Roadmap: **no ○/✓
table flips** (nothing new merged since harness wave 8); §8/§10–12 +
MD AIMD open-residual prose + index GPW residual line refreshed for
wave-9 queue. **Pending = 9**; **in_progress = 1**; **blocked = 5**.

**Harness wave 9 2026-07-26 (parent `dev` @ `dab2b190`):** 8 staffed
slugs merged → `## done` (`ot-uks-smearing`, `ot-uks-dftu`,
`ot-dftu-smearing`, `gpw-d4-forces`, `gpw-mbd-forces`,
`gpw-uks-multigrid`, `dfpt-finite-q-smearing`, `second-gen-cpmd-ot`).
Specs deleted on close. **Pending residual = 0** (staff batch emptied
the landed list — **not** queue exhausted; GapCensus owns next
discovery). **in_progress = 2** unmerged staff `[~]`: `gpw-vv10-ice` on
`wf/gpw-vv10-ice` (CP2K stock VV10 oracle gap; do not demote) +
`dfpt-spin-smearing` on `wf/dfpt-spin-smearing` @ `d7576dea` (feat + reg
on branch; not merged to `dev` — keep incomplete). Roadmap: **no new
○→✓ table flips** (D4/MBD/OT/multigrid rows already ✓; stock VV10 stays
○); residual prose closed UKS OT+smear, multi-k UKS OT+DFT+U,
OT+DFT+U+smear, D4/MBD analytic forces, UKS multigrid, finite-q DFPT
smear, SGCP OT corrector; **DFPT spin+smear** stays open incomplete;
stock VV10 stays open. **Blocked = 5** unchanged. **Queue not
exhausted** — GapCensus owns next discovery.

**GapCensus wave 10 2026-07-26 (parent `dev` @ `94c2401d`):** recover —
incomplete_count=**2** (`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @
`2f1d90e2` docs/block — CP2K 2026.1 has no `TYPE VV10`; do not demote;
`[~] dfpt-spin-smearing` unmerged `wf/dfpt-spin-smearing` @ `d7576dea`
feat+reg — not merged; do not demote). Open CP2K/GPW/PW specs on
`dev`=**0** (PySCF deferred only: `periodic-rhf-gamma` /
`rsh-eri-cam-b3lyp`); isolation worktrees=**0**; leftover builds=**none**.
Tree re-diff vs `docs/roadmap/dft.md` §3/§8/§10–12 + MD AIMD +
`tools/cp2k_oracle` + `include/atomiverse/cpu/{gpw,planewave}`: prior
pending residual **0** (harness emptied list — not capability-closed).
**True residuals re-opened (9)** with offline CP2K / literature twins:
`ot-uks-dftu-smearing` (multi-k UKS OT + DFT+U + metallic smear —
`build` / k-UKS OT path reject; RKS OT+U+smear + UKS OT+smear + UKS
OT+U product),
`second-gen-cpmd-uks-ot` (SGCP RKS OT corrector ✓; UKS OT + ASPC
throws “follow-on”),
`xl-bomd-kpoints` (Γ RKS/UKS XL ✓; multi-k XL product depth / pin after
wave-9 OT depth),
`bands-spin-si` (`Bands::compute` throws spin-polarized foothold),
`dos-metallic-si` (`Dos::compute` throws metallic smearing foothold),
`periodic-wannier-uks` (`PeriodicWannier` throws spin-polarized MLWF),
`berry-spin-si` (`BerryPolarization` throws spin-polarized strings),
`gpw-mbd-kpoint-forces` (Γ MBD ionic forces ✓; multi-k force path folds
d3/d4/ts/xdm but **skips MBD**),
`ugw-uks-rpbe` (`g0w0` UKS RPBE spin XC throws — not product; PBE/revPBE
UKS UGW ✓). **Not added:** BP86/B97/OLYP/HCTH (molecular PySCF path);
stock VV10 (kept `[~]` incomplete — no demote; no CP2K twin); shear /
off-diagonal stress (orthorhombic design); electronic CPMD thermostat
(hold behind SGCP UKS OT / XL×k depth); TS/XDM dedicated force regs
(Γ force assembly already folds ts/xdm — thin after MBD×k); spin DOS /
metallic bands thin siblings (hold behind bands-spin + dos-metallic
first); finite-q DFPT spin+smear (hold behind incomplete Γ spin+smear
merge); XL+OT architectural exclusion (density vs orbital parameter);
GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular post-HF /
classical MW / CUDA ports (other loops). Roadmap: **no ○/✓ table
flips** (nothing new merged since harness wave 9); §8/§10–12 + MD AIMD
open-residual prose + index GPW residual line refreshed for wave-10
queue. **Pending = 9**; **in_progress = 2**; **blocked = 5**.

**Harness wave 10 2026-07-26 (parent `dev` @ `06ecaa66`):** 9 staffed
slugs merged → `## done` (`dfpt-spin-smearing`, `ot-uks-dftu-smearing`,
`second-gen-cpmd-uks-ot`, `xl-bomd-kpoints`, `bands-spin-si`,
`dos-metallic-si`, `periodic-wannier-uks`, `berry-spin-si`,
`gpw-mbd-kpoint-forces`). Specs deleted on close. **Pending residual =
1** (`ugw-uks-rpbe` — not staffed this wave; **not** queue exhausted;
GapCensus owns next discovery). **in_progress = 1** unmerged staff
`[~]`: `gpw-vv10-ice` on `wf/gpw-vv10-ice` (CP2K stock VV10 oracle gap;
do not demote). `dfpt-spin-smearing` incomplete resolved by merge into
`dev`. Roadmap: **no new ○→✓ table flips** (DFPT / Berry / Wannier /
bands-DOS / MBD / OT / SGCP / XL rows already ✓; stock VV10 stays ○);
residual prose closed collinear DFPT spin+smear, multi-k UKS OT+U+smear,
SGCP UKS OT, multi-k XL-BOMD, spin bands, metallic DOS, spin MLWF,
spin Berry, multi-k MBD forces; **UGW@UKS-RPBE** + stock VV10 stay open.
**Blocked = 5** unchanged. **Queue not exhausted** — GapCensus owns next
discovery.

**GapCensus wave 11 2026-07-26 (parent `dev` @ `96a0a7d2`):** recover —
incomplete_count=**1** (`[~] gpw-vv10-ice` unmerged `wf/gpw-vv10-ice` @
`a60c73d5` docs/block — CP2K 2026.1 `TYPE VV10` enum reject, only
DRSLL|LMKLL|RVV10; **resolved → blocked**, no force-merge / no product
`GpwDft::vv10`). Open CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred
only: `periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**. Tree re-diff vs `docs/roadmap/dft.md` §3/§8/
§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending residual **1**
(`ugw-uks-rpbe`) still true (`g0w0` `eval_xc_grid_spin` throws UKS RPBE
despite SCF `xc_grid_rpbe_spin` product). **Added 8** tree-declared
follow-ons with offline CP2K / literature twins:
`xl-bomd-uks-kpoints` (multi-k RKS XL ✓ + Γ-UKS XL ✓; multi-k UKS XL
throws dual-spin `P_k^σ`),
`bands-metallic-si` (`Bands::compute` throws metallic smearing foothold
after `dos-metallic-si` ✓),
`dos-spin-si` (`Dos::compute` throws spin-polarized foothold after
`bands-spin-si` ✓),
`dfpt-finite-q-spin-smearing` (Γ spin+smear ✓ + finite-q smear RKS ✓ +
finite-q spin fixed-occ ✓; finite-q collinear spin+smear throws),
`gpw-kpoint-uks-stress-gga` (multi-k UKS stress product for PADE/PBE/
PBEsol only; Γ UKS stress already RevPBE/RPBE/BLYP — expand k-path
XcKind),
`cpmd-electronic-thermostat` (CPMD + RATTLE + UKS/SGCP depth ✓;
electronic thermostat residual per MD AIMD roadmap),
`gpw-ts-forces` / `gpw-xdm-forces` (force assembly already folds TS/XDM
on Γ + multi-k; dedicated offline pins after MBD×k forces ✓). **Not
added:** BP86/B97/OLYP/HCTH (molecular PySCF path); stock VV10 (moved
**blocked** — no CP2K twin); shear/off-diagonal stress (orthorhombic
design); multi-k SGCP / multi-k CPMD (no explicit throw; unpinned
depth hold); GAPW/PAW/periodic hybrid/meta-GGA (blocked); molecular
post-HF / classical MW / CUDA ports (other loops). Roadmap: **no ○/✓
table flips** (nothing new merged since harness wave 10); §8/§10–12 +
MD AIMD open-residual prose + index GPW residual line refreshed;
**VV10** remains ○ with blocked note. **Pending = 9**; **in_progress
= 0**; **blocked = 6**.

**Harness wave 11 2026-07-26 (parent `dev` @ `09298298`):** 9 staffed
slugs merged → `## done` (`ugw-uks-rpbe`, `gpw-ts-forces`, `gpw-xdm-forces`,
`gpw-kpoint-uks-stress-gga`, `cpmd-electronic-thermostat`,
`xl-bomd-uks-kpoints`, `bands-metallic-si`, `dos-spin-si`,
`dfpt-finite-q-spin-smearing`). Specs deleted on close. **Pending residual
= 0** (staff batch emptied the census pending list — **not** queue
exhausted; GapCensus owns next discovery). **in_progress = 0** (no `[~]`;
stock VV10 stays **blocked**). Roadmap: **no new ○→✓ table flips**
(capability rows already ✓ or VV10 still ○/blocked); residual prose closed
UGW@UKS-RPBE, dedicated TS/XDM force pins, multi-k UKS stress GGA expand,
CPMD electronic thermostat, multi-k UKS XL-BOMD, metallic bands, spin DOS,
finite-q DFPT spin+smear. **Blocked = 6** unchanged. **Queue not
exhausted** — GapCensus owns next discovery.

**GapCensus wave 12 2026-07-26 (parent `dev` @ `087df536`):** recover
clean — incomplete_count=**0**; `## in_progress` empty (no `[~]`); open
CP2K/GPW/PW specs on `dev`=**0** (PySCF deferred only:
`periodic-rhf-gamma` / `rsh-eri-cam-b3lyp`); isolation worktrees=**0**;
leftover builds=**none**; `wf/gpw-vv10-ice` tip remains for blocked stock
VV10 docs only (not incomplete). Tree re-diff vs `docs/roadmap/dft.md`
§3/§8/§10–12 + MD AIMD + `tools/cp2k_oracle` +
`include/atomiverse/cpu/{gpw,planewave}`: prior pending residual **0**
(harness wave 11 emptied list — not capability-closed). **True residuals
re-opened (9)** with offline CP2K / literature twins:
`bands-spin-metallic-si` (`Bands::compute` throws collinear spin +
metallic smear — follow-on after fixed-occ spin bands ✓ + restricted
metallic bands ✓),
`dos-spin-metallic-si` (`Dos::compute` throws spin + metallic smear —
follow-on after spin DOS ✓ + metallic DOS ✓),
`gpw-d3-atm` (D3 zero / BJ two-body product; three-body C9 / ATM off —
CP2K `DISPERSION_FUNCTIONAL` ATM),
`gpw-d4-atm` (D4 two-body BJ product + forces/stress ✓; three-body ATM
off — dftd4 / CP2K s9 path),
`xl-bomd-linearized-forces` (SCF-converged XL product at Γ/multi-k RKS/UKS;
zero-SCF linearized-functional force stays ○ — Niklasson XL residual),
`cpmd-kpoints` (CPMD Γ RKS/UKS product; seed requires real `C` —
multi-k Bloch CPMD unpinned),
`second-gen-cpmd-kpoints` (SGCP ASPC on real Γ \(P\) / dual UKS \(P^\sigma\);
multi-k complex \(P_k\) SGCP unpinned),
`periodic-wannier-gpw` (PW Marzari–Vanderbilt MLWF ✓; GPW-grid MLWF
residual still unqueued per §11),
`gpw-fd-phonon-si` (FD supercell phonons product on SW twin only;
GPW analytic-force supercell vs CP2K `VIBRATIONAL_ANALYSIS` residual).
**Not added:** BP86/B97/OLYP/HCTH (molecular PySCF path); stock VV10
(**blocked** — no CP2K twin); shear/off-diagonal stress (orthorhombic
design); periodic hybrid / meta-GGA / GAPW / PAW (**blocked**); multi-k
UMP2 header “follow-on” is **stale** (product `reg_periodic-ump2-kpoints`);
UKS OT+U+smear header residual is **stale** (`reg_ot-uks-dftu-smearing`);
density-scaled TS already product via `hirshfeld`+`ts`; molecular post-HF /
classical MW / CUDA ports (other loops). Roadmap: **no ○/✓ table flips**
(nothing new merged since harness wave 11); §8/§10–12 + MD AIMD open-
residual prose + index GPW residual line refreshed for wave-12 queue.
**Pending = 9**; **in_progress = 0**; **blocked = 6**.

**Assign wave 12 2026-07-26 (parent `dev` @ `cbe74bdf`):** recover clean —
incomplete_count=**0** (no `[~]`; stock VV10 **blocked**); open CP2K/GPW/PW
specs=**0** (PySCF deferred only); unmerged feature `wf/*`=**0**
(`wf/gpw-vv10-ice` tip is blocked docs only). Staffed **all 9** pending
residuals → `## in_progress` (staff 9 ≤ 10; no open-spec resumes; no
pending held). Surface-disjoint notes for workers: `bands-spin-metallic-si`
∥ `dos-spin-metallic-si` share PW `bands_dos` — serialize edits or
split patches carefully; `gpw-d3-atm` ∥ `gpw-d4-atm` share dispersion
tables / force-assembly fold — careful dual surface; `cpmd-kpoints` ∥
`second-gen-cpmd-kpoints` share AIMD propagator / k-seed paths —
coordinate Bloch \(C_k\) / \(P_k\) design; `xl-bomd-linearized-forces`
free of multi-k CPMD/SGCP surface; `periodic-wannier-gpw` free of
`gpw-fd-phonon-si` (MLWF vs vib FD consumer). All picks write-first: no
open `.claude/specs/<slug>/`, no `wf/<slug>` feature tip, no regression
dir yet. **Pending residual after assign = 0** (staff emptied list —
**not** queue exhausted; GapCensus owns next discovery after harness).

## pending

_(empty — all wave-12 census residuals staffed to `[~]`; GapCensus refills
after harness)_

## in_progress

- [~] `bands-spin-metallic-si` | PW band structure with collinear spin **and** metallic occupation smearing | pw | twin: extend `bands-spin-si` / `bands-metallic-si` (CP2K/QE-class spin metal bands) | `Bands::compute` currently rejects spin+smear | assign wave 12
- [~] `dos-spin-metallic-si` | PW DOS with collinear spin **and** metallic occupation smearing | pw | twin: extend `dos-spin-si` / `dos-metallic-si` | `Dos::compute` currently rejects spin+smear | assign wave 12
- [~] `gpw-d3-atm` | DFT-D3 three-body ATM (C9) on GPW | gpw | twin: CP2K Quickstep D3 + ATM on water dimer / ice | two-body zero/BJ product; ATM off in `dftd3.h` | assign wave 12
- [~] `gpw-d4-atm` | DFT-D4 three-body ATM on GPW | gpw | twin: dftd4 ATM / CP2K D4 with \(s_9\neq 0\) | two-body BJ + forces/stress product; ATM off | assign wave 12
- [~] `xl-bomd-linearized-forces` | Zero-SCF linearized-functional XL-BOMD forces | dynamics | twin: Niklasson XL literature / CP2K XL-BOMD low-SCF path | SCF-converged XL ✓ at Γ/multi-k RKS/UKS; linearized force ○ | assign wave 12
- [~] `cpmd-kpoints` | Multi-k Car–Parrinello MD (Bloch \(C_k\)) | dynamics | twin: CP2K CP with k-points | Γ RKS/UKS CPMD ✓; seed needs real `C` (multi-k unpinned) | assign wave 12
- [~] `second-gen-cpmd-kpoints` | Multi-k second-generation CP (ASPC on \(P_k\)) | dynamics | twin: CP2K SGCP / Kühne ASPC + k-mesh | Γ RKS/UKS + OT corrector ✓; multi-k complex density ASPC unpinned | assign wave 12
- [~] `periodic-wannier-gpw` | Maximally localized Wannier functions on GPW | response | twin: CP2K `LOCALIZE` / Wannier on Quickstep | PW MV MLWF ✓ (`reg_periodic-wannier-si` / uks); GPW-grid residual | assign wave 12
- [~] `gpw-fd-phonon-si` | Finite-displacement Γ phonons with GPW analytic forces | response | twin: CP2K `VIBRATIONAL_ANALYSIS` on small cell / GTH | `cpu::vib::Phonon` product on SW twin only (`reg_fd-phonon-supercell`); GPW force consumer residual | assign wave 12

## done (harness wave 11 2026-07-26)

- [x] `ugw-uks-rpbe` | UGW@UKS with RPBE spin XC | gpw | `regressions/gpw/ugw-uks-rpbe` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-ts-forces` | Analytic ionic forces with free-atom TS | gpw | `regressions/gpw/gpw-ts-forces` | ratio **perf_exempt** (E_disp / force pin) | 2026-07-26
- [x] `gpw-xdm-forces` | Analytic ionic forces with free-atom XDM | gpw | `regressions/gpw/gpw-xdm-forces` | ratio **perf_exempt** (E_disp / force pin) | 2026-07-26
- [x] `gpw-kpoint-uks-stress-gga` | Multi-k UKS analytic stress for RevPBE/RPBE/BLYP | gpw | `regressions/gpw/gpw-kpoint-uks-stress-gga` | ratio **perf_exempt** | 2026-07-26
- [x] `cpmd-electronic-thermostat` | Electronic thermostat on CPMD (GPW) | dynamics | `regressions/dynamics/cpmd-electronic-thermostat` | ratio **perf_exempt** (NVE / thermostat envelope) | 2026-07-26
- [x] `xl-bomd-uks-kpoints` | Multi-k UKS XL-BOMD (dual-spin \(P_k^\sigma\)) | dynamics | `regressions/dynamics/xl-bomd-uks-kpoints` | ratio **perf_exempt** (NVE drift) | 2026-07-26
- [x] `bands-metallic-si` | Metallic occupation band structure (PW) | pw | `regressions/pw/bands-metallic-si` | ratio **perf_exempt** | 2026-07-26
- [x] `dos-spin-si` | Collinear-spin DOS (PW) | pw | `regressions/pw/dos-spin-si` | ratio **perf_exempt** | 2026-07-26
- [x] `dfpt-finite-q-spin-smearing` | Finite-q DFPT collinear spin + metallic smear (PW) | response | `regressions/response/dfpt-finite-q-spin-smearing` | ratio **perf_exempt** | 2026-07-26

## done (harness wave 10 2026-07-26)

- [x] `dfpt-spin-smearing` | DFPT collinear spin + metallic occupation smearing (PW) | response | `regressions/response/dfpt-spin-smearing` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-uks-dftu-smearing` | Multi-k UKS OT + Dudarev DFT+U + metallic smear | gpw | `regressions/gpw/ot-uks-dftu-smearing` | ratio **perf_exempt** | 2026-07-26
- [x] `second-gen-cpmd-uks-ot` | SGCP ASPC + OT corrector on Γ-UKS | dynamics | `regressions/dynamics/second-gen-cpmd-uks-ot` | ratio **perf_exempt** (NVE drift) | 2026-07-26
- [x] `xl-bomd-kpoints` | Multi-k XL-BOMD auxiliary-density AIMD | dynamics | `regressions/dynamics/xl-bomd-kpoints` | ratio **perf_exempt** (NVE drift) | 2026-07-26
- [x] `bands-spin-si` | Spin-polarized band structure (PW) | pw | `regressions/pw/bands-spin-si` | ratio **perf_exempt** | 2026-07-26
- [x] `dos-metallic-si` | Metallic occupation DOS (PW) | pw | `regressions/pw/dos-metallic-si` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-wannier-uks` | Spin-polarized periodic MLWF (PW) | response | `regressions/pw/periodic-wannier-uks` | ratio **perf_exempt** | 2026-07-26
- [x] `berry-spin-si` | Spin-polarized Berry-phase polarization (PW) | response | `regressions/pw/berry-spin-si` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-mbd-kpoint-forces` | Multi-k analytic ionic forces with MBD@rsSCS | gpw | `regressions/gpw/gpw-mbd-kpoint-forces` | ratio **perf_exempt** (E_disp / force pin) | 2026-07-26

## done (harness wave 9 2026-07-26)

- [x] `ot-uks-smearing` | UKS OT + metallic occupation smearing (FD/MP) | gpw | `regressions/gpw/ot-uks-smearing` | ratio **~0.12** | 2026-07-26
- [x] `ot-uks-dftu` | Multi-k UKS OT + Dudarev DFT+U | gpw | `regressions/gpw/ot-uks-dftu` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-dftu-smearing` | OT + DFT+U + metallic smearing (RKS multi-k) | gpw | `regressions/gpw/ot-dftu-smearing` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-d4-forces` | Analytic ionic forces with DFT-D4 | gpw | `regressions/gpw/gpw-d4-forces` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-mbd-forces` | Analytic ionic forces with MBD@rsSCS | gpw | `regressions/gpw/gpw-mbd-forces` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-uks-multigrid` | UKS multigrid product collocation (`ngrids>1`) | gpw | `regressions/gpw/gpw-uks-multigrid` | ratio **~1.54** | 2026-07-26
- [x] `dfpt-finite-q-smearing` | Finite-q DFPT with metallic smearing (PW) | response | `regressions/response/dfpt-finite-q-smearing` | ratio **perf_exempt** | 2026-07-26
- [x] `second-gen-cpmd-ot` | Second-gen CP with OT partial-SCF corrector | dynamics | `regressions/dynamics/second-gen-cpmd-ot` | ratio **perf_exempt** (NVE drift) | 2026-07-26

## done (harness wave 8 2026-07-26)

- [x] `periodic-urpa-kpoints` | Multi-k unrestricted direct RPA residual (URPA@UKS + k-mesh) | gpw | `regressions/gpw/periodic-urpa-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-smearing` | OT orbital minimizer + metallic occupation smearing (FD/MP) | gpw | `regressions/gpw/ot-smearing` | ratio **~0.20** | 2026-07-26
- [x] `ot-uks-kpoints` | Multi-k UKS OT minimizer | gpw | `regressions/gpw/ot-uks-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-dftu-kpoints` | Multi-k OT + Dudarev DFT+U | gpw | `regressions/gpw/ot-dftu-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-d4-stress` | Analytic fixed-P cell stress with DFT-D4 | gpw | `regressions/gpw/gpw-d4-stress` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-vdw-scf-potential` | Nonlocal vdW-DF/DF2/rVV10 SCF potential + forces | gpw | `regressions/gpw/gpw-vdw-scf-potential` | ratio **perf_exempt** | 2026-07-26
- [x] `dfpt-metallic-si` | DFPT phonons with metallic occupation smearing (PW) | response | `regressions/response/dfpt-metallic-si` | ratio **perf_exempt** | 2026-07-26
- [x] `dfpt-spin-incommensurate` | Collinear-spin DFPT at incommensurate q (PW) | response | `regressions/response/dfpt-spin-incommensurate` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-optb88-vdw` | optB88 exchange + nonlocal vdW correlation | gpw | `regressions/gpw/gpw-optb88-vdw` | ratio **perf_exempt** (E_disp / XC combo) | 2026-07-26

## done (harness wave 7 2026-07-26)

- [x] `periodic-urpa-si` | Unrestricted direct RPA (URPA@UKS) | gpw | `regressions/gpw/periodic-urpa-si` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-ugw-kpoints` | Multi-k unrestricted spectral G0W0 residual | gpw | `regressions/gpw/periodic-ugw-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-uks-smearing` | Metallic occupation smearing on GPW UKS (FD/MP) | gpw | `regressions/gpw/gpw-uks-smearing` | ratio **~0.017** | 2026-07-26
- [x] `gpw-hirshfeld-uks` | Hirshfeld volume ratios for UKS TS/MBD | gpw | `regressions/gpw/gpw-hirshfeld-uks` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `cpmd-uks` | Car–Parrinello MD on Γ-UKS GPW | dynamics | `regressions/dynamics/cpmd-uks` | ratio **perf_exempt** (NVE drift) | 2026-07-26
- [x] `second-gen-cpmd-uks` | Second-gen CP (ASPC + partial SCF) on Γ-UKS | dynamics | `regressions/dynamics/second-gen-cpmd-uks` | ratio **perf_exempt** (NVE drift) | 2026-07-26
- [x] `ot-kpoints` | OT orbital minimizer with Monkhorst–Pack k-mesh | gpw | `regressions/gpw/ot-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-dftu` | OT + Dudarev DFT+U (Γ) | gpw | `regressions/gpw/ot-dftu` | ratio **~0.26** | 2026-07-26
- [x] `dfpt-spin-si` | DFPT phonons with collinear spin (PW) | response | `regressions/response/dfpt-spin-si` | ratio **perf_exempt** | 2026-07-26

## done (harness wave 6 2026-07-26)

- [x] `gpw-metallic-smearing` | Metallic occupation smearing on GPW (FD/MP multi-k RKS) | gpw | `regressions/gpw/gpw-metallic-smearing` | ratio **~0.011** (ATV ~1.0 s / CP2K 89.5 s; dual-gate free energy) | 2026-07-26
- [x] `gpw-xdm-ice` | Exchange-hole dipole moment (XDM) free-atom dispersion | gpw | `regressions/gpw/gpw-xdm-ice` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `gpw-uks-rpbe` | Γ-UKS RPBE spin-GGA (`XcKind::RPBE`) | gpw | `regressions/gpw/gpw-uks-rpbe` | ratio **~2.70** | 2026-07-26
- [x] `gpw-dftu-kpoint-stress` | Multi-k Dudarev DFT+U analytic cell stress | gpw | `regressions/gpw/gpw-dftu-kpoint-stress` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-ump2-kpoints` | Multi-k unrestricted GPW MP2 residual | gpw | `regressions/gpw/periodic-ump2-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-ugw-si` | Unrestricted spectral G0W0 (UGW@UKS) | gpw | `regressions/gpw/periodic-ugw-si` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-hirshfeld-kpoints` | Hirshfeld volume ratios on multi-k TS/MBD | gpw | `regressions/gpw/gpw-hirshfeld-kpoints` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `dfpt-gga-si` | DFPT phonons with GGA XC kernel (PW) | response | `regressions/response/dfpt-gga-si` | ratio **perf_exempt** | 2026-07-26
- [x] `xl-bomd-uks` | Niklasson XL-BOMD on Γ-UKS GPW | dynamics | `regressions/dynamics/xl-bomd-uks` | ratio **perf_exempt** (NVE drift envelope) | 2026-07-26

## done (harness wave 5 2026-07-26)

- [x] `gpw-dftu-kpoints` | Multi-k Dudarev DFT+U SCF energy (k-mesh + PLUS_U) | gpw | `regressions/gpw/gpw-dftu-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-gw-finite-q` | Finite-q inter-k \(W\) spectral G0W0 residual | gpw | `regressions/gpw/periodic-gw-finite-q` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-uks-revpbe` | Γ-UKS revPBE spin-GGA (`XcKind::RevPBE`) | gpw | `regressions/gpw/gpw-uks-revpbe` | ratio **~3.22** | 2026-07-26
- [x] `cpmd-rattle-water` | CPMD electronic RATTLE (constraint-preserving) | dynamics | `regressions/dynamics/cpmd-rattle-water` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-ump2-si` | Unrestricted GPW MP2 (UKS reference) | gpw | `regressions/gpw/periodic-ump2-si` | ratio **perf_exempt** | 2026-07-26
- [x] `ot-uks-o2` | OT orbital minimizer for Γ-UKS | gpw | `regressions/gpw/ot-uks-o2` | ratio **~8.50** | 2026-07-26

## done (harness wave 4 2026-07-26)

- [x] `gpw-kpoint-uks-stress` | Multi-k UKS analytic cell stress | gpw | `regressions/gpw/gpw-kpoint-uks-stress` | ratio **perf_exempt** | 2026-07-26
- [x] `periodic-gw-kpoints` | Multi-k k-diagonal spectral G0W0 | gpw | `regressions/gpw/periodic-gw-kpoints` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-vdw-df2-ice` | Nonlocal vdW-DF2 (LMKLL) correlation | gpw | `regressions/gpw/gpw-vdw-df2-ice` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `gpw-rvv10-ice` | Nonlocal rVV10 correlation | gpw | `regressions/gpw/gpw-rvv10-ice` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `gpw-mbd-periodic` | Periodic Ewald + k-mesh MBD@rsSCS | gpw | `regressions/gpw/gpw-mbd-periodic` | ratio **perf_exempt** (E_disp pin) | 2026-07-26
- [x] `pw-spin-gga` | True collinear spin-GGA on plane-wave grid | pw | `regressions/pw/pw-spin-gga` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-uks-blyp` | Γ-UKS BLYP spin-GGA | gpw | `regressions/gpw/gpw-uks-blyp` | ratio **~2.48** | 2026-07-26
- [x] `periodic-rpa-finite-q` | Finite-q inter-k RI-GPW dRPA | gpw | `regressions/gpw/periodic-rpa-finite-q` | ratio **perf_exempt** | 2026-07-26
- [x] `gpw-dftu-stress` | Analytic stress with Dudarev DFT+U | gpw | `regressions/gpw/gpw-dftu-stress` | ratio **perf_exempt** | 2026-07-26

## done (harness wave 3 2026-07-25)

- [x] `gpw-rpbe-water` | GPW RPBE GGA (`XcKind::RPBE`) | gpw | `regressions/gpw/gpw-rpbe-water` | ratio **perf_exempt** (dual-gate \|E_RPBE−E_PBE\|) | 2026-07-25
- [x] `gpw-kpoint-uks-forces` | Multi-k UKS analytic ionic forces | gpw | `regressions/gpw/gpw-kpoint-uks-forces` | ratio **perf_exempt** | 2026-07-25
- [x] `gpw-uks-stress` | Γ-UKS analytic fixed-P cell stress | gpw | `regressions/gpw/gpw-uks-stress` | ratio **perf_exempt** (multi-k UKS stress residual) | 2026-07-25
- [x] `gpw-vdw-df-ice` | Nonlocal vdW-DF (DRSLL) correlation | gpw | `regressions/gpw/gpw-vdw-df-ice` | ratio **perf_exempt** (E_disp pin) | 2026-07-25
- [x] `gpw-mbd-hirshfeld` | Hirshfeld volume ratios for TS/MBD | gpw | `regressions/gpw/gpw-mbd-hirshfeld` | ratio **perf_exempt** (E_disp + V/V_free) | 2026-07-25
- [x] `periodic-rpa-kpoints` | Multi-k RI-GPW dRPA residual | gpw | `regressions/gpw/periodic-rpa-kpoints` | ratio **perf_exempt** | 2026-07-25
- [x] `periodic-gw-si` | Γ-point spectral G0W0 (periodic Si) | gpw | `regressions/gpw/periodic-gw-si` | ratio **perf_exempt** | 2026-07-25
- [x] `dfpt-incommensurate-si` | Incommensurate-q DFPT phonons (PW) | response | `regressions/response/dfpt-incommensurate-si` | ratio **perf_exempt** | 2026-07-25

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

- [!] `gpw-vv10-ice` | blocked — stock VV10 (≠ rVV10): CP2K 2026.1 has no `NON_LOCAL TYPE VV10` (enum DRSLL|LMKLL|RVV10 only; no VV10 kernel table); reconfirmed waves 6–10 on `wf/gpw-vv10-ice` @ `a60c73d5`; product VV-family remains rVV10; no offline CP2K twin — scope_out of pending
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
