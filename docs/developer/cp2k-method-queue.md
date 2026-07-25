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
(GPW `XcKind` still PADE|PBE|PBEsol; `gpw_stress` still strain-FD; no
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

## pending

- [ ] `gpw-variable-cell-molopt` | GPW CELL_OPT MOLOPT re-twin | gpw | ice Ih | DZVP-MOLOPT-GTH | a\* Å | pin-quality follow-up declared in `gpw-variable-cell-ice` done row ("MOLOPT re-twin queued"; GTH-chain pin shipped 2026-07-25); CP2K CELL_OPT same input, MOLOPT basis; low priority; **dep:** `gpw-variable-cell-ice` ✓
- [ ] `gpw-blyp-water` | GPW BLYP GGA | gpw | water | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 BLYP molecular ✓ / GPW `XcKind` still PADE|PBE|PBEsol only; CP2K XC_FUNCTIONAL BLYP
- [ ] `gpw-revpbe-si` | GPW revPBE GGA | gpw | Si | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 revPBE ○; CP2K XC_FUNCTIONAL revPBE textbook bulk
- [ ] `gpw-rpbe-water` | GPW RPBE GGA | gpw | water | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 RPBE ○; CP2K XC_FUNCTIONAL RPBE; shares `XcKind` plumbing — prefer after `gpw-revpbe-si` (merge-conflict avoidance)
- [ ] `pw-uks-fe-smearing` | bulk Fe UKS + metallic smearing | pw | Fe bcc | NC-UPF | e_tot / mag | roadmap §10 follow-on after O₂ k-mesh UKS ✓ (`pw-uks-kmesh-fe`); CP2K/QE-class Fe spin+smearing
- [ ] `pw-variable-cell-si` | PW variable-cell relax | pw | Si | NC-UPF | cell / stress | optimization.md PW stress productization; Nielsen–Martin ✓; CP2K CELL_OPT twin; **dep:** `pw-analytic-stress-si` ✓
- [ ] `gpw-ts-ice` | Tkatchenko–Scheffler dispersion | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp | roadmap §8 TS ○; CP2K VDW_POTENTIAL / PAIR_POTENTIAL; **dep:** D3/D4 pin patterns ✓
- [ ] `gpw-mbd-ice` | Many-body dispersion (MBD) | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp | roadmap §8 MBD ○; CP2K MBD; prefer after TS foothold
- [ ] `gpw-vdw-df-ice` | Nonlocal vdW-DF correlation | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp / e_tot Ha | roadmap §8 vdW-DF ○; CP2K VDW_POTENTIAL NON_LOCAL (DRSLL); **after** dispersion footholds (`gpw-ts-ice` / `gpw-mbd-ice`)
- [ ] `phonon-dispersion-si` | Phonon dispersion (q-path) | response | Si | NC-UPF / GTH | ω(q) cm⁻¹ | roadmap §11 full dispersions ○ (Γ FD ✓, Γ DFPT ✓); CP2K VIBRATIONAL_ANALYSIS / DFPT path; **dep:** `fd-phonon-supercell` ✓ + `dfpt-phonon-si` ✓
- [ ] `second-gen-cpmd-water` | Second-generation Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ 2nd-gen CP; CP2K; **after** `cpmd-water`
- [ ] `periodic-ri-mp2-si` | RI / k-point periodic MP2 productization | gpw | Si | MOLOPT+GTH | e_corr Ha | roadmap §10 Periodic MP2 stays ○ (Γ canonical foothold ✓ 2026-07-25 — foothold ≠ capability closed); CP2K &WF_CORRELATION RI-MP2 (Del Ben 2012); closes the capability row; **dep:** `periodic-mp2-si` ✓
- [ ] `periodic-rpa-si` | Periodic RPA | gpw | Si | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic RPA; CP2K RPA; **after MP2 foothold** (`periodic-mp2-si` ✓)
- [ ] `periodic-gw-si` | Periodic GW | gpw / response | Si (tiny cell) | MOLOPT+GTH | QP gaps | roadmap §10 ○ Periodic GW; CP2K GW; **dep:** solid GPW + bands foothold ✓; not molecular G0W0 (PySCF loop)

## in_progress

<!-- assign wave 3 2026-07-24: 1 resume + 9 new; staff 10 -->
<!-- assign wave (claude) 2026-07-25: staffed 6 of 10 [~] — gpw-kpoints-si,
     gpw-variable-cell-ice, dfpt-phonon-si, xl-bomd-water, berry-polarization,
     periodic-mp2-si; wannier-centers / cpmd-water / mpi-*-scaling held -->
<!-- harness wave 1 (claude) 2026-07-25: 5 of the 6 staffed merged → done;
     xl-bomd-water WIP-checkpointed -->
<!-- assign wave (claude) 2026-07-25 (wave 2): staffed all 5 [~] resumes +
     promoted gpw-analytic-stress from pending; second-gen-cpmd-water skipped
     (dep cpmd-water is a pick this wave) -->
<!-- 2026-07-25: xl-bomd-water verified green + merged → done -->
<!-- 2026-07-25: wannier-centers Foster–Boys foothold on branch wf/wannier-centers -->
- [~] `cpmd-water` | Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ CPMD; CP2K; after solid BOMD pins ✓; wave 3 new
- [~] `mpi-gpw-scaling-4rank` | MPI GPW energy scaling | infra | water | MOLOPT | wall ratio | industrial gate; may start `perf_exempt` until domain decomp lands; wave 3 new
- [~] `mpi-pw-scaling-4rank` | MPI PW k-mesh scaling | infra | Si | NC-UPF | wall ratio | k-point RR ships; multi-rank industrial gate; `perf_exempt` OK; wave 3 new
- [~] `gpw-analytic-stress` | Analytic GPW Pulay stress | gpw | water | DZVP-MOLOPT-GTH | stress Ha/Bohr³ | roadmap §11 productization after strain-FD ✓; CP2K STRESS_TENSOR ANALYTICAL; **dep:** `gpw-stress-water` ✓; promoted assign wave (claude) 2026-07-25 (wave 2)

## blocked

- [!] `periodic-hybrid-gpw` | blocked — no periodic exact exchange (same as periodic HF); HSE/PBE0 periodic ○
- [!] `gapw-dft` | blocked — GAPW not in tree; long-horizon
- [!] `ultrasoft-paw` | blocked — UPF NC + GTH only; PAW/USPP ○
- [!] `cuda-gpw` | blocked — GPW CPU-only; CUDA plane-wave / Kokkos are separate tracks (not this loop)
- [!] `periodic-meta-gga` | blocked — no kinetic-energy density (τ) channel in GPW/PW grid engines

## done

- [x] `wannier-centers` | Maximally localized Wannier (Foster–Boys foothold) | scf | `regressions/scf/wannier-centers` | ratio **perf_exempt** (H₂ midpoint + H₂O rabbit-ear geometry; no CP2K LOCALIZE twin) | 2026-07-25 | molecular FB Jacobi on RHF; periodic MV remains ○
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
