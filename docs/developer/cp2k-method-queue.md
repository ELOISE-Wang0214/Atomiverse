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

## pending

- [ ] `gpw-analytic-stress` | Analytic GPW Pulay stress | gpw | water | DZVP-MOLOPT-GTH | stress Ha/Bohr³ | roadmap §11 productization after strain-FD ✓; CP2K STRESS_TENSOR ANALYTICAL; **dep:** `gpw-stress-water` ✓
- [ ] `gpw-blyp-water` | GPW BLYP GGA | gpw | water | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 BLYP molecular ✓ / GPW `XcKind` still PADE|PBE|PBEsol only; CP2K XC_FUNCTIONAL BLYP
- [ ] `gpw-revpbe-si` | GPW revPBE GGA | gpw | Si | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 revPBE ○; CP2K XC_FUNCTIONAL revPBE textbook bulk
- [ ] `gpw-rpbe-water` | GPW RPBE GGA | gpw | water | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 RPBE ○; CP2K XC_FUNCTIONAL RPBE; shares `XcKind` plumbing — prefer after `gpw-revpbe-si` (merge-conflict avoidance)
- [ ] `pw-uks-fe-smearing` | bulk Fe UKS + metallic smearing | pw | Fe bcc | NC-UPF | e_tot / mag | roadmap §10 follow-on after O₂ k-mesh UKS ✓ (`pw-uks-kmesh-fe`); CP2K/QE-class Fe spin+smearing
- [ ] `pw-variable-cell-si` | PW variable-cell relax | pw | Si | NC-UPF | cell / stress | optimization.md PW stress productization; Nielsen–Martin ✓; CP2K CELL_OPT twin; **dep:** `pw-analytic-stress-si` ✓
- [ ] `gpw-ts-ice` | Tkatchenko–Scheffler dispersion | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp | roadmap §8 TS ○; CP2K VDW_POTENTIAL / PAIR_POTENTIAL; **dep:** D3/D4 pin patterns ✓
- [ ] `gpw-mbd-ice` | Many-body dispersion (MBD) | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp | roadmap §8 MBD ○; CP2K MBD; prefer after TS foothold
- [ ] `gpw-vdw-df-ice` | Nonlocal vdW-DF correlation | gpw | ice Ih | DZVP-MOLOPT-GTH | E_disp / e_tot Ha | roadmap §8 vdW-DF ○; CP2K VDW_POTENTIAL NON_LOCAL (DRSLL); **after** dispersion footholds (`gpw-ts-ice` / `gpw-mbd-ice`)
- [ ] `phonon-dispersion-si` | Phonon dispersion (q-path) | response | Si | NC-UPF / GTH | ω(q) cm⁻¹ | roadmap §11 full dispersions ○ (Γ FD ✓); CP2K VIBRATIONAL_ANALYSIS / DFPT path; **dep:** `fd-phonon-supercell` ✓; prefer after `dfpt-phonon-si`
- [ ] `second-gen-cpmd-water` | Second-generation Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ 2nd-gen CP; CP2K; **after** `cpmd-water`
- [ ] `periodic-rpa-si` | Periodic RPA | gpw | Si | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic RPA; CP2K RPA; **after MP2 foothold** (`periodic-mp2-si`)
- [ ] `periodic-gw-si` | Periodic GW | gpw / response | Si (tiny cell) | MOLOPT+GTH | QP gaps | roadmap §10 ○ Periodic GW; CP2K GW; **dep:** solid GPW + bands foothold ✓; not molecular G0W0 (PySCF loop)

## in_progress

<!-- assign wave 3 2026-07-24: 1 resume + 9 new; staff 10 -->
<!-- assign wave (claude) 2026-07-25: staffed 6 of 10 [~] — gpw-kpoints-si,
     gpw-variable-cell-ice, dfpt-phonon-si, xl-bomd-water, berry-polarization,
     periodic-mp2-si; wannier-centers / cpmd-water / mpi-*-scaling held -->
- [~] `gpw-kpoints-si` | k-point GPW RKS | gpw | Si | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha | **resume** wave 2 leftover; roadmap: GPW Γ-only / no Gaussian k-path; CP2K Quickstep k-points; **dep:** Γ energy/forces pins ✓; **no open spec / no wf merge / no regression** — write-first on resume
- [~] `gpw-variable-cell-ice` | GPW variable-cell relax | gpw | ice Ih | DZVP-MOLOPT-GTH | cell / stress | productization after GPW stress; CP2K CELL_OPT; **dep:** `gpw-stress-water` ✓; wave 3 new
- [~] `dfpt-phonon-si` | DFPT phonons | response | Si | NC-UPF / GTH | ω cm⁻¹ | roadmap §11 DFPT ○; CP2K DFPT; **dep:** `fd-phonon-supercell` ✓; wave 3 new
- [~] `xl-bomd-water` | Extended-Lagrangian BOMD | dynamics | water | DZVP | energy drift | MD roadmap ○ XL-BOMD (ASPC ≠ XL); CP2K; **dep:** `gpw-bomd-water` ✓; wave 3 new
- [~] `berry-polarization` | Berry-phase polarization | response | Si / BaTiO₃ toy | — | P | roadmap §11 ○; CP2K Berry phase; after k-path / `bands-dos-si` ✓; wave 3 new
- [~] `wannier-centers` | Maximally localized Wannier | response | water / Si | — | centers | roadmap §11 ○; CP2K Wannier; **dep:** `bands-dos-si` ✓; wave 3 new
- [~] `cpmd-water` | Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ CPMD; CP2K; after solid BOMD pins ✓; wave 3 new
- [~] `periodic-mp2-si` | Periodic RI-MP2 (GPW) | gpw | Si (tiny cell) | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic MP2; CP2K MP2; **dep:** solid Γ GPW energy ✓; wave 3 new
- [~] `mpi-gpw-scaling-4rank` | MPI GPW energy scaling | infra | water | MOLOPT | wall ratio | industrial gate; may start `perf_exempt` until domain decomp lands; wave 3 new
- [~] `mpi-pw-scaling-4rank` | MPI PW k-mesh scaling | infra | Si | NC-UPF | wall ratio | k-point RR ships; multi-rank industrial gate; `perf_exempt` OK; wave 3 new

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
