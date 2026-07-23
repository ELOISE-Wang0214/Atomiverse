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

Census: 2026-07-23 (parent `dev`). Source: `docs/roadmap/dft.md` §10–12 +
§3/§8 GGA/dispersion ○, `docs/roadmap/index.md` GPW/PW caveats,
`docs/roadmap/molecular-dynamics.md` §AIMD ○ (XL-BOMD / CPMD).
**Harness wave 1 2026-07-23:** 10 staff-batch slugs merged → `## done` (see
below). **6 `[~]` incomplete** remain under `## in_progress` (deps / later
slots). Pending below = open roadmap ○ with offline CP2K twins, **excluding**
slugs already in in_progress / done / blocked and excluding
GAPW/PAW/periodic hybrid/meta-GGA/CUDA (see `## blocked`). Later wave order
(after incomplete): k-point GPW → k-point UKS PW → GPW cell_opt → extra GGA/D4
→ DFPT → Berry/Wannier → XL-BOMD/CPMD → periodic MP2/RPA/GW → MPI industrial.

## pending

- [ ] `gpw-kpoints-si` | k-point GPW RKS | gpw | Si | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha | roadmap: GPW Γ-only / no Gaussian k-path; CP2K Quickstep k-points; **dep:** Γ energy/forces pins ✓ (`gpw-water-pbe-mt`, `gpw-forces-water`)
- [ ] `pw-uks-kmesh-fe` | k-point UKS plane-wave | pw | Fe (or O₂ cell) | NC-UPF | e_tot / mag | roadmap ○ k-point periodic UKS; **dep:** `pw-spin-polarized-o2` ✓ + `pw-smearing-al` ✓
- [ ] `gpw-variable-cell-ice` | GPW variable-cell relax | gpw | ice Ih | DZVP-MOLOPT-GTH | cell / stress | productization after GPW stress; CP2K CELL_OPT; **dep:** `gpw-stress-water`
- [ ] `gpw-pbesol-si` | GPW PBEsol GGA | gpw | Si | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 GGA ○ PBEsol; CP2K XC_FUNCTIONAL PBEsol textbook bulk
- [ ] `gpw-d4-water` | GPW DFT-D4 periodic | gpw | water | DZVP-MOLOPT-GTH | e_tot / E_disp | roadmap §8 D4 ○; CP2K DFTD4; **dep:** GPW D3 pin pattern ✓ (`gpw-pbe-d3bj-ice`)
- [ ] `dfpt-phonon-si` | DFPT phonons | response | Si | NC-UPF / GTH | ω cm⁻¹ | roadmap §11 DFPT ○; CP2K DFPT; **dep:** `fd-phonon-supercell`
- [ ] `berry-polarization` | Berry-phase polarization | response | Si / BaTiO₃ toy | — | P | roadmap §11 ○; CP2K Berry phase; after k-path / `bands-dos-si`
- [ ] `wannier-centers` | Maximally localized Wannier | response | water / Si | — | centers | roadmap §11 ○; CP2K Wannier; **dep:** `bands-dos-si`
- [ ] `xl-bomd-water` | Extended-Lagrangian BOMD | dynamics | water | DZVP | energy drift | MD roadmap ○ XL-BOMD (ASPC ≠ XL); CP2K; **dep:** `gpw-bomd-water`
- [ ] `cpmd-water` | Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ CPMD; CP2K; after solid BOMD pins
- [ ] `periodic-mp2-si` | Periodic RI-MP2 (GPW) | gpw | Si (tiny cell) | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic MP2; CP2K MP2; **dep:** solid Γ GPW energy ✓
- [ ] `periodic-rpa-si` | Periodic RPA | gpw | Si | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic RPA; CP2K RPA; after MP2 foothold
- [ ] `periodic-gw-si` | Periodic GW | gpw / response | Si (tiny cell) | MOLOPT+GTH | QP gaps | roadmap §10 ○ Periodic GW; CP2K GW; **dep:** solid GPW + bands foothold; not molecular G0W0 (PySCF loop)
- [ ] `mpi-gpw-scaling-4rank` | MPI GPW energy scaling | infra | water | MOLOPT | wall ratio | industrial gate; may start `perf_exempt` until domain decomp lands
- [ ] `mpi-pw-scaling-4rank` | MPI PW k-mesh scaling | infra | Si | NC-UPF | wall ratio | k-point RR ships; multi-rank industrial gate; `perf_exempt` OK

## in_progress

- [~] `gpw-stress-water` | GPW cell stress | gpw | water | DZVP-MOLOPT-GTH | stress Ha/Bohr³ | PW strain-FD stress ✓; **GPW stress open**; CP2K cell stress twin; **dep:** energy+forces pins ✓
- [~] `pw-analytic-stress-si` | PW Nielsen–Martin analytic stress | pw | Si | NC-UPF | stress eV/Å³ | strain-FD stress ships; analytic form is follow-on; CP2K STRESS_TENSOR; **dep:** `pw-lda-si` ✓
- [~] `gpw-bomd-water` | GPW BOMD short traj | dynamics | water | DZVP | energy drift | GPWEvaluator forces + MD marked ✓; need CP2K MD reference trajectory pin; **dep:** `gpw-forces-water` ✓
- [~] `fd-phonon-supercell` | FD phonons (supercell) | response | Si or H2O ice | — | ω cm⁻¹ | molecular `vib::Hessian` ≠ phonons; CP2K VIBRATIONAL_ANALYSIS / supercell
- [~] `gpw-dftu-nio` | Periodic DFT+U | gpw | NiO | MOLOPT+GTH | e_tot / mag | roadmap ○ DFT+U; CP2K PLUS_U textbook antiferromagnet
- [~] `bands-dos-si` | Band structure / DOS | pw | Si | NC-UPF | ε(k) / DOS | roadmap ○ post-processing; CP2K BAND_STRUCTURE / DOS; **enables** wannier later

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
