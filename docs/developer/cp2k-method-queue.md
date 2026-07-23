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
**Recover reconfirm 2026-07-23:** **16 `[~]` incomplete** under `## in_progress`
still own assign priority (INDEX Live empty; no CP2K-scope open specs; local
`wf/*` are PySCF molecular only; no `regressions/gpw/`; PW footholds are
`regressions/pw/h_*` only). Do not promote pending until incomplete clear; do
not demote incomplete. Pending below = open roadmap ○ with offline CP2K twins,
**excluding** slugs already in in_progress / done / blocked and excluding
GAPW/PAW/periodic hybrid/meta-GGA/CUDA (see `## blocked`). Later wave order
(after incomplete): k-point GPW → k-point UKS PW → GPW cell_opt → extra GGA/D4
→ DFPT → Berry/Wannier → XL-BOMD/CPMD → periodic MP2/RPA/GW → MPI industrial.

## pending

- [ ] `gpw-kpoints-si` | k-point GPW RKS | gpw | Si | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha | roadmap: GPW Γ-only / no Gaussian k-path; CP2K Quickstep k-points; **dep:** Γ energy/forces pins (`gpw-water-pbe-mt`, `gpw-forces-water`)
- [ ] `pw-uks-kmesh-fe` | k-point UKS plane-wave | pw | Fe (or O₂ cell) | NC-UPF | e_tot / mag | roadmap ○ k-point periodic UKS; **dep:** `pw-spin-polarized-o2` + `pw-smearing-al`
- [ ] `gpw-variable-cell-ice` | GPW variable-cell relax | gpw | ice Ih | DZVP-MOLOPT-GTH | cell / stress | productization after GPW stress; CP2K CELL_OPT; **dep:** `gpw-stress-water`
- [ ] `gpw-pbesol-si` | GPW PBEsol GGA | gpw | Si | DZVP-MOLOPT-GTH | e_tot Ha | roadmap §3 GGA ○ PBEsol; CP2K XC_FUNCTIONAL PBEsol textbook bulk
- [ ] `gpw-d4-water` | GPW DFT-D4 periodic | gpw | water | DZVP-MOLOPT-GTH | e_tot / E_disp | roadmap §8 D4 ○; CP2K DFTD4; **dep:** GPW D3 pin pattern (`gpw-pbe-d3bj-ice`)
- [ ] `dfpt-phonon-si` | DFPT phonons | response | Si | NC-UPF / GTH | ω cm⁻¹ | roadmap §11 DFPT ○; CP2K DFPT; **dep:** `fd-phonon-supercell`
- [ ] `berry-polarization` | Berry-phase polarization | response | Si / BaTiO₃ toy | — | P | roadmap §11 ○; CP2K Berry phase; after k-path / `bands-dos-si`
- [ ] `wannier-centers` | Maximally localized Wannier | response | water / Si | — | centers | roadmap §11 ○; CP2K Wannier; **dep:** `bands-dos-si`
- [ ] `xl-bomd-water` | Extended-Lagrangian BOMD | dynamics | water | DZVP | energy drift | MD roadmap ○ XL-BOMD (ASPC ≠ XL); CP2K; **dep:** `gpw-bomd-water`
- [ ] `cpmd-water` | Car–Parrinello MD | dynamics | water | DZVP | energy drift | MD roadmap ○ CPMD; CP2K; after solid BOMD pins
- [ ] `periodic-mp2-si` | Periodic RI-MP2 (GPW) | gpw | Si (tiny cell) | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic MP2; CP2K MP2; **dep:** solid Γ GPW energy
- [ ] `periodic-rpa-si` | Periodic RPA | gpw | Si | MOLOPT+GTH | e_tot Ha | roadmap §10 ○ Periodic RPA; CP2K RPA; after MP2 foothold
- [ ] `periodic-gw-si` | Periodic GW | gpw / response | Si (tiny cell) | MOLOPT+GTH | QP gaps | roadmap §10 ○ Periodic GW; CP2K GW; **dep:** solid GPW + bands foothold; not molecular G0W0 (PySCF loop)
- [ ] `mpi-gpw-scaling-4rank` | MPI GPW energy scaling | infra | water | MOLOPT | wall ratio | industrial gate; may start `perf_exempt` until domain decomp lands
- [ ] `mpi-pw-scaling-4rank` | MPI PW k-mesh scaling | infra | Si | NC-UPF | wall ratio | k-point RR ships; multi-rank industrial gate; `perf_exempt` OK

## in_progress

- [~] `gpw-water-pbe-mt` | Γ-RKS GPW energy (MT water) | gpw | water | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha | product in tutorial ch09; **no** `regressions/gpw/` yet; CP2K getting-started MT pin
- [~] `gpw-forces-water` | GPW analytic forces | gpw | water | DZVP-MOLOPT-GTH | forces Ha/Bohr | Hellmann–Feynman + Pulay ✓ in engine; e2e vs CP2K FORCE_EVAL; watch `local_pp_on_grid`
- [~] `gpw-uks-o2-gamma` | Γ-point UKS GPW | gpw | O2 | DZVP-MOLOPT-GTH | e_tot / ⟨S²⟩ | roadmap ○ Γ-UKS; closed-shell only today; CP2K Quickstep UKS
- [~] `gpw-stress-water` | GPW cell stress | gpw | water | DZVP-MOLOPT-GTH | stress Ha/Bohr³ | PW strain-FD stress ✓; **GPW stress open**; CP2K cell stress twin; **dep:** energy+forces pins
- [~] `gpw-pbe-d3bj-ice` | GPW PBE-D3(BJ) periodic | gpw | ice Ih / water dimer cell | MOLOPT+GTH | e_tot / E_disp | zero-damping D3 on GPW ✓; BJ on molecular ✓; BJ on GPW still gap
- [~] `pw-lda-si` | PW LDA Si k-mesh energy | pw | Si | NC-UPF | e_tot Ry/eV | product PW LDA + MP mesh; re-pin offline CP2K (or ABINIT-class) under this loop
- [~] `pw-pbe-gga-si` | PW PBE (true GGA) | pw | Si | NC-UPF | e_tot / forces | PW XC still LDA-limit (`sigma=0`); CP2K/QE PBE twin
- [~] `pw-smearing-al` | metallic smearing PW | pw | Al | NC-UPF | e_tot / occupations | roadmap ○ smearing for metals; CP2K Fermi/Methfessel–Paxton
- [~] `pw-spin-polarized-o2` | PW spin-polarized KS | pw | O2 | NC-UPF | e_tot / mag | PW spin-unpolarized only; k-point UKS ○; CP2K LSD
- [~] `pw-analytic-stress-si` | PW Nielsen–Martin analytic stress | pw | Si | NC-UPF | stress eV/Å³ | strain-FD stress ships; analytic form is follow-on; CP2K STRESS_TENSOR; **dep:** `pw-lda-si`
- [~] `gpw-bomd-water` | GPW BOMD short traj | dynamics | water | DZVP | energy drift | GPWEvaluator forces + MD marked ✓; need CP2K MD reference trajectory pin; **dep:** `gpw-forces-water`
- [~] `pw-bomd-h` | PW BOMD foothold | dynamics | H / Si | UPF | energy drift | PW analytic forces + MD ✓; short NVE vs CP2K
- [~] `fd-phonon-supercell` | FD phonons (supercell) | response | Si or H2O ice | — | ω cm⁻¹ | molecular `vib::Hessian` ≠ phonons; CP2K VIBRATIONAL_ANALYSIS / supercell
- [~] `ot-minimizer-water` | OT vs DIIS large cell | gpw | water box | MOLOPT | e_tot / iters | OT minimizer exists; productize + CP2K OT parity case
- [~] `gpw-dftu-nio` | Periodic DFT+U | gpw | NiO | MOLOPT+GTH | e_tot / mag | roadmap ○ DFT+U; CP2K PLUS_U textbook antiferromagnet
- [~] `bands-dos-si` | Band structure / DOS | pw | Si | NC-UPF | ε(k) / DOS | roadmap ○ post-processing; CP2K BAND_STRUCTURE / DOS; **enables** wannier later

## blocked

- [!] `periodic-hybrid-gpw` | blocked — no periodic exact exchange (same as periodic HF); HSE/PBE0 periodic ○
- [!] `gapw-dft` | blocked — GAPW not in tree; long-horizon
- [!] `ultrasoft-paw` | blocked — UPF NC + GTH only; PAW/USPP ○
- [!] `cuda-gpw` | blocked — GPW CPU-only; CUDA plane-wave / Kokkos are separate tracks (not this loop)
- [!] `periodic-meta-gga` | blocked — no kinetic-energy density (τ) channel in GPW/PW grid engines

## done

(seed later from harness; existing product not yet re-pinned under this loop:
Γ-RKS GPW energy/forces + D3 zero, PW LDA k-mesh, PW analytic forces, PW
strain-FD stress, MT Poisson, OT scaffolding — document when each gets a
`regressions/{gpw,pw}/<slug>/` offline CP2K pin)
