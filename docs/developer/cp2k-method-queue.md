# CP2K capability queue

Durable worklist for workflow `cp2k-method-loop`. Do not delete done rows.
Census rewrites **pending** from `docs/roadmap/dft.md` §10–12 (periodic / GPW /
PW / response / infra) and related MD ○ footholds that **CP2K can oracle
offline** (Quickstep GPW, plane-wave DFT, forces/stress, OT, D3, BOMD, phonons).
**done** / **blocked** history is preserved.

Out of scope for this queue (use `pyscf-method-queue.md`): molecular Gaussian
post-HF, pure PySCF continuum/charge methods without CP2K twin, classical MD/FF,
CUDA-Kokkos product tracks.

Census: 2026-07-23 (parent `dev`). Source: `docs/roadmap/dft.md` §10–12,
`docs/roadmap/index.md` GPW/PW notes, `docs/roadmap/molecular-dynamics.md` §AIMD.
Priority: official e2e re-pins of shipping Γ-RKS GPW / PW H·Si footholds →
open-shell GPW → GPW stress → PW true GGA → metallic smearing → BOMD / phonons
→ OT product depth → DFT+U / bands → MPI industrial later.

## pending

- [ ] `gpw-water-pbe-mt` | Γ-RKS GPW energy (MT water) | gpw | water | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha | product in tutorial ch09; **no** `regressions/gpw/` yet; CP2K getting-started MT pin
- [ ] `gpw-forces-water` | GPW analytic forces | gpw | water | DZVP-MOLOPT-GTH | forces Ha/Bohr | Hellmann–Feynman + Pulay ✓ in engine; e2e vs CP2K FORCE_EVAL; watch `local_pp_on_grid`
- [ ] `gpw-uks-o2-gamma` | Γ-point UKS GPW | gpw | O2 | DZVP-MOLOPT-GTH | e_tot / ⟨S²⟩ | roadmap ○ Γ-UKS; closed-shell only today; CP2K Quickstep UKS
- [ ] `gpw-stress-water` | GPW cell stress | gpw | water | DZVP-MOLOPT-GTH | stress Ha/Bohr³ | PW strain-FD stress ✓; **GPW stress open**; CP2K cell stress twin
- [ ] `gpw-pbe-d3bj-ice` | GPW PBE-D3(BJ) periodic | gpw | ice Ih / water dimer cell | MOLOPT+GTH | e_tot / E_disp | zero-damping D3 on GPW ✓; BJ on molecular ✓; BJ on GPW still gap
- [ ] `pw-lda-si` | PW LDA Si k-mesh energy | pw | Si | NC-UPF | e_tot Ry/eV | product PW LDA + MP mesh; re-pin offline CP2K (or ABINIT-class) under this loop
- [ ] `pw-pbe-gga-si` | PW PBE (true GGA) | pw | Si | NC-UPF | e_tot / forces | PW XC still LDA-limit (`sigma=0`); CP2K/QE PBE twin
- [ ] `pw-smearing-al` | metallic smearing PW | pw | Al | NC-UPF | e_tot / occupations | roadmap ○ smearing for metals; CP2K Fermi/Methfessel–Paxton
- [ ] `pw-spin-polarized-o2` | PW spin-polarized KS | pw | O2 | NC-UPF | e_tot / mag | PW spin-unpolarized only; k-point UKS ○; CP2K LSD
- [ ] `pw-analytic-stress-si` | PW Nielsen–Martin analytic stress | pw | Si | NC-UPF | stress eV/Å³ | strain-FD stress ships; analytic form is follow-on; CP2K STRESS_TENSOR
- [ ] `gpw-bomd-water` | GPW BOMD short traj | dynamics | water | DZVP | energy drift | GPWEvaluator forces + MD marked ✓; need CP2K MD reference trajectory pin
- [ ] `pw-bomd-h` | PW BOMD foothold | dynamics | H / Si | UPF | energy drift | PW analytic forces + MD ✓; short NVE vs CP2K
- [ ] `fd-phonon-supercell` | FD phonons (supercell) | response | Si or H2O ice | — | ω cm⁻¹ | molecular `vib::Hessian` ≠ phonons; CP2K VIBRATIONAL_ANALYSIS / supercell
- [ ] `ot-minimizer-water` | OT vs DIIS large cell | gpw | water box | MOLOPT | e_tot / iters | OT minimizer exists; productize + CP2K OT parity case
- [ ] `gpw-dftu-nio` | Periodic DFT+U | gpw | NiO | MOLOPT+GTH | e_tot / mag | roadmap ○ DFT+U; CP2K PLUS_U textbook antiferromagnet
- [ ] `bands-dos-si` | Band structure / DOS | pw | Si | NC-UPF | ε(k) / DOS | roadmap ○ post-processing; CP2K BAND_STRUCTURE / DOS
- [ ] `mpi-gpw-scaling-4rank` | MPI GPW energy scaling | infra | water | MOLOPT | wall ratio | industrial gate; may start `perf_exempt` until domain decomp lands
- [ ] `berry-polarization` | Berry-phase polarization | response | Si / BaTiO3 toy | — | P | roadmap ○; CP2K Berry phase; after robust k-path infra
- [ ] `wannier-centers` | Maximally localized Wannier | response | water / Si | — | centers | roadmap ○; CP2K Wannier; after bands foothold

## in_progress

(empty)

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
