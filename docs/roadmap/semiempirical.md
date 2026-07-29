# Semiempirical and Tight-Binding Methods

Legend and caveats: [Roadmap index](index.md).

Atomiverse ships two **concrete NDDO methods**: **AM1** and **RM1**
(`cpu::semiempirical::AM1` / `RM1`; evaluators `AM1Evaluator` / `RM1Evaluator`
of kind `cpu::Semiempirical`). The shared NDDO engine is **internal**
(`NDDOSolver` is not a public bag). Methods compute energies, Mulliken
charges, the density matrix, and **nuclear forces** via fixed-density central
finite difference of the SCF energy (product path; matches re-SCF FD at SCF
stationarity). Closed-form multipole / integral derivatives remain ○ (see
queue `nddo-closed-form-forces`). Evaluators write forces for Evaluate / MD /
Minimize by default. Launch unwraps an isolated fragment for molecular AM1
under PBC wrap (same pattern as Gaussian SCF). Regressions: `am1_water_force`,
`am1_water_minimize`, `am1_water_bomd`, `rm1_water_energy`. Tutorial:
[Chapter 17](../tutorial/ch17-semiempirical-am1.md) + `examples/15_am1.cpp`.
Heat of formation is MOPAC-style. Hard OpenMOPAC energy pins and PM6/PM7 remain
open. AM1-BCC charge assignment still uses the energy/charge path (see
[Multiscale and Properties](multiscale.md)).

**SCC-DFTB2** multi-element H/C/N/O (`DFTB2`) ships with static in-tree
**mio-style sketch** SK tables (`reg_h2_dftb2`, `reg_h2o_dftb2_mio`),
**mio-1-1 H–H** radial grids for the absolute DFTB+ energy pin
(`dftb_mio_hh_sk.h`, `reg_dftb2-dftbplus-energy-pin`), and production
**3ob-3-1** organic CHNO tables (`DftbSkSet::ThreeOb`, `reg_dftb2-3ob-organic`).
**Nuclear forces** are analytic on MioSketch for closed shell **and collinear
open shell** (Hellmann–Feynman + Pulay + SCC γ + repulsive;
`reg_dftb2-analytic-forces`, `reg_dftb2-spin-analytic-forces`) and on the
**ThreeOb** organic path (`reg_dftb2-analytic-forces-threeob`). **Collinear
spin-DFTB2** ships (`multiplicity`, dual α/β, Mulliken m; `reg_spin-dftb2-o2`).
**Γ-point periodic** SCC-DFTB2 lattice sums (H⁰/S, screened γ, repulsive) ship
on MioSketch (`reg_dftb2-periodic-gamma`) with **lattice-image-consistent
nuclear forces** under PBC (`reg_dftb2-periodic-forces`). **Monkhorst–Pack
k-mesh** beyond Γ ships for closed-shell MioSketch **energy**
(`k_mesh` / `k_time_reversal` members; `reg_dftb2-kpoints`); multi-k analytic
forces/stress remain ○ (`dftb2-multi-k-forces`). Opt-in **Grimme D3(BJ)**
dispersion ships as an additive **molecular** correction with its own analytic
gradient (`DFTB2::d3`, `reg_dftb2-dispersion-d3`); a periodic D3 lattice sum
remains ○ (`dftb2-dispersion-d3-periodic`). Non-SCC **DFTB0** ships
(`reg_dftb0-h2`). **GFN0-xTB toy** (`XtbGfn0`, `reg_h2o_gfn0`) remains a
sketch. **DFTB3** energy ships (`reg_dftb3-water-mio`) after its third-order
potential was corrected to be the derivative of its own energy; DFTB3 analytic
forces remain ○ (`dftb3-analytic-forces`). **Analytic stress** ships as the
closed-form pair virial on restricted MioSketch (`reg_dftb2-analytic-stress`)
**and on the ThreeOb 3ob-3-1 path**
(`reg_dftb2-analytic-stress-threeob`, 2026-07-29 — the multi-SCC strain-FD
product path for ThreeOb is **deleted**); collinear open-shell stress still uses
multi-SCC strain-FD (`dftb2-spin-analytic-stress`). Wave-3 TB staff closed 2026-07-29; GapCensus
wave-4 (2026-07-29) refilled Phase A deepens — TD depth and GFN1/2 remain ○.
Loop queue: `docs/developer/tb-method-queue.md`. Family oracles: OpenMOPAC
(NDDO) / DFTB+ (DFTB) / xtb (GFN) — see `tools/tb_oracle/README.md`.

## 1. Semiempirical Quantum Methods

| Method                          | CPU | CUDA | Kokkos |
| ------------------------------- | :-: | :--: | :---: |
| AM1                             |  ✓  |   ○  |   ○  |
| RM1                             |  ✓  |   ○  |   ○  |
| PM3                             |  ○  |   ○  |   ○  |
| PM6                             |  ○  |   ○  |   ○  |
| PM7                             |  ○  |   ○  |   ○  |
| OM3                             |  ○  |   ○  |   ○  |
| Semiempirical nuclear forces (product path) |  ✓  |   ○  |   ○  |
| Semiempirical closed-form integral-derivative forces |  ○  |   ○  |   ○  |
| Semiempirical analytic Hessian  |  ○  |   ○  |   ○  |
| AM1 heat of formation (MOPAC-style ΔHf) |  ✓  |   ○  |   ○  |

The NDDO core (overlap, one-electron core Hamiltonian, two-electron integrals,
Fock build, density, SCF, forces) is parameterization-agnostic in structure.
**AM1** and **RM1** (H/C/N/O/P/S; Rocha et al. 2006 / MOPAC tables) are
transcribed. PM3/PM6/PM7 remain tables away (oracle = **OpenMOPAC**). OM3 is
not in the OpenMOPAC primary surface (blocked in the TB queue until a
family-locked pin exists). Product-path forces are fixed-density FD; closed-form
multipole derivatives stay ○. Analytic Hessian remains finite-difference via
`cpu::vib::Hessian` over SE forces (VIB off under SE-only loop cmake).

## 2. Density-Functional Tight Binding

| Method                    | CPU | CUDA | Kokkos |
| ------------------------- | :-: | :--: | :---: |
| Non-SCC DFTB / DFTB0      |  ✓  |   ○  |   ○  |
| SCC-DFTB / DFTB2          |  ✓  |   ○  |   ○  |
| DFTB3                     |  ✓  |   ○  |   ○  |
| Spin-polarized DFTB       |  ✓  |   ○  |   ○  |
| Noncollinear-spin DFTB    |  ○  |   ○  |   ○  |
| DFTB+U                    |  ○  |   ○  |   ○  |
| Dispersion-corrected DFTB |  ✓  |   ○  |   ○  |
| Spin-orbit DFTB           |  ○  |   ○  |   ○  |
| Periodic DFTB             |  ✓  |   ○  |   ○  |
| DFTB analytic gradient    |  ✓  |   ○  |   ○  |
| DFTB analytic stress      |  ✓  |   ○  |   ○  |
| TD-DFTB                   |  ✓  |   ○  |   ○  |

DFTB2 is multi-element H/C/N/O with static SK tables (`dftb_mio_tables.h` sketch
+ `dftb_mio_hh_sk.h` mio-1-1 H–H + `dftb_3ob_tables.h` 3ob-3-1 organics).
Collinear spin (`multiplicity`) and **analytic nuclear forces**
(HF+Pulay+γ+rep; **open shell included** on MioSketch; **ThreeOb** closed-form
F ships) are product paths. **Analytic stress** ships as the closed-form pair
virial on restricted MioSketch (`reg_dftb2-analytic-stress`) **and on ThreeOb**
(`reg_dftb2-analytic-stress-threeob`; all six Voigt components vs strain-FD at
5e-13 eV/Å³, plus the isotropic tr(σ)/3 = −dE/dΩ identity); only collinear
open shell still falls back to multi-SCC affine strain-FD. Offline **DFTB+** absolute energy
identity on mio-1-1 H–H ships (`reg_dftb2-dftbplus-energy-pin`); 3ob CH4 has a
DFTB+ pin under pure SCC-DFTB2. **Γ-point PBC** lattice-sum energy **and**
lattice-image nuclear forces ship; **multi-k energy** (Monkhorst–Pack) ships
for closed-shell MioSketch. **TD-DFTB** ships as a linear-response Casida
foothold on the restricted MioSketch ground state
(`DFTB2::td_casida`; Niehaus et al., PRB **63**, 085108 (2001);
`reg_td-dftb-h2`): **energies only** (no oscillator strengths, no excited-state
gradients), singlet + triplet, molecular Γ only. It has **no absolute oracle** —
no `dftb+` binary on the pin host and no published H₂/mio-1-1 excitation energy
exists, so it is gated on the λ→0 single-particle limit, Ω symmetry, the
transition-charge sum rule and analytic two-level identities; see its
`REFERENCE.md` `Gaps` (which also records a **pre-existing γ_AB units defect**
shared with the DFTB3 ground state — see `.claude/notes/`). Wave 4 also landed
**multi-k analytic forces + stress**, **DFTB3 analytic forces**, **periodic
lattice-sum D3**, and **ThreeOb / open-shell analytic stress**, so no
strain-FD product path is left on DFTB2. Still ○: **TD-DFTB oscillator
strengths / excited-state gradients**, **long-range-corrected TD-DFTB**,
**multi-k open shell**, **ThreeOb multi-k**, and **DFTB3 beyond molecular
min-image** (no PBC / k-mesh / stress / open shell / ThreeOb / D3 on the DFTB3
path). Oracle for all DFTB rows: **DFTB+**
(never OpenMOPAC). Ladder: DFTB0 → DFTB2 depth → DFTB3 F → TD → NDDO PM →
GFN late.

## 3. Extended Tight-Binding Methods

| Method                    | CPU | CUDA | Kokkos |
| ------------------------- | :-: | :--: | :---: |
| GFN0-xTB                  |  ✓  |   ○  |   ○  |
| GFN1-xTB                  |  ○  |   ○  |   ○  |
| GFN2-xTB                  |  ○  |   ○  |   ○  |
| Spin-polarized GFN-xTB    |  ○  |   ○  |   ○  |
| IPEA-xTB                  |  ○  |   ○  |   ○  |
| Periodic GFN-xTB          |  ○  |   ○  |   ○  |
| GFN-xTB analytic gradient |  ○  |   ○  |   ○  |
| GFN-xTB analytic Hessian  |  ○  |   ○  |   ○  |
| GFN-FF                    |  ○  |   ○  |   ○  |
| sTDA-xTB                  |  ○  |   ○  |   ○  |

GFN0 is a **toy** H/C/O skeleton (`XtbGfn0`) — not production GFN0 parameters
or xtb parity. GFN-xTB is **not** a DFTB3 successor; schedule after the DFTB
ladder (oracle = **xtb**). GFN-FF is classical force-field territory (scope-out
of the electronic SE/TB capability loop).
