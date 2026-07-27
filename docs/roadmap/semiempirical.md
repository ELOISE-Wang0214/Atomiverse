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
**mio-style sketch** SK tables (not DFTB.org mio-1-1 bit-identical;
`reg_h2_dftb2`, `reg_h2o_dftb2_mio`) and production **3ob-3-1** organic CHNO
tables (`DftbSkSet::ThreeOb`, `reg_dftb2-3ob-organic`, offline DFTB+ pin).
**Closed-shell** nuclear forces on MioSketch are **analytic** (Hellmann–Feynman
+ Pulay + SCC γ + repulsive; `reg_dftb2-analytic-forces`); **open-shell** and
**ThreeOb** paths still use re-SCC FD (deepen: `dftb2-spin-analytic-forces`,
`dftb2-analytic-forces-threeob`). **Collinear spin-DFTB2** ships
(`multiplicity`, dual α/β, Mulliken m; `reg_spin-dftb2-o2`). **GFN0-xTB toy**
(`XtbGfn0`, `reg_h2o_gfn0`) remains a sketch. Non-SCC DFTB0, DFTB3,
D3/PBC/k-mesh/stress/TD depth, and GFN1/2 remain ○ (several on unmerged
`wf/*` tips — see queue `## in_progress`). Loop queue:
`docs/developer/tb-method-queue.md`. Family oracles: OpenMOPAC (NDDO) / DFTB+
(DFTB) / xtb (GFN) — see `tools/tb_oracle/README.md`.

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
| Non-SCC DFTB / DFTB0      |  ○  |   ○  |   ○  |
| SCC-DFTB / DFTB2          |  ✓  |   ○  |   ○  |
| DFTB3                     |  ○  |   ○  |   ○  |
| Spin-polarized DFTB       |  ✓  |   ○  |   ○  |
| Noncollinear-spin DFTB    |  ○  |   ○  |   ○  |
| DFTB+U                    |  ○  |   ○  |   ○  |
| Dispersion-corrected DFTB |  ○  |   ○  |   ○  |
| Spin-orbit DFTB           |  ○  |   ○  |   ○  |
| Periodic DFTB             |  ○  |   ○  |   ○  |
| DFTB analytic gradient    |  ✓  |   ○  |   ○  |
| DFTB analytic stress      |  ○  |   ○  |   ○  |
| TD-DFTB                   |  ○  |   ○  |   ○  |

DFTB2 is multi-element H/C/N/O with static SK tables (`dftb_mio_tables.h` sketch
+ `dftb_3ob_tables.h` 3ob-3-1 organics). Collinear spin (`multiplicity`) and
**restricted** analytic nuclear forces (HF+Pulay+γ+rep) ship on MioSketch;
open-shell and ThreeOb forces remain re-SCC FD. Offline **DFTB+** absolute
energy identity on mio H–H is an incomplete tip (`dftb2-dftbplus-energy-pin`);
3ob CH4 has a DFTB+ pin under pure SCC-DFTB2. Analytic **stress**, **k-mesh**
beyond Γ, **DFTB3** forces, and **TD-DFTB** remain ○. Oracle for all DFTB
rows: **DFTB+** (never OpenMOPAC). Ladder: DFTB0 → DFTB2 depth → DFTB3 →
spin product depth → …

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
