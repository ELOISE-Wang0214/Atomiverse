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

**SCC-DFTB2** multi-element H/C/N/O (`DFTB2`, static in-tree **mio-style sketch**
SK tables — not DFTB.org mio-1-1 bit-identical; `reg_h2_dftb2`,
`reg_h2o_dftb2_mio`) and **GFN0-xTB toy** (`XtbGfn0`, `reg_h2o_gfn0`) ship as TB
families beside NDDO. DFTB2 forces are re-SCC central FD of E (`with_forces`),
not analytic. Non-SCC DFTB0, production 3ob/mio pins, DFTB3, and GFN1/2 remain ○.
Loop queue: `docs/developer/tb-method-queue.md`. Family oracles:
OpenMOPAC (NDDO) / DFTB+ (DFTB) / xtb (GFN) — see `tools/tb_oracle/README.md`.

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
| Spin-polarized DFTB       |  ○  |   ○  |   ○  |
| Noncollinear-spin DFTB    |  ○  |   ○  |   ○  |
| DFTB+U                    |  ○  |   ○  |   ○  |
| Dispersion-corrected DFTB |  ○  |   ○  |   ○  |
| Spin-orbit DFTB           |  ○  |   ○  |   ○  |
| Periodic DFTB             |  ○  |   ○  |   ○  |
| DFTB analytic gradient    |  ○  |   ○  |   ○  |
| DFTB analytic stress      |  ○  |   ○  |   ○  |
| TD-DFTB                   |  ○  |   ○  |   ○  |

DFTB2 is multi-element H/C/N/O with static SK tables (`basis_data/dftb_mio_tables.h`).
Energy foothold only (no offline **DFTB+** pin yet). Forces exist via re-SCC
central FD of energy (`with_forces`) — the analytic-gradient row stays ○ until
Hellmann–Feynman / SK derivatives ship. Oracle for all DFTB rows: **DFTB+**
(never OpenMOPAC). Ladder: DFTB0 → DFTB2 depth → DFTB3 → spin → …

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
