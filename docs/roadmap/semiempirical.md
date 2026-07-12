# Semiempirical and Tight-Binding Methods

Legend and caveats: [Roadmap index](index.md).

Atomiverse ships one semiempirical method: **AM1**, on the NDDO framework
(`cpu::semiempirical::NDDOSolver` + `AM1Evaluator`). It computes energies,
Mulliken charges, and the density matrix — and it has **no gradients**. The MD
driver rejects it: binding `AM1Evaluator` to a force phase throws. Its main
production use today is as the base for AM1-BCC charge assignment (see
[Multiscale and Properties](multiscale.md)).

There is no DFTB and no extended tight binding (xTB/GFN-n) anywhere in the tree.

## 1. Semiempirical Quantum Methods

| Method                          | CPU | CUDA |
| ------------------------------- | :-: | :--: |
| AM1                             |  ✓  |   ○  |
| RM1                             |  ○  |   ○  |
| PM6                             |  ○  |   ○  |
| PM7                             |  ○  |   ○  |
| OM3                             |  ○  |   ○  |
| Semiempirical analytic gradient |  ○  |   ○  |
| Semiempirical analytic Hessian  |  ○  |   ○  |

The NDDO core (overlap, one-electron core Hamiltonian, two-electron integrals,
Fock build, density, SCF) is parameterization-agnostic in structure, but only
the AM1 parameter set is transcribed. MNDO/PM3/PM6 are parameter tables away in
principle — the missing gradient is the larger gap, since it is what blocks
semiempirical MD and optimization.

## 2. Density-Functional Tight Binding

| Method                    | CPU | CUDA |
| ------------------------- | :-: | :--: |
| SCC-DFTB / DFTB2          |  ○  |   ○  |
| DFTB3                     |  ○  |   ○  |
| Spin-polarized DFTB       |  ○  |   ○  |
| Noncollinear-spin DFTB    |  ○  |   ○  |
| DFTB+U                    |  ○  |   ○  |
| Dispersion-corrected DFTB |  ○  |   ○  |
| Spin-orbit DFTB           |  ○  |   ○  |
| Periodic DFTB             |  ○  |   ○  |
| DFTB analytic gradient    |  ○  |   ○  |
| DFTB analytic stress      |  ○  |   ○  |
| TD-DFTB                   |  ○  |   ○  |

## 3. Extended Tight-Binding Methods

| Method                    | CPU | CUDA |
| ------------------------- | :-: | :--: |
| GFN0-xTB                  |  ○  |   ○  |
| GFN2-xTB                  |  ○  |   ○  |
| Spin-polarized GFN-xTB    |  ○  |   ○  |
| IPEA-xTB                  |  ○  |   ○  |
| Periodic GFN-xTB          |  ○  |   ○  |
| GFN-xTB analytic gradient |  ○  |   ○  |
| GFN-xTB analytic Hessian  |  ○  |   ○  |
| GFN-FF                    |  ○  |   ○  |
| sTDA-xTB                  |  ○  |   ○  |
