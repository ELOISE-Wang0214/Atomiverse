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
| Hückel molecular orbital theory|  ⊘  |   ⊘  |
| Extended Hückel theory|  ⊘  |   ⊘  |
| CNDO/2|  ⊘  |   ⊘  |
| INDO|  ⊘  |   ⊘  |
| MINDO/3|  ⊘  |   ⊘  |
| MNDO|  ⊘  |   ⊘  |
| AM1                             |  ✓  |   ○  |
| RM1                             |  ○  |   ○  |
| PM3|  ⊘  |   ⊘  |
| PM6                             |  ○  |   ○  |
| PM7                             |  ○  |   ○  |
| OM1|  ⊘  |   ⊘  |
| OM2|  ⊘  |   ⊘  |
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
| Non-self-consistent DFTB0|  ⊘  |   ⊘  |
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
| GFN1-xTB|  ⊘  |   ⊘  |
| GFN2-xTB                  |  ○  |   ○  |
| Spin-polarized GFN-xTB    |  ○  |   ○  |
| IPEA-xTB                  |  ○  |   ○  |
| Periodic GFN-xTB          |  ○  |   ○  |
| GFN-xTB analytic gradient |  ○  |   ○  |
| GFN-xTB analytic Hessian  |  ○  |   ○  |
| GFN-FF                    |  ○  |   ○  |
| sTDA-xTB                  |  ○  |   ○  |

## Superseded (⊘)

The pre-NDDO approximations are all dominated by the NDDO/AM1 core that is
already implemented — they neglect strictly more two-electron integrals for no
gain in cost that matters today.

| Method | Dominated by |
| ------ | ------------ |
| Hückel, Extended Hückel | Any self-consistent method; these carry no two-electron terms at all. |
| CNDO/2, INDO, MINDO/3 | NDDO — retains the differential overlap these discard. |
| MNDO | AM1 (implemented) — AM1 *is* MNDO plus the corrected core repulsion. |
| PM3 | PM6/PM7 — the same NDDO Hamiltonian, better parameterized and far broader element coverage. |
| OM1, OM2 | OM3 — the same orthogonalization-corrected family, latest parameterization. |
| Non-self-consistent DFTB0 | SCC-DFTB (DFTB2) — charge self-consistency is what makes DFTB usable for anything polar. |
| GFN1-xTB | GFN2-xTB — same cost tier, better Hamiltonian (multipole electrostatics, D4). |

PM6, PM7, RM1, OM3, DFTB2/DFTB3 and GFN0/GFN2-xTB stay ○: each is the *top* of
its family, and each occupies a distinct cost/accuracy niche.
