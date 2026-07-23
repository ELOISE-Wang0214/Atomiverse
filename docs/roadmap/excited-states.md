# Excited States and Many-Body Methods

Legend and caveats: [Roadmap index](index.md).

**Ground-state engine with a closed-shell TDHF foothold.** `RHF::tda` (CIS /
TDA) and `RHF::tdhf` (full RPA) return excitation energies and length-gauge
transition dipoles on a converged RHF reference
(`regressions/scf/ch2o_tdhf`). **`RKS::tda`** is the CIS@RKS skeleton on a
converged KS reference; LSDA adds a **Dirac ALDA diagonal** on dense A
(`reg_h2_alda_vs_cis`). Full grid MO-product ALDA and hybrid TDDFT remain ○.
**FD TDA forces** and **NAC strength** ship (`RHF::tda_forces` /
`tda_nac_strength`); analytic ES gradients and true NAC vectors are still ○.

The same orbital-Hessian / A,B machinery also feeds Seeger–Pople stability
analysis on the ground-state HF path (see [Hartree–Fock](hartree-fock.md)).

## 1. SCF and Linear-Response Methods

| Method                             | CPU | CUDA | Kokkos |
| ---------------------------------- | :-: | :--: | :---: |
| Delta-SCF excited states           |  ○  |   ○  |   ○  |
| MOM excited states                 |  ○  |   ○  |   ○  |
| IMOM excited states                |  ○  |   ○  |   ○  |
| ROKS excited states                |  ○  |   ○  |   ○  |
| CIS excitation energies            |  ✓  |   ○  |   ○  |
| CIS(D)                             |  ○  |   ○  |   ○  |
| SCS-CIS(D)                         |  ○  |   ○  |   ○  |
| SOS-CIS(D)                         |  ○  |   ○  |   ○  |
| Time-dependent Hartree–Fock, TDHF  |  ✓  |   ○  |   ○  |
| Tamm–Dancoff TDDFT                 |  ○  |   ○  |   ○  |
| Full linear-response TDDFT         |  ○  |   ○  |   ○  |
| Spin-flip TDA                      |  ○  |   ○  |   ○  |
| Spin-flip TDDFT                    |  ○  |   ○  |   ○  |
| Restricted excitation-window TDDFT |  ○  |   ○  |   ○  |
| Core–valence-separation TDDFT      |  ○  |   ○  |   ○  |
| Real-time TDDFT                    |  ○  |   ○  |   ○  |
| TDDFT analytic gradient            |  ○  |   ○  |   ○  |
| TDDFT nonadiabatic coupling        |  ○  |   ○  |   ○  |

## 2. Algebraic-Diagrammatic Construction

| Method                | CPU | CUDA | Kokkos |
| --------------------- | :-: | :--: | :---: |
| ADC(1)                |  ○  |   ○  |   ○  |
| ADC(2)                |  ○  |   ○  |   ○  |
| ADC(2)-x              |  ○  |   ○  |   ○  |
| ADC(3)                |  ○  |   ○  |   ○  |
| CVS-ADC(2)            |  ○  |   ○  |   ○  |
| CVS-ADC(2)-x          |  ○  |   ○  |   ○  |
| IP-ADC(2)             |  ○  |   ○  |   ○  |
| EA-ADC(2)             |  ○  |   ○  |   ○  |
| ADC analytic gradient |  ○  |   ○  |   ○  |

## 3. Equation-of-Motion Coupled Cluster

| Method                     | CPU | CUDA | Kokkos |
| -------------------------- | :-: | :--: | :---: |
| EOM-EE-CCSD                |  ○  |   ○  |   ○  |
| EOM-IP-CCSD                |  ○  |   ○  |   ○  |
| EOM-EA-CCSD                |  ○  |   ○  |   ○  |
| EOM-SF-CCSD                |  ○  |   ○  |   ○  |
| EOM-DIP-CCSD               |  ○  |   ○  |   ○  |
| EOM-DEA-CCSD               |  ○  |   ○  |   ○  |
| CVS-EOM-CCSD               |  ○  |   ○  |   ○  |
| EOM-CC2                    |  ○  |   ○  |   ○  |
| EOM-MP2                    |  ○  |   ○  |   ○  |
| STEOM-CCSD                 |  ○  |   ○  |   ○  |
| EOM-CCSD analytic gradient |  ○  |   ○  |   ○  |

Blocked upstream: there is no coupled cluster to build EOM on
(see [Post-Hartree–Fock](post-hf.md#3-coupled-cluster)).

## 4. Simplified Excited-State Methods

| Method                       | CPU | CUDA | Kokkos |
| ---------------------------- | :-: | :--: | :---: |
| Simplified TDA, sTDA         |  ○  |   ○  |   ○  |
| Simplified TDDFT, sTDDFT     |  ○  |   ○  |   ○  |
| sTDA-xTB                     |  ○  |   ○  |   ○  |
| TD-DFTB                      |  ○  |   ○  |   ○  |
| Long-range-corrected TD-DFTB |  ○  |   ○  |   ○  |

Blocked upstream: no DFTB, no xTB
(see [Semiempirical and Tight Binding](semiempirical.md)).

## 5. Green's-Function and Many-Body Methods

| Method                                        | CPU | CUDA | Kokkos |
| --------------------------------------------- | :-: | :--: | :---: |
| Second-order Green's function, GF2            |  ○  |   ○  |   ○  |
| Self-consistent GF2                           |  ○  |   ○  |   ○  |
| G0W0                                          |  ○  |   ○  |   ○  |
| Eigenvalue-self-consistent GW                 |  ○  |   ○  |   ○  |
| Quasiparticle self-consistent GW              |  ○  |   ○  |   ○  |
| Fully self-consistent GW                      |  ○  |   ○  |   ○  |
| Bethe–Salpeter equation, BSE                  |  ○  |   ○  |   ○  |
| Direct random-phase approximation             |  ○  |   ○  |   ○  |
| RPA with exchange                             |  ○  |   ○  |   ○  |
| SOSEX                                         |  ○  |   ○  |   ○  |
| Renormalized second-order perturbation theory |  ○  |   ○  |   ○  |
| ACFDT-RPA                                     |  ○  |   ○  |   ○  |
