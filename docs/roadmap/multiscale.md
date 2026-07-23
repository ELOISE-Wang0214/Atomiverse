# Multiscale Methods and Properties

Legend and caveats: [Roadmap index](index.md).

Everything that couples two levels of theory, or extracts a chemical property
from a converged wavefunction. Charge schemes still anchor the product surface
(Mulliken observer; AM1-BCC via [AM1](semiempirical.md) + molrs typing). Beyond
that, several *environment* hooks now exist on the SCF mean-field methods:

- **C-PCM / COSMO** electrostatic free energy in the Fock build
  (`reg_h2o_rhf_631gs_cpcm`) — non-electrostatic terms absent.
- **Fixed point-charge embedding** + **uniform external electric field**, with
  analytic forces (`reg_embed_h2o_grad`, `reg_field_h2o_grad`).
- **ECP** (semi-numerical nuclear forces) and **x2c-1e** (energy only; forces
  throw).

**Textbook QM/MM** ships as `QmmmMechanical` (DFTB2 QM + LJ),
`QmmmElectrostatic` (RHF + MM point charges + LJ, optional **H link atoms**),
and `QmmmPolarizable` (RHF + induced MM dipoles, SCC). Still no ONIOM /
subtractive product path, no SMD/GB, and no picture-change / DKH / SOC stack.
Plane-wave UPF and GPW GTH remain frozen-core pseudopotentials, separate from
Gaussian ECPs.

## 1. Atomic-Charge Methods

| Method                                 | CPU | CUDA | Kokkos |
| -------------------------------------- | :-: | :--: | :---: |
| Mulliken charges                       |  ✓  |   ✓  |   ○  |
| Löwdin charges                         |  ○  |   ○  |   ○  |
| Natural population analysis            |  ○  |   ○  |   ○  |
| Hirshfeld charges                      |  ○  |   ○  |   ○  |
| Iterative Hirshfeld charges            |  ○  |   ○  |   ○  |
| CM5 charges                            |  ○  |   ○  |   ○  |
| Bader charges                          |  ○  |   ○  |   ○  |
| DDEC charges                           |  ○  |   ○  |   ○  |
| Merz–Kollman ESP charges               |  ○  |   ○  |   ○  |
| CHELPG charges                         |  ○  |   ○  |   ○  |
| RESP charges                           |  ○  |   ○  |   ○  |
| RESP2 charges                          |  ○  |   ○  |   ○  |
| AM1-BCC charges                        |  ✓  |   ○  |   ○  |
| Gasteiger charges                      |  ○  |   ○  |   ○  |
| Electronegativity-equalization charges |  ○  |   ○  |   ○  |
| QEq charges                            |  ○  |   ○  |   ○  |
| EQeq charges                           |  ○  |   ○  |   ○  |

Löwdin is ○ deliberately: the Löwdin **orthogonalization** (S^(−1/2)) is used
inside the SCF, but no Löwdin *population analysis* is exposed. Charges are
Observer / method outputs written to Frame, never Drivers — there is no
"charge driver".

## 2. Implicit-Solvent Methods

| Method                             | CPU | CUDA | Kokkos |
| ---------------------------------- | :-: | :--: | :---: |
| Polarizable continuum model, PCM   |  ○  |   ○  |   ○  |
| Integral-equation-formalism PCM    |  ○  |   ○  |   ○  |
| Conductor-like PCM                 |  ✓  |   ○  |   ○  |
| COSMO                              |  ✓  |   ○  |   ○  |
| SMD                                |  ○  |   ○  |   ○  |
| SM12                               |  ○  |   ○  |   ○  |
| Generalized Born                   |  ○  |   ○  |   ○  |
| GBSA                               |  ○  |   ○  |   ○  |
| ALPB                               |  ○  |   ○  |   ○  |
| Poisson–Boltzmann                  |  ○  |   ○  |   ○  |
| PBSA                               |  ○  |   ○  |   ○  |
| Nonequilibrium PCM                 |  ○  |   ○  |   ○  |
| State-specific excited-state PCM   |  ○  |   ○  |   ○  |
| Implicit-solvent analytic gradient |  ○  |   ○  |   ○  |
| Implicit-solvent analytic Hessian  |  ○  |   ○  |   ○  |

C-PCM / COSMO share one conductor-like cavity implementation (`internal::Pcm`,
x=0 vs x=0.5). Energy contribution is validated; surface-switching / SWIG and
analytic PCM gradients are out of scope for the current path.

## 3. QM/MM Methods

| Method                           | CPU | CUDA | Kokkos |
| -------------------------------- | :-: | :--: | :---: |
| Mechanical embedding QM/MM       |  ✓  |   ○  |   ○  |
| Electrostatic embedding QM/MM    |  ✓  |   ○  |   ○  |
| Polarizable embedding QM/MM      |  ✓  |   ○  |   ○  |
| Additive QM/MM                   |  ○  |   ○  |   ○  |
| Subtractive QM/MM                |  ○  |   ○  |   ○  |
| ONIOM two-layer method           |  ○  |   ○  |   ○  |
| ONIOM multilayer method          |  ○  |   ○  |   ○  |
| Link-atom QM/MM                  |  ✓  |   ○  |   ○  |
| Localized-orbital boundary QM/MM |  ○  |   ○  |   ○  |
| Pseudobond QM/MM                 |  ○  |   ○  |   ○  |
| Periodic QM/MM                   |  ○  |   ○  |   ○  |
| Long-range PME QM/MM             |  ○  |   ○  |   ○  |
| Adaptive QM/MM                   |  ○  |   ○  |   ○  |
| Multi-region QM/MM               |  ○  |   ○  |   ○  |
| QM/MM geometry optimization      |  ○  |   ○  |   ○  |
| QM/MM molecular dynamics         |  ○  |   ○  |   ○  |
| Excited-state QM/MM              |  ○  |   ○  |   ○  |

Fixed point-charge embedding on the SCF method (external MM charges in the
one-electron Hamiltonian, with analytic forces) remains a lower-level hook.
The product QM/MM components are:

- `QmmmMechanical` — DFTB2 on QM tags + LJ environment (`reg_qmmm_h2_ar`)
- `QmmmElectrostatic` — RHF QM + MM charges + LJ, with optional H link atoms
  (`reg_qmmm_es_h2`, `reg_qmmm_link_ch3oh`)
- `QmmmPolarizable` — RHF + isotropic induced MM dipoles (`reg_qmmm_pol_h2o`)

No automatic region partitioner, no ONIOM, no PME long-range QM/MM.

## 4. Quantum Embedding Methods

| Method                            | CPU | CUDA | Kokkos |
| --------------------------------- | :-: | :--: | :---: |
| Fixed point-charge embedding      |  ✓  |   ○  |   ○  |
| Distributed-multipole embedding   |  ○  |   ○  |   ○  |
| Polarizable embedding             |  ○  |   ○  |   ○  |
| Effective fragment potential      |  ○  |   ○  |   ○  |
| Frozen-density embedding          |  ○  |   ○  |   ○  |
| Projection-based embedding        |  ○  |   ○  |   ○  |
| Density-matrix embedding theory   |  ○  |   ○  |   ○  |
| Subsystem DFT                     |  ○  |   ○  |   ○  |
| Fragment molecular orbital method |  ○  |   ○  |   ○  |

## 5. Relativistic Methods

| Method                                        | CPU | CUDA | Kokkos |
| --------------------------------------------- | :-: | :--: | :---: |
| Scalar-relativistic effective core potentials |  ✓  |   ○  |   ○  |
| Spin-orbit effective core potentials          |  ○  |   ○  |   ○  |
| Douglas–Kroll–Hess second order               |  ○  |   ○  |   ○  |
| Scalar ZORA                                   |  ○  |   ○  |   ○  |
| Spin-orbit ZORA                               |  ○  |   ○  |   ○  |
| Scalar X2C                                    |  ✓  |   ○  |   ○  |
| Spin-orbit X2C                                |  ○  |   ○  |   ○  |
| Exact two-component molecular mean field      |  ○  |   ○  |   ○  |
| Four-component Dirac–Hartree–Fock             |  ○  |   ○  |   ○  |
| Four-component Dirac–Kohn–Sham                |  ○  |   ○  |   ○  |
| Breit interaction                             |  ○  |   ○  |   ○  |
| Spin-orbit mean-field approximation           |  ○  |   ○  |   ○  |
| Relativistic analytic gradient                |  ○  |   ○  |   ○  |

Gaussian ECPs are available with semi-numerical forces (`reg_ecp_i2_grad`,
`reg_ecp_heavy_atom`). x2c-1e replaces the one-electron Hamiltonian (energy
path validated on Ar); analytic forces with `x2c=true` throw. Stacking x2c
on an ECP atom is refused. Periodic calculations still use UPF/GTH
pseudopotentials, separate from the Gaussian ECP tables.

## 6. Intermolecular Interaction Methods

### 6.1 Symmetry-adapted perturbation theory

| Method                 | CPU | CUDA | Kokkos |
| ---------------------- | :-: | :--: | :---: |
| SAPT0                  |  ○  |   ○  |   ○  |
| Scaled SAPT0           |  ○  |   ○  |   ○  |
| SAPT2                  |  ○  |   ○  |   ○  |
| SAPT2+                 |  ○  |   ○  |   ○  |
| SAPT2+(3)              |  ○  |   ○  |   ○  |
| SAPT2+3                |  ○  |   ○  |   ○  |
| SAPT(DFT)              |  ○  |   ○  |   ○  |
| XSAPT                  |  ○  |   ○  |   ○  |
| XSAPT+MBD              |  ○  |   ○  |   ○  |
| SAPT analytic gradient |  ○  |   ○  |   ○  |

### 6.2 Energy-decomposition methods

| Method                                   | CPU | CUDA | Kokkos |
| ---------------------------------------- | :-: | :--: | :---: |
| Morokuma energy decomposition            |  ○  |   ○  |   ○  |
| ALMO energy-decomposition analysis       |  ○  |   ○  |   ○  |
| Natural energy-decomposition analysis    |  ○  |   ○  |   ○  |
| Fragment molecular orbital decomposition |  ○  |   ○  |   ○  |
| Many-body expansion                      |  ○  |   ○  |   ○  |
| Counterpoise correction                  |  ✓  |   ○  |   ○  |

## 7. Composite and Multiscale Methods

| Method                                       | CPU | CUDA | Kokkos |
| -------------------------------------------- | :-: | :--: | :---: |
| Delta-learning QM correction                 |  ○  |   ○  |   ○  |
| ML/MM hybrid potential                       |  ○  |   ○  |   ○  |
| ML/QM delta correction                       |  ○  |   ○  |   ○  |
| Fragment-based quantum calculation           |  ○  |   ○  |   ○  |
| Many-body expansion potential                |  ○  |   ○  |   ○  |
| Adaptive-resolution simulation               |  ○  |   ○  |   ○  |
| Concurrent atomistic–continuum coupling      |  ○  |   ○  |   ○  |
| Quantum embedding molecular dynamics         |  ○  |   ○  |   ○  |
| Multiple-time-scale mixed-potential dynamics |  ○  |   ○  |   ○  |

Blocked upstream where product ML potentials are still ○ (see
[Force Fields](force-fields.md#6-machine-learning-potentials)). Textbook
mechanical and electrostatic QM/MM skeletons ship; full product QM/MM does not.
