# Hartree–Fock

Tutorial: [HF](../tutorial/ch02-rhf.md) · [open shell](../tutorial/ch03-open-shell.md) · [basis](../tutorial/ch04-basis-sets.md) · [SCF](../tutorial/ch18-scf-techniques.md).

Shares one SCF solver with [Kohn–Sham DFT](dft.md): the 2×2 product
{restricted, unrestricted} × {HF, KS}. Correlated post-HF lives in
[Post-Hartree–Fock](post-hf.md). Legend: [Roadmap index](index.md).

## 1. Molecular Hartree–Fock

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| RHF | ✓ | ✓ | CUDA: `cuda::Driver<scf::Phase>` (fock, mix, diagonalize, density, check); device Fock/density/CDIIS; STO-3G s/p (l ≤ 1) native ERI+forces. PySCF energy: H₂/STO-3G on device; H₂O/6-31G on CPU (6-31G is not a device ERI kernel) |
| UHF | ✓ | ○ | Default open-shell workhorse |
| ROHF | ✓ | ○ | Guest–Saunders; energy + analytic G + ROMP2 + FD Hessian |
| GHF | ✓ | ○ | Real two-component; energy + analytic G |
| Complex / spinor HF | ○ | ○ | |
| ΔSCF / MOM / IMOM | ✓ | ○ | `UHF::delta_scf` |
| RHF analytic gradient | ✓ | ✓ | CUDA same envelope as RHF energy |
| UHF / ROHF / GHF analytic gradient | ✓ | ○ | |
| RHF / UHF analytic Hessian | ✓ | ○ | CPHF (`RHF::hessian` / `UHF::hessian`); FD skeleton is not a mean-field method |
| ROHF Hessian | ✓ | ○ | FD of analytic forces (not full CPHF) |
| RHF static polarizability | ✓ | ○ | CPHF |
| Stability analysis | ✓ | ○ | Seeger–Pople |

### Gaps

- CUDA is a **restricted-HF engine**: no UHF/KS; forces outside l ≤ 1 throw —
  use `cpu::scf::RHF`.
- ROHF: SOSCF / level shift still throw. GHF: `MP2` / `FCI` / soscf throw.
- Post-HF and continuum solvent are not `MeanFieldMethod` members (`MP2` / `CCSD` / `WithSolvent`). SCF there is `run` / `forces` only.

## 2. Periodic Hartree–Fock

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| Γ-point periodic RHF | ✓ | ○ | Truncated lattice; `RHF.periodic`; short-range erfc gate |
| k-point / UHF periodic HF | ○ | ○ | Need Ewald for long-range 1/r |

PW/GPW product paths are semi-local DFT — see [DFT § Periodic](dft.md#4-periodic-kohnsham).

## 3. SCF infrastructure

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| DIIS | ✓ | ✓ | Device CDIIS |
| Linear damping | ✓ | ○ | |
| Level shift (Saunders–Hillier) | ✓ | ✓ | |
| SOSCF | ✓ | ○ | |
| Canonical orthogonalization | ✓ | ○ | |
| Hcore / MINAO guess | ✓ | ✓ / ○ | MINAO: H, C, N, O, F only; else Hcore |
| ASPC density predictor | ✓ | ○ | AIMD warm start |
| ⟨S²⟩ analysis | ✓ | ○ | |
| Spherical (pure) AO | ✓ | ○ | e.g. def2-SVP |
| RI-J / RI-JK | ✓ | ○ | def2-universal JKFIT; SCF only |
| CFMM-J | ✓ | ○ | `Coulomb::CFMM` |
| Native Rys ERI (l ≤ 3) | ✓ | ○ | s/p/d/f |
| Device ERI | — | ✓ | s/p, l ≤ 1 |
| ECP | ✓ | ○ | Semi-numerical forces |
| x2c-1e (energy) | ✓ | ○ | Forces throw |
| Diagonalization | ✓ | ✓ | LAPACK / cuSOLVER |

**Bases:** STO-3G, 6-31G, 6-31G\*, def2-SVP, def2-TZVP. Cartesian default;
spherical available where the set defines it.
