# Semiempirical and Tight-Binding Methods

Tutorial: [Semiempirical models](../tutorial/ch17-semiempirical-am1.md).

CPU-only domain. Oracles: **OpenMOPAC** (NDDO), **DFTB+** (DFTB), **xtb
6.7.1** (GFN/IPEA). Legend: [Roadmap index](index.md).

## 1. NDDO

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| AM1 | ✓ | ○ | `cpu::semiempirical::AM1`; OpenMOPAC ΔHf pins |
| RM1 | ✓ | ○ | Same NDDO path; MOPAC `calpar` free-atom repair |
| PM3 / PM6 / PM7 | ○ | ○ | Tables away |
| OM3 | ○ | ○ | **Blocked** — no OpenMOPAC primary surface |
| Nuclear forces (product) | ✓ | ○ | Fixed-density FD of SCF energy |
| Closed-form integral-derivative forces | ○ | ○ | |
| FD Hessian via `vib::Hessian` | ✓ | ○ | |
| Heat of formation (MOPAC-style) | ✓ | ○ | |
| NDDO UHF / COSMO / Γ-PBC | ○ | ○ | Molecular AM1 under wrap unwraps an isolated fragment |

## 2. DFTB

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| DFTB0 | ✓ | ○ | Forces still 6N FD |
| DFTB2 (SCC) | ✓ | ○ | mio-1-1 + 3ob-3-1 (H/C/N/O) |
| DFTB3 | ✓ | ○ | Energy + closed-form F + Γ-PBC F; no k / stress / spin / 3ob / D3 |
| Collinear spin-DFTB2 | ✓ | ○ | MioSketch; ThreeOb open shell throws |
| DFTB+U | ✓ | ○ | Analytic F still ○ |
| D3 dispersion (mol + Γ-PBC) | ✓ | ○ | DFTB2; D4 ○ |
| Γ-point / multi-k periodic DFTB2 | ✓ | ○ | Analytic F/σ on MioSketch; ThreeOb is molecular only |
| DFTB2 analytic gradient / stress | ✓ | ○ | Pair virial; all product paths |
| TD-DFTB | ✓ | ○ | See [Excited States](excited-states.md) |
| LC-DFTB / Fermi filling / GBSA | ○ | ○ | |

## 3. Extended tight binding (xTB)

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| GFN0 / GFN1 / GFN2 | ✓ | ○ | H/C/N/O; analytic G; FD Hessian |
| GFN2 AES (quadrupole) | ✓ | ○ | In SCC Fock |
| IPEA-xTB | ✓ | ○ | Energy + G + FD Hessian + ΔSCC VIP/VEA |
| Periodic / spin / ALPB / Fermi / sTDA-xTB | ○ | ○ | |
| GFN-FF | ○ | ○ | Classical FF scope-out |

### Gaps

- Elements beyond H/C/N/O in shipped SK / xTB tables.
- Analytic (closed-form) Hessians anywhere in this chapter — all FD via
  `cpu::vib::Hessian`.
- No CUDA path.
