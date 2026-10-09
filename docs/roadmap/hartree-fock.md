# Hartree–Fock

Tutorial: [HF](../tutorial/ch02-rhf.md) · [open shell](../tutorial/ch03-open-shell.md) · [basis](../tutorial/ch04-basis-sets.md) · [SCF](../tutorial/ch18-scf-techniques.md).

Shares one SCF solver with [Kohn–Sham DFT](dft.md): the 2×2 product
{restricted, unrestricted} × {HF, KS}. Correlated post-HF lives in
[Post-Hartree–Fock](post-hf.md). Legend: [Roadmap index](index.md).

## 1. Molecular Hartree–Fock

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| RHF | ✓ | ✓ | `cpu::scf::RHF::run(Frame&)` / `run(const Frame&)` / `forces(Frame&)`. ADC, direct RPA, BSE, GF2, and EOM are not members. CUDA `cuda::scf::RHF` : `cuda::Driver<scf::Phase>` (fock, mix, diagonalize, density, check); device Fock/density/CDIIS; STO-3G s/p (l ≤ 1) native ERI+forces. PySCF energy: H₂/STO-3G on device; H₂O/6-31G on CPU (6-31G is not a device ERI kernel) |
| UHF | ✓ | ✓ | `cpu::scf::UHF{basis, n_alpha, n_beta, charge}` : `run(Frame&)` / `run(const Frame&)` / `forces(Frame&)`. CUDA `cuda::scf::UHF{basis, n_alpha, n_beta, charge}` : `cuda::Driver<scf::Phase>` with `run(Frame&)` / `forces(Frame&)` only (no const overload). `run` writes `kEnergy` [eV]; `forces` returns `GradientResult` and writes `kAtom_f*` [eV/Å]. Device Fock, density, CDIIS, and `launch_uhf_forces_sp` for STO-3G l ≤ 1. The device start is the spin-averaged MINAO density (`UHF.init_guess_by_minao`; H, C, N, O, F), so triplet O₂ matches PySCF's default rather than the core-Hamiltonian local minimum; any other element keeps that core start. Outside STO-3G l ≤ 1 the device call throws (`launch_eri_sp`) — no host-ERI fill. `CCSD<UHF>::eom_sf_ccsd`, `SfTDA<UHF>::run`, and `SfTDHF<UHF>::run` are not members |
| ROHF | ✓ | ○ | Guest–Saunders. `ROHF::run` / `forces`. Correlation is `MP2<ROHF>::run`, not a member |
| GHF | ✓ | ○ | Real two-component. `GHF::run` / `forces` |
| Complex / spinor HF | ○ | ○ | |
| ΔSCF / MOM / IMOM | ✓ | ○ | `UHF::delta_scf(Frame&, from, to, MomStyle)` |
| RHF analytic gradient | ✓ | ✓ | `RHF::forces(Frame&)` and `cuda::scf::RHF::forces`. Same angular-momentum envelope as the energy |
| UHF / ROHF / GHF analytic gradient | ✓ | UHF ✓ | `forces(Frame&)`. CUDA UHF uses `launch_uhf_forces_sp` on STO-3G l ≤ 1 (same envelope as the energy). ROHF and GHF gradients stay on the CPU |
| RHF / UHF analytic Hessian | ✓ | ○ | CPHF: `Hessian<RHF>::run(Frame&)` / `Hessian<UHF>::run(Frame&)` (`pyscf.hessian.rhf` / `uhf`) |
| ROHF Hessian | ✓ | ○ | `Hessian<ROHF>::run(Frame&)` — `central_difference_nuclear_hessian` of analytic `forces` (not CPHF; PySCF has no `hessian.ROHF`) |
| RHF static polarizability | ✓ | ○ | `Polarizability<RHF>::run()` after `RHF::run` (no Frame). Not an SCF method |
| Stability analysis | ✓ | ○ | `Stability<RHF>::run` after `RHF::run` (Seeger–Pople). Not an SCF method |

### Gaps

- CUDA Hartree–Fock is **RHF and UHF** on STO-3G l ≤ 1. No KS, ROHF, or GHF on device. `cuda::scf::UHF` throws outside that envelope (`launch_eri_sp`, `launch_uhf_forces_sp`) and does not fill ERIs on the host — use `cpu::scf::UHF`. `cuda::scf::RHF::forces` throws outside the same envelope; RHF energy for other bases still builds integrals on the host.
- ROHF: SOSCF / level shift still throw. GHF: `MP2` / `FCI` / soscf throw.
- Post-HF, continuum solvent, and CPHF response are not `MeanFieldMethod` members. SCF entries are `run(Frame&)` / `run(const Frame&)` / `forces(Frame&)`. Beside that: `ADC<RHF>::adc2` / `adc2x` / `ip_adc2` / `ea_adc2`, `DirectRPA<RHF>::run`, `BSE<RHF>::run`, `GF2<RHF>::run`, `G0W0<RHF>::run`, `Stability<RHF>::run`, `TDA<RHF>::run`, `TDHF<RHF>::run`, `CISD<RHF>::run`, `CCSD<RHF>::eom_ee_ccsd` / `eom_ip_ccsd` / `eom_ea_ccsd`, `CC2<RHF>::eom_cc2`, `CCSD<UHF>::eom_sf_ccsd`, `SfTDA<UHF>::run`, `SfTDHF<UHF>::run`, `MP2<Method>::run`, `WithSolvent<Method>::run`, `Hessian<RHF|UHF|ROHF>::run`, `Polarizability<RHF>::run`.

## 2. Periodic Hartree–Fock

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| Γ-point periodic RHF | ✓ | ○ | Truncated lattice; `RHF.periodic`; short-range erfc gate |
| k-point periodic RHF | ✓ | ✓ | `cpu::scf::KRHF` and `cuda::scf::KRHF`: `{basis, k_mesh}.run(Frame&)`. Orthorhombic cell, s functions, Monkhorst–Pack mesh, FFT Coulomb/exchange with Ewald exchange divergence. Device step is J/K, Fock, `zheevd`, density. `energy()` [Ha], `kEnergy` [eV]. l > 0 throws (`krhf_ao_lgt0`). Analytic forces throw (`krhf_nuclear_forces`). Tutorial: [k-point HF](../tutorial/methods/periodic-kpoint-hf.md) |
| k-point UHF periodic HF | ○ | ○ | |

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
