# Excited States and Many-Body Methods

Tutorial: [Electronic excitations](../tutorial/electronic-excitations.md).

Response and many-body methods on HF/KS references. Legend:
[Roadmap index](index.md).

## 1. SCF and linear response

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| ΔSCF / MOM / IMOM | ✓ | ○ | `UHF::delta_scf` |
| CIS / TDA (RHF) | ✓ | ○ | `TDA<RHF>::run` after `RHF::run` |
| TDHF / RPA (RHF) | ✓ | ○ | `RHF::tdhf` |
| CIS(D) | ✓ | ○ | On RHF TDA |
| TDA-TDDFT (RKS) | ✓ | ○ | Full MO-product LSDA ALDA + hybrid B3LYP/PBE0 |
| TDA-TDDFT (UKS) | ✓ | ○ | Spin-conserving LSDA |
| Full LR-TDDFT (RKS) | ✓ | ○ | Casida A/B; pure KS + hybrid B3LYP |
| Full LR-TDDFT (UKS) | ✓ | ○ | Pure LSDA |
| Spin-flip TDA / TDDFT | ✓ | ○ | `UHF::sf_tda` / `sf_tddft` |
| CIS/TDA analytic forces | ✓ | ○ | Default analytic; FD fallback |
| TDA NAC vectors | ✓ | ○ | True NAC; no FSSH consumer yet |
| Real-time / CVS / REW TDDFT | ○ | ○ | |

## 2. ADC

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| ADC(2) | ✓ | ○ | EE on RHF |
| ADC(2)-x | ✓ | ○ | |
| IP-ADC(2) / EA-ADC(2) | ✓ | ○ | |
| ADC(3) / CVS-ADC / ADC gradient | ○ | ○ | |

## 3. EOM coupled cluster

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| EOM-EE-CCSD | ✓ | ○ | `CCSD<RHF>::eom_ee_ccsd` (singlet). Not an `RHF` member |
| EOM-IP / EA-CCSD | ✓ | ○ | |
| EOM-SF-CCSD | ✓ | ○ | On UHF+UCCSD |
| EOM-CC2 | ✓ | ○ | |
| EOM-DIP / DEA / CVS / STEOM / gradients | ○ | ○ | |

## 4. Simplified / tight-binding response

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| TD-DFTB | ✓ | ○ | `DFTB2::td_casida`; energies + length oscillator strengths; Γ, closed shell; no absolute oracle |
| TD-DFTB excited gradient / NAC | ○ | ○ | |
| sTDA / sTDDFT / sTDA-xTB | ○ | ○ | |

## 5. Green’s function and many-body

| Capability | CPU | CUDA | Notes |
| ---------- | :-: | :--: | ----- |
| G0W0 (spectral) | ✓ | ○ | On RHF; also periodic on GPW |
| Direct RPA | ✓ | ○ | `RHF::direct_rpa`; periodic on GPW |
| BSE@G0W0 | ✓ | ○ | |
| One-shot GF2 | ✓ | ○ | |
| Self-consistent GW / RPAX / SOSEX | ○ | ○ | |

Periodic MP2 / RPA / GW depth is tabulated under
[DFT § Periodic](dft.md#4-periodic-kohnsham).
