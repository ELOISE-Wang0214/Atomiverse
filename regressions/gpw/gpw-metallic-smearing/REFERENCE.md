# Reference — metallic Fermi–Dirac smearing (BCC Na, GPW multi-k)

## Method

Spin-unpolarized **GPW** Kohn–Sham DFT with **Fermi–Dirac** electronic
smearing (Mermin free energy \(F = E_{\mathrm{KS}} - TS\)). Public API:
`cpu::gpw::GpwDft` members `smearing`, `degauss` [Ry], `extra_bands`,
`mp_order`. Occupations reuse `cpu::pw::internal` (FD / Methfessel–Paxton).

## System (Atomiverse pin)

| Quantity | Value |
| -------- | ----- |
| Geometry | BCC Na conventional cubic, \(a = 4.225\) Å, 2 atoms |
| Basis | `DZVP-MOLOPT-SR-GTH` (Na) — in-tree table |
| Pseudopotential | `GTH-PBE-q9` (9 valence e⁻ / atom → \(N_e = 18\)) |
| XC | PADE (LDA) |
| Density cutoff | 150 Ry |
| k-mesh | Monkhorst–Pack \(2\times2\times2\), time-reversal fold |
| Smearing | Fermi–Dirac, `degauss = 0.02` Ry |
| SCF | DIIS 6, `tol = 1e-6` Ha, `EPS_DEFAULT = 1e-7` |

### Offline free-energy pin (frozen ATV)

| Observable | Pin | Tol |
| ---------- | --- | --- |
| Free energy \(F\) | **−92.126330608378 Ha** | \(10^{-4}\) Ha |
| \(N_{\mathrm{occ}}=\sum w f\) | **18.0** | \(10^{-5}\) |
| \(E_F\) | **0.055826 Ha** (finite) | \(\lvert E_F\rvert < 10^3\) Ha |
| Entropy term \(-\sigma S\) | **−2.26×10⁻⁴ Ha** | finite |

Recorded 2026-07-26 on Arrhenius (x86_64) first green of
`reg_gpw-metallic-smearing` (7 SCF iters, wall ≈ 1.9 s single rank).

## Offline CP2K twin (algorithm / SMEAR documentation)

**No CP2K binary on the worker host** at pin time. Absolute CP2K digits are
not invented. Documented product counterpart:

| Item | Value |
| ---- | ----- |
| Input (Al twin) | `tools/cp2k_oracle/inputs/al_fermi_dirac.inp` |
| Method | Quickstep GPW, `&SMEAR METHOD FERMI_DIRAC` |
| Notes | Al primitive FCC is non-orthorhombic (ATV GPW is ortho-only); Na BCC uses in-tree tables. Same SMEAR + multi-k algorithm class. |

Na CP2K twin (to run offline when CP2K is available): same BCC cell, PADE,
`GTH-PBE-q9` + `DZVP-MOLOPT-SR-GTH`, MP \(2\times2\times2\),
`ELECTRONIC_TEMPERATURE [K] 3157.75` (≈ 0.02 Ry \(k_B T\)).

## Perf gate

| Engine | Wall |
| ------ | ---- |
| Atomiverse `reg_gpw-metallic-smearing` | **1.90 s** (single rank, first green) |
| CP2K twin | **n/a** (no `cp2k` / module on worker host) |

`perf_exempt = true` — no fair offline CP2K wall on the same Na setup.
`within_10x` not scored.
