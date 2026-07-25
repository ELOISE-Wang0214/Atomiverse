# Reference — bulk Fe collinear UKS + metallic smearing (plane-wave)

## Method

Collinear spin-polarized plane-wave KS (`cpu::pw::PW` with
`spin_polarized = true`) **together with** metallic Fermi–Dirac smearing
(`smearing = fermi_dirac`, `degauss` [Ry], `extra_bands`). Shared chemical
potential \(\mu\) over spin-up and spin-down spin-orbitals; occupations
\(f_{\sigma kb}\in[0,1]\); magnetization \(M = N_\uparrow - N_\downarrow\) is
**self-consistent** (multiplicity only seeds the initial \(\rho_\uparrow/
\rho_\downarrow\) split). Total energy is the Mermin free energy
\(F = E_{\mathrm{KS}} - TS\).

Closes the product gap left by:

| Foothold | What it lacked |
| -------- | -------------- |
| `pw-smearing-al` | spin |
| `pw-spin-polarized-o2` / `pw-uks-kmesh-fe` | metallic smearing |

`PW::run` previously threw if both were requested.

## System

| Quantity | Value |
| -------- | ----- |
| Geometry | BCC Fe primitive cell, \(a_{\mathrm{cubic}} = 2.87\) Å, 1 atom at origin |
| PP | SG15 ONCV PBE Fe `tests/data/pw/Fe_ONCV_PBE-1.2.upf` (\(Z_{\mathrm{val}}=16\)) |
| XC | PBE |
| \(E_{\mathrm{cut}}\) | 12 Ry (algorithm pin, not production cutoff) |
| k-mesh | Monkhorst–Pack \(2\times2\times2\) |
| Smearing | Fermi–Dirac, `degauss = 0.02` Ry |
| Extra bands | 4 (per spin; Davidson states \(=\max(\lceil N_e/2\rceil+\mathrm{extra},\,N_\uparrow^{\mathrm{seed}})\)) |
| Multiplicity seed | 5 (initial \(M\approx 4\); free under smear) |
| SCF | `tol=2\times10^{-3}`, `max_iter=30`, `mix_beta=0.3` |

## Pins

### Electron count

\[
N_{\mathrm{occ}} = \sum_{\sigma kb} w_k f_{\sigma kb} = 16
\quad (| \Delta N | \le 10^{-4}).
\]

### Magnetization (ferro seed)

Self-consistent \(M > 0.1\). Recorded value on the settings above:
**\(M \approx 0.991\)**. Low \(E_{\mathrm{cut}}\) undershoots experimental
BCC Fe (~2.2 \(\mu_B\)); the gate is the ferro sign + stability under smear.

### Free energy (offline ATV freeze)

\[
F = -205.8365363098~\mathrm{Ry}
\quad (| \Delta F | \le 0.05~\mathrm{Ry}).
\]

Recorded 2026-07-25 on Arrhenius (x86_64), OMP=8, single rank, 11 SCF
iterations, wall ≈ 20 s.

### Fermi level

Finite \(E_F\) with \(\lvert E_F\rvert < 10^3\) Ry (recorded ≈ 2.16 Ry).

## Offline CP2K twin (algorithm / documentation)

**Not an absolute-energy twin** (GTH-GPW ≠ NC-UPF PW). Documents that CP2K
UKS + `&SMEAR METHOD FERMI_DIRAC` is the product counterpart.

| Item | Value |
| ---- | ----- |
| Note | `tools/cp2k_oracle/inputs/pw-uks-fe-smearing.note` |
| Host CP2K | not loadable as fair NC-UPF PW twin at pin time |
| Method class | Quickstep UKS + FERMI_DIRAC + Fe GTH (if re-run offline) |

No live CP2K at ctest.

## Perf

| Field | Value |
| ----- | ----- |
| perf_exempt | true |
| reason | no fair single-rank NC-UPF PW twin wall vs CP2K GTH-GPW; same class as `pw-uks-kmesh-fe` |
| T_atv (OMP=8) | **~20 s** |
| ratio T_atv/T_cp2k | n/a |

## Dual gate

1. Correctness: converge + \(N_{\mathrm{occ}}\) + free-energy pin + \(M>0.1\) + finite \(E_F\).
2. Perf: exempt (no fair CP2K NC-UPF twin).
