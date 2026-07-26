# Acceptance — gpw-d4-forces

| # | Gate | Pass criterion |
| - | ---- | -------------- |
| 1 | Offline E_disp | \|E_disp − (−9.9458842745846774e−4 Ha)\| ≤ 1e−8 Ha (periodic water dimer, dftd4 3.7.0) |
| 2 | Offline ∇E pin | max\|F_atv + g_dftd4\| ≤ 1e−7 Ha/Bohr (ATV force = −dftd4 gradient; ATV atom order OHH OHH) |
| 3 | Self FD | max\|F_product − F_energy_fd\| ≤ 1e−7 Ha/Bohr on standalone D4 |
| 4 | Wiring | max\|(F_gpw_d4 − F_gpw_plain) − F_d4\| ≤ 1e−8 Ha/Bohr (additive) |
| 5 | Dual-gate product | max\|F_an − F_fd_multiSCF\| ≤ 5e−5 Ha/Bohr with `d4=true` (same functional) |
| 6 | Public surface | no new `*Config` / `make_*`; Frame forces still eV/Å on write-back if exposed |
| 7 | Perf | within_10x **or** perf_exempt (no CP2K DFTD4 build) |

`acceptance_met=true` only when 1–6 green and 7 within_10x or exempt.
