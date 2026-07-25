# Spec — `gpw-kpoints-uks`

## Goal

Collinear unrestricted Kohn–Sham GPW on a Monkhorst–Pack k mesh (CP2K
`&KPOINTS` + `LSD` / UKS), building on:

- `gpw-kpoints-si` — multi-k RKS Bloch machinery
- `gpw-uks-pbe` / `gpw-uks-o2-gamma` — Γ-UKS dual H / spin XC

## Public API

Concrete `cpu::gpw::GpwDft` members already present:

- `n_alpha` / `n_beta` — spin occupations (UKS)
- `k_mesh[3]` / `k_points` + `k_point_weights` / `k_time_reversal`
- `XcKind::PADE | PBE | PBEsol` for UKS (including multi-k)

No public `*Config`, no factories. Frame energy still eV via evaluators;
internal Ha/Ry OK when cited.

## Physics

1. Per-k complex Hermitian S_k / T_k / V_nl^k (existing k-point helpers).
2. Dual spin densities \(P^\alpha_k\), \(P^\beta_k\) with occupation factors 1
   (not closed-shell 2).
3. k-summed \(\rho_\alpha,\rho_\beta\) → Poisson on \(\rho_\alpha+\rho_\beta-n_\mathrm{ion}\)
   + spin XC (`xc_grid_*_spin`) → dual \(V_\mathrm{eff}^{\alpha/\beta}\).
4. Dual \(H_k^\sigma\), dual Pulay DIIS, fixed fillings \(N_\alpha,N_\beta\) at every k.
5. Time-reversal fold exact for collinear UKS without SOC (same as RKS density).

## Constraints (no-fallback)

Multi-k still rejects: OT, DFT+U, multigrid, nonperiodic MT, analytic local PP
off, level shift. BLYP not on the spin-GGA path.

## Acceptance (summary)

See `acceptance.md`. Dual gate: correctness identities + absolute pin;
perf ≤ 10× CP2K or `perf_exempt`.

## Out of scope

Analytic k-point forces (`gpw-kpoint-forces`), k-point stress, Γ-UKS forces,
metallic smearing, noncollinear spin.
