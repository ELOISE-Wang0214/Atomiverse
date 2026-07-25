# Spec — `gpw-uks-forces`

## Goal

Ship **Γ-point UKS analytic ionic forces** for GPW (Quickstep), matching offline
CP2K `ENERGY_FORCE` on open-shell O₂ (triplet, PBE / DZVP-MOLOPT-GTH / MT).

## Context

- Γ-UKS **energy** is productized (`gpw-uks-o2-gamma` PADE, `gpw-uks-pbe` PBE).
- Closed-shell RKS analytic forces are productized (`gpw-forces-water`).
- Force assembly in `src/cpu/gpw/internal/gpw_forces.cpp` is **closed-shell only**:
  energy-weighted density uses \(W = \sum_i 2\varepsilon_i c_i c_i^\top\) and a
  single spin-unpolarized \(V_\mathrm{eff}\). UKS needs \(W = W_\alpha + W_\beta\)
  (occupancy 1 per spin) and spin-resolved grid potentials
  \(V_\mathrm{eff}^\alpha, V_\mathrm{eff}^\beta\).

## Public surface

- `cpu::gpw::GPWEvaluator` with `n_alpha`/`n_beta` + `compute_forces=true` writes
  forces to Frame (`kAtom_fx/y/z`) in **eV/Å** (internal Ha/Bohr).
- No new public `*Config` bags; no `make_*` factories.
- Internal: `gpw_forces_analytic` accepts UKS `GpwDft` (`is_uks()`).

## Science

At the UKS SCF fixed point:

\[
F_A = -\mathrm{Tr}\big((P_\alpha+P_\beta)\,\partial(T+V_\mathrm{pp})/\partial R_A\big)
      + \mathrm{Tr}\big((W_\alpha+W_\beta)\,\partial S/\partial R_A\big)
      - \sum_{\sigma\in\{\alpha,\beta\}} 2\sum_{\mu\in A}\int
         (\partial\chi_\mu/\partial R_A)\, V_\mathrm{eff}^\sigma\,
         (P_\sigma\chi)_\mu\,dr
      + \text{ionic HF} + \text{local-PP HF (grid route)}.
\]

with \(W_\sigma = \sum_{i\in\mathrm{occ}_\sigma} \varepsilon_{\sigma i}
c_{\sigma i} c_{\sigma i}^\top\). Electrostatics / local PP spin-independent;
XC from `xc_grid_pbe_spin` / `xc_grid_pade_lsd` matching the SCF.

## Out of scope

- k-point UKS forces (`gpw-kpoint-forces` / `gpw-kpoints-uks`)
- UKS analytic stress (`gpw_stress` still throws UKS)
- ROKS, OT+UKS, multigrid product UKS
- Molecular Gaussian UKS forces (PySCF loop)

## Acceptance summary

See `acceptance.md`. Dual gate: CP2K force pin + optional fair wall ≤ 10×.
