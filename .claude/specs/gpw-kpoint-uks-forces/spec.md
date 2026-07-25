# Spec — `gpw-kpoint-uks-forces`

## Goal

Ship **analytic ionic forces for collinear multi-k UKS GPW** — the residual
named after Γ-UKS forces (`gpw-uks-forces` ✓), multi-k RKS forces
(`gpw-kpoint-forces` ✓), and multi-k UKS energy (`gpw-kpoints-uks` ✓).

Today `gpw_forces_from_scf` / `gpw_forces_analytic` throw
`"k-point UKS forces are not implemented"` when `uses_kpath(cfg) && cfg.is_uks()`.
k-UKS SCF also discards occupied Bloch MOs (stores only α `eps_k`), so forces
cannot rebuild \(P_k^\sigma / W_k^\sigma\).

## Scope

**In**

- Analytic \(F_A = -\partial E/\partial R_A\) for periodic collinear UKS GPW
  with `k_mesh` / explicit `k_points` (same foothold constraints as k-energy:
  single grid, `local_pp_on_grid`, PADE/PBE/PBEsol spin XC).
- k-summed Pulay / core / nonlocal with \(P_\mathrm{tot}^k = P_\alpha^k + P_\beta^k\)
  and \(W^k = W_\alpha^k + W_\beta^k\) (occupancy 1 per spin orbital).
- Grid Hellmann–Feynman with spin-resolved \(V_\mathrm{eff}^{\alpha/\beta}\)
  and Bloch \(\chi^k\) derivatives (same structure as Γ UKS + k-RKS).
- Store α/β occupied MO coefficients and β band energies per k on `GpwResult`
  so `gpw_forces_from_scf` rebuilds densities without a second SCF.
- Unit: multi-k UKS forces finite, Newton-3 residual, optional analytic↔FD.
- Regression under `regressions/gpw/gpw-kpoint-uks-forces/` with offline pin
  (FD identity and/or frozen first-green; CP2K ENERGY_FORCE when host allows).

**Out**

- UKS analytic stress (`gpw-uks-stress` residual).
- Metallic smearing / variable occupations at k.
- Hybrid / meta-GGA UKS, OT+UKS, multigrid product UKS.
- Molecular Gaussian UKS forces (PySCF loop).
- CUDA ports.

## Public surface

- `gpw_forces_analytic(atoms, cfg)` / `gpw_forces_from_scf` — already public;
  remove k-UKS throw; k-path UKS works when `n_alpha`/`n_beta` + k mesh set.
- `GpwResult`: add `C_occ_k_beta` + `eps_k_beta` (α remains `C_occ_k` / `eps_k`).
- No public `*Config`, no `make_*` / `build_*`.

## Physics

At the multi-k UKS SCF fixed point (integer fillings \(N_\alpha, N_\beta\) at every k):

\[
\begin{aligned}
F_A &= \sum_k w_k\Bigl[
  \mathrm{Tr}\big((W_\alpha^k+W_\beta^k)\,\partial S_k/\partial R_A\big)
  - \mathrm{Tr}\big((P_\alpha^k+P_\beta^k)\,\partial H_k^\mathrm{core}/\partial R_A\big)
\Bigr] \\
&\quad - \sum_k w_k\sum_{\sigma\in\{\alpha,\beta\}}
  2\int V_\mathrm{eff}^{\sigma}(r)\,
  \mathrm{Re}\big[(\partial\chi_\mu^k/\partial R_A)\,(P_\sigma^k\chi^k)_\mu\big]\,dr \\
&\quad + \text{ionic HF} + \text{local-PP HF (grid route)},
\end{aligned}
\]

with \(P_\sigma^k = \sum_{i\in\mathrm{occ}_\sigma} c_{\sigma i}^k (c_{\sigma i}^k)^\dagger\)
and \(W_\sigma^k = \sum_i \varepsilon_{\sigma i}^k\, c_{\sigma i}^k (c_{\sigma i}^k)^\dagger\).
Complex bra-side center derivatives use \(2\,\mathrm{Re}(\cdots)\) as in
k-RKS forces. Electrostatics / local PP spin-independent; XC from the same
spin-GGA/LSD kernels as k-UKS SCF.

Refs: VandeVondele et al., CPC **167**, 103 (2005); ATV Γ UKS forces;
k-RKS forces `assemble_gpw_forces_kpoints`; k-UKS energy
`run_gpw_scf_kpoints_uks`.

## Acceptance (summary)

See `acceptance.md`. Dual gate: unit + regression identities/pins.
Perf: \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) or `perf_exempt` if no fair wall.
