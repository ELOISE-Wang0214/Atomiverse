# gpw-kpoint-stress

## Goal

Ship **analytic (Hellmann–Feynman fixed-\(P\)) cell stress for closed-shell
GPW RKS at Monkhorst–Pack / explicit k-points** — the natural follow-on after
Γ analytic stress (`gpw-analytic-stress` ✓) and k-point energy
(`gpw-kpoints-si` ✓).

Today `gpw_stress` / `energy_fixed_P` throw on multi-k
(`"k-point GPW not supported in this foothold"`). Multi-SCF strain-FD
(`gpw_stress_fd`) already works because it only re-runs SCF energy.

## Scope

**In**

- Analytic product stress \(\sigma\) [Ha/Bohr³] Voigt for periodic closed-shell
  RKS GPW with `k_mesh` / explicit `k_points` (same foothold constraints as
  the k energy path: single grid, `local_pp_on_grid`, periodic, no UKS/OT).
- k-summed fixed-\(P\) Lagrangian under orthorhombic diagonal strain at fixed
  fractional coordinates:
  \[
  E_L[\{P_k\},\{W_k\};\varepsilon]
  = \sum_k w_k\bigl[\mathrm{Re}\,\mathrm{Tr}(P_k T_k)
    + \mathrm{Re}\,\mathrm{Tr}(P_k V_{\mathrm{nl}}^k)
    - \mathrm{Re}\,\mathrm{Tr}(W_k S_k)\bigr]
  + E_{\mathrm{locgrid}}[\rho]
  + E_{\mathrm{es}}[\rho-n_c]
  + E_{\mathrm{xc}}[\rho]
  + E_{\mathrm{disp}}.
  \]
  \(\sigma_{aa} = -\Omega^{-1}\partial E_L/\partial\varepsilon_{aa}\) (CP2K
  `STRESS_TENSOR` sign). Off-diagonal shear remains 0 (orthorhombic only).
- Store occupied MO coefficients per k on `GpwResult` (`C_occ_k`) so stress
  rebuilds \(P_k\) / \(W_k\) after the reference SCF without a second solve.
- Unit: analytic vs multi-SCF FD on a cheap k-mesh H₂ toy.
- Regression under `regressions/gpw/gpw-kpoint-stress/`: Si (bct 4-atom)
  DZVP-MOLOPT-GTH + GTH-PBE-q4, offline CP2K pin or dual-gate + frozen band
  if no CP2K host.

**Out**

- UKS stress, triclinic shear, hybrid/meta-GGA, DFT+U/D4 stress at k.
- Pure closed-form Pulay AO stress tables without fixed-\(P\) energy FD.
- k-point forces (separate slug `gpw-kpoint-forces`).

## Public surface

- `gpw_stress` / `gpw_stress_fd` — already public; analytic gains multi-k.
- `GpwResult`: add `C_occ_k` (complex occupied MOs per k; also useful for
  forces).
- No public `*Config`, no `make_*` / `build_*`.

## Physics

At the SCF fixed point (closed shell, occupancy 2):

\[
W_k = \sum_i^{n_\mathrm{occ}} 2\,\varepsilon_i^k\, c_i^k (c_i^k)^\dagger,
\quad
P_k = \sum_i^{n_\mathrm{occ}} 2\, c_i^k (c_i^k)^\dagger.
\]

Density \(\rho(\mathbf r)=\sum_k w_k\sum_{\mu\nu}P^k_{\mu\nu}
\chi_\mu^k(\mathbf r)\,\overline{\chi_\nu^k}(\mathbf r)\) via
`collocate_density_k`. Strain FD holds \(\{P_k,W_k\}\) fixed while rebuilding
geometry-dependent \(S_k,T_k,V_{\mathrm{nl}}^k\), grid, Poisson, XC.

Refs: Nielsen & Martin, PRB **32**, 3780 (1985); VandeVondele et al., CPC
**167**, 103 (2005); ATV Γ path `internal/gpw_stress.cpp`; k machinery
`internal/kpoints.h`.

## Acceptance (summary)

See `acceptance.md`. Dual gate: unit analytic↔FD + regression pin/identities.
Perf: \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) or `perf_exempt` if no fair
CP2K wall.
