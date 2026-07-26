# Spec — `gpw-vdw-scf-potential`

## Goal

Productize **self-consistent nonlocal vdW potential** for GPW:

| Functional | Flag | CP2K twin |
| ---------- | ---- | --------- |
| vdW-DF (DRSLL / DF1) | `GpwDft::vdw_df` | `NON_LOCAL TYPE DRSLL` |
| vdW-DF2 (LMKLL / DF2) | `GpwDft::vdw_df2` | `NON_LOCAL TYPE LMKLL` |
| rVV10 | `GpwDft::rvv10` | `NON_LOCAL TYPE RVV10` |

Beyond post-SCF \(E_{\mathrm{disp}}\) (`gpw-vdw-df-ice` / `gpw-vdw-df2-ice` /
`gpw-rvv10-ice`): each SCF iteration adds \(V_{\mathrm{nl}}[\rho] =
\delta E_c^{\mathrm{nl}}/\delta\rho\) to the KS effective potential; analytic
forces include \(V_{\mathrm{nl}}\) through the existing \(\int V_{\mathrm{eff}}
\partial\chi/\partial R\) assembly.

## Algorithm (Román–Pérez–Soler + CP2K `qs_dispersion_nonloc.F`)

1. \(q_0(\rho,|\nabla\rho|^2)\) with soft saturation + \(\partial q_0/\partial\rho\),
   \(\partial q_0/\partial|\nabla\rho|\).
2. \(\theta_\alpha = \rho\,p_\alpha(q_0)\) (DRSLL/LMKLL) or
   \(c\,\rho^{3/4}p_\alpha(q_0)\) (rVV10); FFT → \(k\).
3. \(u_\alpha(k)=\sum_\beta\phi_{\alpha\beta}(k)\,\theta_\beta(k)\);
   \(E_c^{\mathrm{nl}}=\frac12\sum w\,\mathrm{Re}(\theta^*u)\) (rVV10: extra \(\frac12\) + \(\beta\int\rho\)).
4. IFFT \(u_\alpha\to r\); Soler eq. 10 potential + divergence of gradient flux
   (same spectral `divergence` as GGA XC).
5. SCF: `veff += v_nl` every `build_H_and_energy`; final \(E_c^{\mathrm{nl}}\) on
   SC density → `e_disp`. Forces: same `v_nl` on force `veff`.

## Scope

| In | Out |
| -- | --- |
| Γ RKS single-grid (product path) | Multigrid (already throws) |
| DRSLL / LMKLL / rVV10 potential + energy | Stock VV10 (separate slug) |
| Analytic forces via \(V_{\mathrm{nl}}\) in `veff` | Analytic cell stress (virial \(k\)-terms) |
| UKS/multi-k: total-\(\rho\) \(V_{\mathrm{nl}}\) on both spins / shared grid | Full spin-polarized kernel |

## API

No new public flags. Existing `vdw_df` / `vdw_df2` / `rvv10` become **SC**
(potential in SCF + forces). Post-SCF energy-only ice regressions update pins
if SC density shifts \(E_{\mathrm{disp}}\) beyond tol.

Internal: `vdw_df_drsll` / `vdw_df_lmkll` / `rvv10_energy` gain optional
`v_r` / full `*_evaluate` returning energy + potential.

## Acceptance

See `acceptance.md`.

## Deps

`gpw-vdw-df-ice` ✓ · `gpw-vdw-df2-ice` ✓ · `gpw-rvv10-ice` ✓
