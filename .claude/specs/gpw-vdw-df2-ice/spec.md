# Spec — gpw-vdw-df2-ice

## Intent

Ship **nonlocal vdW-DF2 correlation (LMKLL / DF2)** on GPW as a post-SCF
additive, matching CP2K `VDW_POTENTIAL NON_LOCAL TYPE LMKLL` / libvdwxc DF2
energy on the same density. Closes roadmap §8 vdW-DF2 ○ after DRSLL ✓
(`reg_gpw-vdw-df-ice`).

## Public API

- `GpwDft::vdw_df2` / `GPWEvaluator::vdw_df2` (bool members; **no** `*Config` bag)
- Mutually exclusive with `d3` / `d3bj` / `d4` / `ts` / `mbd` / `vdw_df`
- Energy on `GpwResult::e_disp` [Hartree], folded into `e_total`
- Periodic single-grid only (throws if nonperiodic / multigrid)

## Algorithm

Same Román–Pérez–Soler path as DRSLL; only \(Z_{ab}\) changes:

1. \(q_0(\rho,\sigma)\) with DF2 \(Z_{ab}=-1.887\), soft \(q_{\mathrm{cut}}=5\)
2. \(\theta_\alpha(\mathbf r)=\rho\,p_\alpha(q_0)\) cubic splines
3. FFT → reciprocal convolution with **shared** community RPS kernel
   \(\phi_{\alpha\beta}(k)\) (libvdwxc default; identical to DF1 tables)
4. \(E_c^{\mathrm{nl}}=\frac12\sum_k w_k \mathrm{Re}[\theta^*\phi\theta]\)

Lee et al., *Phys. Rev. B* **82**, 081101 (2010); Dion et al. PRL **92**,
246401 (2004); Soler RPS PRL **103**, 096102 (2009).

## Out of scope (follow-ons)

- Self-consistent \(V_c^{\mathrm{nl}}\) in the KS potential
- optB88-vdW / rVV10 / XDM
- Multigrid / UKS / k-point density path
- Analytic ionic forces / stress of \(E_c^{\mathrm{nl}}\)

## Acceptance

See `acceptance.md`. Dual gate: unit test vs libvdwxc DF2 synthetic pin +
regression `reg_gpw-vdw-df2-ice` vs offline E_disp pin; perf_exempt.
