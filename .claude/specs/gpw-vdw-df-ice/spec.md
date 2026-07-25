# Spec — gpw-vdw-df-ice

## Intent

Ship **nonlocal vdW-DF correlation (DRSLL / DF1)** on GPW as a post-SCF
additive, matching CP2K `VDW_POTENTIAL NON_LOCAL TYPE DRSLL` / libvdwxc DF1
energy on the same density. Closes roadmap §8 vdW-DF ○ foothold (after TS ✓ /
free-atom MBD ✓).

## Public API

- `GpwDft::vdw_df` / `GPWEvaluator::vdw_df` (bool members; **no** `*Config` bag)
- Mutually exclusive with `d3` / `d3bj` / `d4` / `ts` / `mbd`
- Energy on `GpwResult::e_disp` [Hartree], folded into `e_total`
- Periodic single-grid only (throws if nonperiodic / multigrid)

## Algorithm

1. \(q_0(\rho,\sigma)\) with DF1 \(Z_{ab}=-0.8491\), soft \(q_{\mathrm{cut}}=5\)
2. \(\theta_\alpha(\mathbf r)=\rho\,p_\alpha(q_0)\) cubic splines
3. FFT → reciprocal convolution with community RPS kernel \(\phi_{\alpha\beta}(k)\)
4. \(E_c^{\mathrm{nl}}=\frac12\sum_k w_k \mathrm{Re}[\theta^*\phi\theta]\)

Kernel tables: float32 zlib of libvdwxc 0.5.0 / QE community tables
(Soler 2009; Dion 2004). Inflate on first use (no runtime file I/O).

## Out of scope (follow-ons)

- Self-consistent \(V_c^{\mathrm{nl}}\) in the KS potential
- vdW-DF2 / optB88-vdW / rVV10
- Multigrid / UKS / k-point density path
- Analytic ionic forces / stress of \(E_c^{\mathrm{nl}}\)

## Acceptance

See `acceptance.md`. Dual gate: unit test vs libvdwxc synthetic pin +
regression `reg_gpw-vdw-df-ice` vs offline E_disp pin; perf_exempt.
