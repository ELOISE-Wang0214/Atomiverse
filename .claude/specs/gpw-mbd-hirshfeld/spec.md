# Spec — `gpw-mbd-hirshfeld`

## Goal

Ship **Hirshfeld volume ratios** \(V/V_{\mathrm{free}}\) from the GPW SCF
valence density and use them to **density-scale** free-atom TS / MBD@rsSCS
parameters (Tkatchenko–Scheffler 2009; Ambrosetti et al. 2014). Closes the
roadmap §8 residual after free-atom TS ✓ (`gpw-ts-ice`) and free-atom MBD ✓
(`gpw-mbd-ice`, volume_ratio = 1).

## Method

Hirshfeld weight and effective volume (atomic units):

\[
w_A(\mathbf r)=\frac{\rho_A^{\mathrm{free}}(|\mathbf r-\mathbf R_A|)}
{\sum_B\rho_B^{\mathrm{free}}(|\mathbf r-\mathbf R_B|)},\quad
V_A=\int w_A(\mathbf r)\,\rho(\mathbf r)\,|\mathbf r-\mathbf R_A|^3\,d^3r
\]

\[
v_A=\frac{V_A}{V_A^{\mathrm{free}}},\quad
V_A^{\mathrm{free}}=\int\rho_A^{\mathrm{free}}(r)\,r^3\,d^3r
\]

Free-atom valence densities are in-tree **exponential** models
\(\rho_Z^{\mathrm{free}}(r)=A_Z\exp(-\alpha_Z r)\) with \(A_Z\) fixed by
GTH valence count \(N_Z\) and static \(\alpha_Z\) (no file parsing). Scaling of
TS/MBD free-atom parameters matches the existing internal API:

\[
\alpha\leftarrow\alpha_{\mathrm{free}}\,v,\quad
C_6\leftarrow C_{6,\mathrm{free}}\,v^2,\quad
R\leftarrow R_{\mathrm{free}}\,v^{1/3}.
\]

## Public API

- `cpu::gpw::GpwDft::hirshfeld = true` (and `GPWEvaluator::hirshfeld`)
- Requires `ts` **or** `mbd` (mutually exclusive flags unchanged)
- Γ-point RKS product path only in v1 (k-points / UKS throw)
- Dispersion energy is **post-SCF** when `hirshfeld` (density-dependent);
  free-atom path (`hirshfeld=false`) remains geometry-only
- Ratios on `GpwResult::volume_ratios` (empty when free-atom)
- Energy on `GpwResult::e_disp` [Hartree], folded into `e_total`
- No public `*Config`; no `make_*` / `build_*` factories

## Out of scope

- Analytic volume-force contributions (\(d v_A / d R\))
- Periodic Ewald MBD (still isolated MBD@rsSCS geometry)
- UKS / multi-k Hirshfeld
- All-electron free-atom densities / iterative Hirshfeld-I
- vdW-DF / XDM / VV10
- Live CP2K or live libMBD in ctest

## Offline oracle

Hard pin is **libMBD** MBD@rsSCS energy with the ATV Hirshfeld ratios (dual
pure-Python Hirshfeld on the same free-atom exponential model + same water-dimer
geometry). CP2K stock has no TS/MBD `PAIR_POTENTIAL`; notes under
`tools/cp2k_oracle/inputs/gpw-mbd-hirshfeld/`.

## Acceptance

See `acceptance.md`.
