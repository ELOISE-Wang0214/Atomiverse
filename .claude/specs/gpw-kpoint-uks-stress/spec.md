# Spec — `gpw-kpoint-uks-stress`

## Goal

Productize **analytic (Hellmann–Feynman fixed-\(P_k^\sigma\)) cell stress** for
**collinear multi-k UKS GPW**, closing the residual after:

- Γ-UKS stress ✓ (`reg_gpw-uks-stress`)
- multi-k RKS stress ✓ (`reg_gpw-kpoint-stress`)
- multi-k UKS forces ✓ (`reg_gpw-kpoint-uks-forces`)
- multi-k UKS SCF ✓ (`reg_gpw-kpoints-uks`)

Today `gpw_stress` throws `k-point UKS stress not implemented`.

## Science

Orthorhombic diagonal Voigt stress from the non-self-consistent Lagrangian

\[
\begin{aligned}
E_L[\{P_k^\sigma\},\{W_k^\sigma\};\varepsilon]
&= \sum_k w_k\Big[
  \mathrm{Re}\,\mathrm{Tr}(P_k^{\mathrm{tot}} T_k)
  + \mathrm{Re}\,\mathrm{Tr}(P_k^{\mathrm{tot}} V_{\mathrm{nl}}^k)
  - \mathrm{Re}\,\mathrm{Tr}(W_k^{\mathrm{tot}} S_k)
\Big] \\
&\quad + E_{\mathrm{locgrid}}[\rho]
  + E_{\mathrm{es}}[\rho-n_c]
  + E_{\mathrm{xc}}[\rho_\alpha,\rho_\beta]
  + E_{\mathrm{disp}}
\end{aligned}
\]

with \(P_k^{\mathrm{tot}}=P_k^\alpha+P_k^\beta\), \(W_k^{\mathrm{tot}}=W_k^\alpha+W_k^\beta\),
occupancy 1 per spin channel, densities from `collocate_density_k` on
\(P_k^\sigma\), and spin XC (PADE LSD / PBE / PBEsol). Strain FD at fixed
fractional coordinates; \(\sigma_{aa}=-\Omega^{-1}\partial E_L/\partial\varepsilon_{aa}\).

Oracle twin: multi-SCF strain-FD (`gpw_stress_fd`) on the same UKS+k method.

## System (regression)

| Item | Value |
|------|--------|
| Molecule | High-spin H₂ (M=3 → Nα=2, Nβ=0) |
| Positions [Bohr] | H (5.5,6,6), H (6.5,6,6) compressed bond |
| Cell | cubic 12 Bohr, PERIODIC XYZ |
| Basis / PP | DZV-GTH-PADE / GTH-PADE-q1 (in-tree) |
| XC | PADE LSD |
| Cutoff | 100 Ry, local_pp_on_grid |
| k mesh | MONKHORST-PACK 1×1×2 (TR-folded) |
| Stress step | \(h=2\times10^{-3}\) |

## API (no Config / no factories)

```cpp
GpwDft m;
m.n_alpha = 2; m.n_beta = 0;
m.k_mesh = {1,1,2};
m.xc = XcKind::PADE;
// ...
auto s = gpw_stress(atoms, m, /*h=*/2e-3); // Voigt Ha/Bohr³
```

Public energy still eV on Frame via `GPWEvaluator`; internal stress path remains
Ha/Bohr³ (same as siblings).

## Out of scope

- Shear / off-diagonal stress (orthorhombic only, = 0)
- DFT+U + k-UKS stress (`gpw-dftu-stress` separate)
- CUDA GPW
- Absolute same-basis CP2K NUMERICAL digits if module absent (dual-gate decisive)

## Deps

`gpw-kpoints-uks` ✓, `gpw-kpoint-stress` ✓, `gpw-uks-stress` ✓,
`gpw-kpoint-uks-forces` ✓
