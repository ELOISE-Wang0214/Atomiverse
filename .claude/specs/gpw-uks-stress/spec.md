---
slug: gpw-uks-stress
status: acceptance_met
created: 2026-07-25
domain: gpw
queue: cp2k-method-queue
depends_on: [gpw-uks-forces, gpw-analytic-stress]
---

# Spec: `gpw-uks-stress` — Γ-point UKS analytic cell stress

**Slug:** `gpw-uks-stress`  
**Domain:** GPW (Quickstep) unrestricted Kohn–Sham  
**Queue:** CP2K capability loop  
**Depends on:** `gpw-uks-forces` ✓, `gpw-analytic-stress` ✓ (closed-shell fixed-P)

## Summary

Product **Hellmann–Feynman fixed-\(P\)** analytic cell stress for **Γ-point UKS**
GPW. Today `gpw_stress` / `energy_fixed_P` throw on `is_uks()` at Γ and multi-k.
Closed-shell Γ and k-point RKS stress and Γ-UKS forces are product; this closes
the UKS stress residual for orthorhombic diagonal Voigt components.

## Domain basis

Reference SCF at the unstrained cell yields spin densities \(P_\alpha,P_\beta\)
and energy-weighted densities

\[
W_\sigma = \sum_i^{n_\sigma} \varepsilon_i^\sigma\,
  c_i^\sigma (c_i^\sigma)^\top
  \qquad (\sigma=\alpha,\beta).
\]

The fixed-\(P\) Lagrangian under orthorhombic strain \(\varepsilon\) (fixed
fractional coordinates; rebuild geometry-dependent \(T\), \(V_{\mathrm{pp}}\),
\(S\), grid, Poisson, spin XC; **do not** re-converge SCF):

\[
\begin{aligned}
E_L[P_\alpha,P_\beta;W_\alpha,W_\beta;\varepsilon]
&= \mathrm{Tr}\!\big((P_\alpha+P_\beta)\,T\big)
 + \mathrm{Tr}\!\big((P_\alpha+P_\beta)\,V_{\mathrm{pp}}\big) \\
&\quad + E_{\mathrm{es}}[\rho_\alpha+\rho_\beta-n_c]
 + E_{\mathrm{xc}}[\rho_\alpha,\rho_\beta]
 + E_{\mathrm{disp}} \\
&\quad - \mathrm{Tr}(W_\alpha S) - \mathrm{Tr}(W_\beta S).
\end{aligned}
\]

Stress (CP2K `STRESS_TENSOR` sign):

\[
\sigma_{aa} = -\Omega^{-1}\partial E_L/\partial\varepsilon_{aa},
\quad a\in\{x,y,z\},
\quad \sigma_{\mathrm{shear}}=0
\quad\text{(orthorhombic v1)}.
\]

Multi-SCF strain-FD (`gpw_stress_fd`) is the dual-gate twin (already works for
UKS because it only re-runs SCF energy). Absolute offline CP2K pin is energy
from the same O₂ twin as `gpw-uks-forces` (and optional NUMERICAL stress when
CP2K is available).

## Scope

| In | Out |
| -- | --- |
| Γ-point UKS, PADE / PBE / PBEsol spin XC | k-point UKS stress (separate residual) |
| Product `gpw_stress` + evaluator `compute_stress` | Closed-form Pulay AO tables without fixed-P FD |
| Diagonal orthorhombic Voigt (xx,yy,zz) | Triclinic shear; DFT+U / D4 / hybrid stress |
| Dual-gate vs `gpw_stress_fd` + energy pin | Live CP2K in ctest |

## Public API

No new public type. Existing:

- `cpu::gpw::gpw_stress(...)` — stop throwing on Γ UKS
- `cpu::gpw::GPWEvaluator` with `n_alpha`/`n_beta`, `compute_stress=true`

Parameters remain members / ctor args (no `*Config`, no `make_*`).

## Acceptance (summary)

See `acceptance.md`. Green = dual-gate analytic vs multi-SCF FD on O₂ UKS
within tol + energy twin vs offline CP2K pin + unit gate + scoped build.

## Twin system

O₂ triplet, same geometry as `gpw-uks-forces` / `gpw-uks-pbe`:

- R(O–O) = 1.2075 Å, centered in 12 Å cubic cell, `PERIODIC NONE` / MT
- DZVP-MOLOPT-GTH + GTH-PBE-q6, PBE, cutoff 400 Ry, Nα=7, Nβ=5
