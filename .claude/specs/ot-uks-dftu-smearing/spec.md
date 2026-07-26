# Spec — ot-uks-dftu-smearing

## Goal

Productize **multi-k collinear UKS OT + Dudarev DFT+U + metallic occupation
smearing** (Fermi–Dirac / Methfessel–Paxton) on `cpu::gpw::GpwDft`.

Prior product paths (each landed separately):

| Combo | Slug | Status |
| ----- | ---- | ------ |
| multi-k RKS OT + U + smear | `ot-dftu-smearing` | ✓ |
| multi-k UKS OT + smear | `ot-uks-smearing` | ✓ |
| multi-k UKS OT + U (integer) | `ot-uks-dftu` | ✓ |
| multi-k UKS OT + U + smear | **this** | residual |

## Current tree (before)

- `GpwDft::build` rejects
  `is_uks() && use_ot && !plus_u.empty() && smearing != none`.
- Multi-k UKS OT runtime throws when `do_dftu && use_smear`.
- Integer multi-k UKS OT+U uses OUTER_SCF + fixed-H OT.
- UKS OT+smear (no U) uses free-energy PR+ Armijo; `build_uks_k_energy`
  already folds \(V_U\) / \(E_U\) when `plus_u` is set.

## Design

Mirror multi-k **RKS** OT+U+smear:

1. Remove build-time reject and runtime throw for UKS OT+U+smear.
2. Keep OUTER_SCF only for **integer** multi-k UKS OT+U
   (`do_dftu && !use_smear`).
3. Route metallic U+smear through the existing free-energy PR+ Armijo path
   (full AO manifold, shared \(\mu\) over \(\alpha+\beta\), \(f_\sigma\in[0,1]\));
   \(V_U\) via `build_uks_k_energy` (already product).
4. Soften g-norm / energy-plateau policy when `do_dftu` (same as RKS).

Fixed point: multi-k UKS DIIS + PLUS_U + FD/MP (CP2K 2026.1 cannot OT+k or
OT+SMEAR; offline pin is the DIIS twin).

## API

Public members only — no `*Config`, no factories:

```cpp
GpwDft m;
m.n_alpha = 1; m.n_beta = 1; // equal-spin seed; free under smear
m.use_ot = true;
m.plus_u.push_back(PlusU{molrs::Element::H, 0, 4.0});
m.smearing = Smearing::fermi_dirac;
m.degauss = 0.02; // Ry
m.extra_bands = 2;
m.k_points = {{0.25, 0.0, 0.0}};
m.k_point_weights = {1.0};
```

## Out of scope

- Fair OT wall twin vs CP2K (CP2K aborts OT+k / OT+SMEAR).
- Analytic multi-k UKS DFT+U forces/stress changes (existing paths).
- Molecular post-HF / MetalWalls / CUDA.

## Acceptance

See `acceptance.md`. Green dual gate + offline pin in
`regressions/gpw/ot-uks-dftu-smearing/REFERENCE.md`.
