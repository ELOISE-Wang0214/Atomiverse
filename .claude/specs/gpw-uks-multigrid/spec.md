---
slug: gpw-uks-multigrid
status: approved
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [gpw-uks-o2-gamma]
oracle: cp2k
---

# Spec: `gpw-uks-multigrid` — Γ-UKS multigrid product collocation

**Slug:** `gpw-uks-multigrid`  
**Domain:** GPW / Quickstep-style multigrid density collocation  
**Queue:** `cp2k-method-loop` (wave 9)  
**Oracle:** offline CP2K (NGRIDS / REL_CUTOFF UKS twin)

## Summary

Wire **CP2K-style multigrid product collocation** on the **Γ-point UKS** SCF
path. RKS already collocates primitive-pair products per level when
`ngrids > 1 && rel_cutoff > 0` (`build_products` → prolong FFT → Hartree/XC on
the fine grid → restrict → `integrate_products`). UKS currently **throws**
`"multigrid product collocation is not wired"` and forces `ngrids=1`. This
slug removes that gate and routes dual spin densities \(P_\alpha, P_\beta\)
through the same product multigrid assembly used for closed-shell \(P\).

## Gap (current tree)

| Surface | Status |
|---------|--------|
| RKS `ngrids>1` product collocation | ✓ product |
| UKS single-grid (`CompactBasisGrid`) | ✓ product (`gpw-uks-o2-gamma`, …) |
| UKS `ngrids>1` | throws fail-fast |
| k-point + multigrid | still rejected (`build`: k-path does not support multigrid) — **out of scope** |

## Physics / algorithm

For each spin channel \(\sigma \in \{\alpha,\beta\}\):

1. Collocate products of \(P_\sigma\) on every multigrid level.
2. Prolong each coarse \(\rho_{\sigma,\ell}\) to the finest grid (spectral FFT;
   adjoint of restrict).
3. Form \(\rho = \rho_\alpha + \rho_\beta\) on the fine grid; Poisson + spin XC
   (`xc_grid_*_spin` / PADE LSD) as today.
4. Restrict \(v_{\mathrm{eff},\sigma}\) to each level; integrate products into
   \(H_\sigma\).

`ngrids==1` remains the bit-identical single-grid CompactBasisGrid path.

## API surface

- `GpwDft::{ngrids, rel_cutoff}` already public (solver).
- Expose matching **members** on `GPWEvaluator` (`ngrids`, `rel_cutoff`) and
  forward in `build()` — parameters as members, no `*Config` bag, no factories.
- Frame energy still **eV** on `frame.meta.get(kEnergy)`; internal pin in Ha.

## Acceptance (see `acceptance.md`)

1. Unit: Γ-UKS with `ngrids=4`, `rel_cutoff>0` converges (toy H₂ high-spin);
   closed-shell UKS multigrid energy matches RKS multigrid to tight tol.
2. Regression `regressions/gpw/gpw-uks-multigrid`: O₂ triplet UKS PADE +
   DZVP-MOLOPT-GTH + MT, **CUTOFF 400 / NGRIDS 4 / REL_CUTOFF 60**, offline
   CP2K pin (same twin family as `gpw-uks-o2-gamma`).
3. Dual gate: correctness pin + \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\)
   or documented `perf_exempt`.

## Out of scope

- Multi-k + multigrid
- Multigrid analytic forces/stress special path (forces may still use solver
  settings; no new force pin required for this energy-foothold slug)
- OT-only multigrid (OT UKS already goes through `build_uks_H` once wired)

## Domain basis

Lippert–Hutter–Parrinello GPW; CP2K Quickstep `&MGRID` (CUTOFF, REL_CUTOFF,
NGRIDS). Adjoint prolong/restrict keeps the KS matrix the gradient of the
collocated energy (existing `multigrid.h` contract).
