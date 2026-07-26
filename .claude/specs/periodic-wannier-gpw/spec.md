---
slug: periodic-wannier-gpw
status: acceptance_met
created: 2026-07-26
domain: gpw / response
queue: cp2k-method-loop
depends_on: [periodic-wannier-si]
unblocks: [periodic-wannier-gpw]
---

# Spec: `periodic-wannier-gpw` — multi-k MLWF on GPW

**Slug:** `periodic-wannier-gpw`  
**Domain:** GPW response (Maximally Localized Wannier Functions)  
**Queue:** `cp2k-method-loop` wave 12  
**Twin:** CP2K Quickstep `&LOCALIZE` / Wannier centers on GTH+MOLOPT  
**Prior art:** PW Marzari–Vanderbilt (`cpu::pw::PeriodicWannier`,
`reg_periodic-wannier-si` / `uks`); molecular Boys (`cpu::scf::WannierCenters`).

## Summary

Deliver **periodic multi-k Marzari–Vanderbilt MLWF** on the CPU GPW method
(`GpwDft` / Gaussian AO + plane-wave density). Residual after PW MV product:
GPW-grid / AO-basis localization with the same physical pin — diamond Si
valence WFs sit on tetrahedral bond midpoints (Marzari & Vanderbilt, Phys.
Rev. B 56, 12847 (1997) §IV).

Public API: concrete `cpu::gpw::PeriodicWannier` with members (`k_mesh`,
`max_iter`, `tol`, `step`, `use_projection`, `proj_alpha_ang`,
`trial_centers`); `compute(GpwDft&, atoms)` runs multi-k RKS SCF (insulator)
then localizes. No public `*Config`, no factories. Centers [Å], spreads [Å²]
on `PeriodicWannierResult`.

## Domain basis

Occupied manifold only (fixed integer occupations; reject metallic smearing
and UKS in this foothold). Neighbor overlaps on an \(N^3\) Monkhorst–Pack
mesh (`k_time_reversal = false` for a complete neighbor graph):

\[
M_{mn}^{(\mathbf k,\mathbf b)}
  = \sum_{\mu\nu}
    C^*_{\mu m}(\mathbf k)\,
    S^{\mathbf k,\mathbf k+\mathbf b}_{\mu\nu}\,
    C_{\nu n}(\mathbf k+\mathbf b),
\]

with the phase-weighted AO cross-overlap

\[
S^{\mathbf k,\mathbf q}_{\mu\nu}
  = \sum_{\mathbf T}
    e^{i\mathbf q\cdot\mathbf T}
    \langle\varphi_\mu|
      e^{-i(\mathbf q-\mathbf k)\cdot\mathbf r}|
      \varphi_\nu(\cdot-\mathbf T)\rangle
\]

evaluated analytically on contracted Cartesian Gaussians (product-center
phase × Gaussian form factor; exact for s, product-phase for higher \(\ell\)).

Projection initialization onto bond-centered real-space Gaussians (BCT
diamond auto midpoints when 4 Si + 8 WFs) and optional MV steepest-descent
polish. Centers from the Berry-phase formula
\(\mathbf r_n = -N_k^{-1}\sum_{\mathbf k,\mathbf b} w_b\,\mathbf b\,
\operatorname{Im}\ln M_{nn}\).

## System pin (regression)

Diamond Si in the **body-centred-tetragonal 4-atom** orthorhombic cell
(GPW is orthorhombic-only):

| | |
| --- | --- |
| Lattice | \(a=5.43\) Å; cell \((a/\sqrt2,\,a/\sqrt2,\,a)\) |
| Atoms | 4 Si; DZVP-MOLOPT-GTH + GTH-PBE-q4 |
| Valence | 16 e⁻ → `n_occ = 8` → 8 MLWFs |
| XC | PADE (LDA) |
| Density cutoff | 100 Ry (foothold) |
| MLWF mesh | \(2\times2\times2\), TR off |
| Pin | literature bond-midpoint geometry (not a snapshot of ATV output) |

## Acceptance

See `acceptance.md`. Hard gates in `reg_periodic-wannier-gpw`:

1. SCF converges; `n_wann == 8`.
2. After greedy assignment + min-image, every center lies within **0.90 Å**
   of a unique Si–Si bond midpoint (same coarse SMV window as PW Si reg).
3. Spread reduction: \(\Omega_\mathrm{final} < 0.5\,\Omega_\mathrm{init}\).
4. Unit tests cover pure math (midpoints, SMV weight, phase→center) and
   argument validation — no full SCF in `tests/`.

## Out of scope (follow-on)

- Collinear UKS / spin MLWF on GPW
- Metallic disentanglement
- Non-orthorhombic cells
- Analytic Wannier forces / spreads derivatives
- Dense-mesh sub-0.15 Å absolute center placement

## Perf

`perf_exempt` when no offline CP2K wall twin is available on the pin host
(same policy as `periodic-wannier-si` / `gpw-kpoints-si`).
