---
slug: ugw-uks-rpbe
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: gpw / periodic GW
queue: cp2k-method-loop
deferred: false
depends_on: [periodic-ugw-si, gpw-uks-rpbe]
unblocks: [ugw-uks-rpbe]
grilled: false
---

# Spec — `ugw-uks-rpbe`

## Intent

Close the **UGW@UKS RPBE** residual: spectral G0W0 on a converged
**Γ-point UKS GPW RPBE** reference must evaluate spin-resolved
\(v_{\mathrm{xc}}^\sigma\) via the same `xc_grid_rpbe_spin` product path
already used by UKS SCF / forces / stress. Today `GpwDft::g0w0` throws in
`eval_xc_grid_spin` for `XcKind::RPBE` while PBE / PBEsol / RevPBE / BLYP
UKS UGW are product.

## Gap (binding)

| Piece | Status on parent |
| ----- | ---------------- |
| UKS RPBE SCF (`XcKind::RPBE`, `xc_grid_rpbe_spin`) | ✓ `reg_gpw-uks-rpbe` |
| Γ UGW@UKS spectral G0W0 | ✓ `reg_periodic-ugw-si` (PBE) |
| Multi-k UGW | ✓ `reg_periodic-ugw-kpoints` |
| `g0w0` → `eval_xc_grid_spin(RPBE)` | **throws** (this slug) |

## Science

Spectral G0W0 quasiparticle equation (Hybertsen–Louie form; same as
`GpwDft::g0w0` product):

\[
\omega = \varepsilon_p + \mathrm{Re}\,\Sigma_{c,p}(\omega)
         + \Sigma_{x,p} - v_{\mathrm{xc},p}.
\]

For UKS: spin-conserving dRPA over α→α ⊕ β→β, same-spin \(G_\sigma W\),
MT \(\Sigma_x^\sigma\), spin \(v_{\mathrm{xc}}^\sigma\) from **RPBE**
(Hammer, Hansen, Nørskov 1999: exponential exchange enhancement + PBE
correlation; libxc `GGA_X_RPBE` + `GGA_C_PBE`).

When \(n_\alpha = n_\beta\) and \(C^\alpha \approx C^\beta\), **UGW ≡ RGW**
— dual-gate against closed-shell `GpwDft::g0w0` on the same RPBE deck.

**Not** revPBE (Zhang–Yang κ-revision of the rational PBE form).

## Public API (no Config / no factories)

```cpp
cpu::gpw::GpwDft m;
m.n_alpha = 4; m.n_beta = 4;   // UKS
m.xc = cpu::gpw::XcKind::RPBE;
// ... cutoff, cell, build(atoms); run(atoms);
const cpu::gpw::G0W0Result& gw = m.g0w0(atoms);
// gw.gap_ev [eV] α HOMO→LUMO QP gap
```

Parameters remain members of `GpwDft`. Energy on Frame (if any outer
Evaluate path) stays eV; internal QP orbitals in Ha as today.

## Implementation sketch

1. In `src/cpu/gpw/g0w0.cpp` `eval_xc_grid_spin`, replace the RPBE throw
   with `return xc_grid_rpbe_spin(grid, rho_a, rho_b);` (already used by
   SCF). Closed-shell `eval_xc_grid` already maps RPBE → `xc_grid_rpbe`.
2. No new public types. No change to `G0W0Result` schema.
3. Multi-k UGW path reuses the same `eval_xc_grid_spin` — one fix covers
   Γ and multi-k.

## Regression

`regressions/gpw/ugw-uks-rpbe/` — artificial sc-Si₂ twin of
`periodic-ugw-si` with `XcKind::RPBE`. Offline CP2K input twin under
`tools/cp2k_oracle/inputs/` (RI-G0W0 diagnostic only; not a hard
spectral digit twin — same contract as `periodic-ugw-si`).

## Out of scope

- New XC functionals beyond RPBE spin in G0W0
- Frequency-integrated RI-GW matching CP2K absolute QP digits
- Periodic HF exchange / hybrid GW
- CUDA
- Live CP2K in ctest
- Molecular post-HF G0W0 (PySCF loop)

## Perf

`perf_exempt` — CP2K RI-G0W0 ≠ spectral plasmon UGW + MT exchange (same
as `periodic-ugw-si` / `periodic-ugw-kpoints`). Wall budget ≤ 60 s.
