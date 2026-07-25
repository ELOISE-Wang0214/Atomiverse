# Spec — `gpw-rpbe-water`

## Goal

Ship **RPBE** (Hammer–Hansen–Nørskov, *Phys. Rev. B* **59**, 7413 (1999)) as a
pure GGA on the GPW grid path, with a textbook isolated-water dual-gate
regression against PBE. Distinct from **revPBE** (Zhang & Yang κ-revision):
RPBE keeps PBE κ/μ but uses the exponential enhancement

\[
F_x^{\mathrm{RPBE}}(s) = 1 + \kappa\bigl(1 - e^{-\mu s^2/\kappa}\bigr),
\]

with PBE correlation (`GGA_C_PBE`). Public surface: `XCType::RPBE` /
`XcKind::RPBE`; no `*Config`, no factories. Frame energy remains eV; pins
may cite Ha.

## Scope

**In**

- Pointwise `evaluate_xc(XCType::RPBE, …)` (scalar + batch).
- `XcKind::RPBE` on closed-shell Γ GPW RKS (`GpwDft` / `GPWEvaluator`).
- Grid helpers `xc_grid_rpbe` via shared `xc_grid_gga`.
- Unit tests (metadata, σ=0 ≡ PBE, finite-σ dual-gate) + e2e regression
  `regressions/gpw/gpw-rpbe-water/`.
- Offline CP2K twin input `tools/cp2k_oracle/inputs/gpw-rpbe-water.inp`
  (`XC_FUNCTIONAL RPBE`).

**Out**

- Spin-resolved (UKS) RPBE; k-point RPBE; forces/stress special-casing
  beyond the existing RKS GGA dispatch.
- Live CP2K in ctest; molecular pure-PySCF paths; CUDA.

## Acceptance

1. `XCType::RPBE` / `XcKind::RPBE` compile and are pure GGA (`is_gga`,
   `hybrid_exchange_fraction == 0`).
2. At σ = 0, \(E_{\mathrm{xc}}^{\mathrm{RPBE}} = E_{\mathrm{xc}}^{\mathrm{PBE}}\)
   (shared PW92 LDA base).
3. At finite σ, \(|E_{\mathrm{xc}}^{\mathrm{RPBE}} - E_{\mathrm{xc}}^{\mathrm{PBE}}| > 0\).
4. Regression: isolated H₂O MT / DZVP-MOLOPT-GTH / 400 Ry converges;
   \(|E_{\mathrm{RPBE}} - E_{\mathrm{PBE}}| > 10^{-4}\) Ha; soft absolute
   band vs first-green ATV (or offline CP2K when available).
5. Public Frame energy finite [eV]; no external process spawn.
6. `perf_exempt` if no CP2K wall twin; else \(r = T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\).

## Science notes

- CP2K name: `XC_FUNCTIONAL RPBE` (libxc `GGA_X_RPBE` + `GGA_C_PBE`).
- Not the same as `revPBE` / `GGA_X_PBE_R` (κ = 1.245, rational form).
- GTH-PBE PP + MOLOPT basis reused (standard CP2K water twin family).
