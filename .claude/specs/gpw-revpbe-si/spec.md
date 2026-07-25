# Spec — `gpw-revpbe-si`

**Slug:** `gpw-revpbe-si`
**Loop:** CP2K method (GPW / PW / periodic)
**Oracle:** CP2K offline pin (`tools/cp2k_oracle`)
**Branch:** `wf/gpw-revpbe-si`
**Parent:** `dev`

## Goal

Ship **revPBE** pure GGA as a first-class GPW `XcKind` on a textbook Si₂
fragment, dual-gated against PBE and soft-pinned to an offline CP2K twin
(or first-green ATV band when CP2K is unavailable on the pin host).

## Science

- **revPBE** (Zhang & Yang, *Phys. Rev. Lett.* **80**, 890 (1998)): same
  PBE functional form with revised exchange enhancement parameter
  \(\kappa = 1.245\) (PBE uses \(\kappa = 0.804\)); \(\mu\) and correlation
  \(\beta\) match PBE. libxc: `GGA_X_PBE_R` + `GGA_C_PBE`. CP2K:
  `XC_FUNCTIONAL revPBE`.
- Shared pointwise engine: extend `cpu::scf::XCType` + `evaluate_xc`; GPW
  grid path reuses `xc_grid_gga` with the new type.

## Scope (in)

1. `XCType::RevPBE` in `cpu/scf/xc_functional.h` (+ `is_gga`, hybrid
   fractions 0).
2. `pbe_exchange_spin` accepts \(\kappa\); revPBE uses \(1.245\).
3. `XcKind::RevPBE` on GPW public surface; RKS closed-shell path (forces /
   stress dispatch included for consistency).
4. Thin wrappers `xc_grid_revpbe` (+ unit test constant-density LDA limit
   and finite-gradient difference vs PBE).
5. Regression `regressions/gpw/gpw-revpbe-si/` — Si₂, DZVP-MOLOPT-GTH,
   GTH-PBE-q4, MT isolated (same budget pattern as `gpw-pbesol-si`),
   dual-gate \(|E_{\mathrm{revPBE}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha.
6. Offline twin input `tools/cp2k_oracle/inputs/gpw-revpbe-si.inp`.

## Scope (out)

- RPBE (Hammer form) — separate slug `gpw-rpbe-water`.
- UKS revPBE spin path (UKS currently PADE/PBE/PBEsol only).
- Periodic MOLOPT Si e2e (budget; twin deck may document periodic for
  offline re-pin).
- Roadmap table flip / queue done row (harness on merge).

## Public API

```cpp
cpu::gpw::GpwDft m;
m.xc = cpu::gpw::XcKind::RevPBE;
// members: cutoff_ry, n_occ, cell, nonperiodic, … — no *Config, no make_*
m.build(atoms);
auto res = m.run(atoms); // res.e_total [Ha]
```

## Acceptance

See `acceptance.md`.
