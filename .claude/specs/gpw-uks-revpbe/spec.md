# Spec — gpw-uks-revpbe

## Goal

Γ-point **UKS** GPW single-point energy with spin-resolved **revPBE** GGA
(`XcKind::RevPBE`), offline-pinned to CP2K 2026.1 Quickstep LSD on O₂ triplet
(same geometry as `gpw-uks-pbe` / `gpw-uks-blyp`).

## Scope

**In**

- `xc_grid_revpbe_spin` thin wrapper over `xc_grid_gga_spin(..., XCType::RevPBE)`.
- Γ-UKS SCF path accepts / dispatches `XcKind::RevPBE` (dual H, spin XC, ⟨S²⟩).
- Matching UKS force / Γ-UKS stress / k-UKS energy+force gates so RevPBE is not
  a silent “energy-only” trap (mirrors BLYP surface).
- Public `GPWEvaluator` / `GpwDft` members only — no `*Config`, no factories.
- Regression `regressions/gpw/gpw-uks-revpbe/` vs offline pin + PBE dual-gate.
- Unit tests: spin collapse, build accept, toy H₂ high-spin converges.

**Out**

- UKS **RPBE** (separate residual if needed).
- Multi-k UKS stress for RevPBE (k-UKS stress still PADE/PBE/PBEsol only — same
  as BLYP today).
- Molecular post-HF, classical MW, CUDA ports.

## Science

- Zhang–Yang revPBE: PBE form with \(\kappa = 1.245\) (exchange enhancement).
- CP2K 2026.1 does **not** accept short `XC_FUNCTIONAL revPBE`; pin uses nested
  `&PBE` / `PARAMETRIZATION revPBE`.
- Public Frame energy remains eV; internal pin cites Hartree.

## Acceptance (see acceptance.md)

1. Unit tests green (scoped `test_cpu_gpw_uks`).
2. Regression green vs CP2K pin |ΔE| ≤ 1e-3 Ha and |ΔS²| ≤ 0.15.
3. Dual-gate |E_revPBE − E_PBE| > 1e-4 Ha.
4. Perf: T_atv / T_cp2k ≤ 10 (or documented exempt).
