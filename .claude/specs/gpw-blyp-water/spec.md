# Spec — `gpw-blyp-water`

## Goal

Ship **GPW RKS BLYP** (Becke 1988 exchange + Lee–Yang–Parr correlation) on the
closed-shell grid GGA path, dual-gated against offline CP2K when available and
against in-tree PBE on the same geometry. Textbook case: isolated water,
DZVP-MOLOPT-GTH + GTH-PBE, Martyna–Tuckerman (same family as
`gpw-water-pbe-mt`).

## Scope

| In | Out |
| -- | --- |
| `XcKind::BLYP` on GPW RKS Γ | UKS BLYP / spin GGA |
| `xc_grid_gga` via `XCType::BLYP` (already in molecular `evaluate_xc`) | Hybrid B3LYP on GPW |
| Energy dual-gate + regression `reg_gpw-blyp-water` | Forces/stress BLYP-specific e2e (share grid path) |
| Offline CP2K twin input under `tools/cp2k_oracle/inputs/` | Live CP2K in ctest |

## Design

1. **Extend** `enum class XcKind { PADE, PBE, PBEsol, BLYP }`.
2. **Wire** all RKS XC grid dispatches (SCF, forces, stress) to call
   `xc_grid_gga(..., XCType::BLYP)` when `xc == BLYP` (thin `xc_grid_blyp`
   wrapper, same shape as `xc_grid_pbe` / `xc_grid_pbesol`).
3. **Reuse** molecular pure-GGA engine — no second B88/LYP implementation.
4. **UKS** still rejects non-`PADE` (unchanged; follow-on `gpw-uks-pbe` owns
   spin GGA).
5. **Public API:** `GPWEvaluator(..., XcKind::BLYP, …)` / `GpwDft::xc =
   XcKind::BLYP` — ctor args / members only; no `*Config`, no factories.
6. **Units:** internal Ha; Frame energy still eV via evaluator.

## Acceptance (summary)

See `acceptance.md`. Dual gate: SCF converges; finite \(E_{\mathrm{BLYP}}\);
\(|E_{\mathrm{BLYP}}-E_{\mathrm{PBE}}| > 10^{-4}\) Ha on the water MT twin;
absolute band vs offline CP2K pin when available, else first-green ATV soft
band until re-pin. Perf: \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) or
`perf_exempt` when no CP2K binary on host.

## Non-goals

- revPBE / RPBE / OLYP (separate queue slugs)
- GPW hybrids / meta-GGA
- Molecular all-electron BLYP (already on PySCF path)
