# Spec — `ot-uks-dftu`

**Slug:** `ot-uks-dftu`  
**Method:** gpw  
**Branch:** `wf/ot-uks-dftu`  
**Oracle:** CP2K offline (DIAG multi-k UKS + PLUS_U twin; no OT+k on CP2K 2026.1)

## Gap

Multi-k collinear **UKS OT** is product (`reg_ot-uks-kpoints`) and multi-k
**RKS OT + Dudarev DFT+U** is product (`reg_ot-dftu-kpoints`). Combining
`use_ot` + open-shell `n_alpha`/`n_beta` + `plus_u` on the complex multi-k
path throws:

```
GPW k-UKS OT: DFT+U is not wired on the multi-k OT path
```

## Product target

Wire multi-k UKS OT + Dudarev PLUS_U via the same **OUTER_SCF + fixed-H OT**
pattern as multi-k RKS OT+U:

1. Rebuild \(H_k^\sigma\) (including \(V_U^\sigma\)) from the k-summed spin
   densities.
2. Per \((k,\sigma)\): `ot_minimize_fixed_k` on fixed \(H_k^\sigma\).
3. Outer loop on total energy \(E_{\mathrm{KS}}+E_U\) until \(\Delta E <\) tol.

Fixed point must match multi-k UKS DIIS + PLUS_U (shared energy assembly).

## Out of scope

- UKS OT + metallic smearing (`ot-uks-smearing`)
- OT + DFT+U + smear (`ot-dftu-smearing`)
- Forces / stress with multi-k OT+U (energy SCF only)
- Multigrid / nonperiodic multi-k

## API

Public members only on `cpu::gpw::GpwDft`: `use_ot`, `n_alpha`/`n_beta`,
`k_mesh` / `k_points`, `plus_u`. No `*Config`, no factories.

## Regression

`regressions/gpw/ot-uks-dftu/`:

| Gate | Check | Tol |
| ---- | ----- | ---: |
| A1 | multi-k DIIS+U at Γ ≡ multi-k OT+U at Γ (UKS) | 1e-4 Ha |
| A2 | multi-k OT+U ≡ multi-k DIIS+U (κ=¼) | 1e-4 Ha |
| A3 | MP mesh OT+U converges; \(\|E(U)-E(0)\|\) activation; e_dftu | floors |
| A4 | absolute vs first-green / DIIS twin pin | 1e-3 Ha |

## Perf

CP2K 2026.1 aborts OT+KPOINTS → **perf_exempt** (same as
`ot-dftu-kpoints` / `ot-uks-kpoints`). Decisive gates are OT≡DIIS identities.

## Acceptance

- Spec + offline pin / REFERENCE.md
- Multi-k UKS OT + `plus_u` no longer throws
- `reg_ot-uks-dftu` green ≤ 60 s Release
- Commit on `wf/ot-uks-dftu`
