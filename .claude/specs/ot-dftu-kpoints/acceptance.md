# Acceptance — ot-dftu-kpoints

## Correctness

| Gate | Check | Tol |
|------|-------|-----|
| A1 | multi-k DIIS+U at k=Γ ≡ multi-k OT+U at k=Γ | 1e-4 Ha |
| A1 | `e_dftu` non-trivial; OT vs DIIS e_dftu | 1e-3 Ha |
| A2 | multi-k OT+U ≡ multi-k DIIS+U at κ=(¼,0,0) | 1e-4 Ha |
| A3 | MP 2×2×2 OT+U converges; n_kpts=4 (TR) | — |
| A3 | \|E(U=4)−E(U=0)\| activation floor | ≥ 1e-4 Ha |
| A3 | absolute E vs frozen first-green pin | 1e-3 Ha |
| unit | `GpwOt.MultiKOtPlusUMatchesDiisEnergy` | 1e-4 Ha |

## Perf

- Prefer \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\).
- **`perf_exempt` (accepted):** CP2K 2026.1 hard-aborts
  `OT not possible with kpoint calculations`. No fair OT wall twin exists;
  decisive gates are OT≡DIIS multi-k identities (same pattern as
  `ot-kpoints`). Document in REFERENCE.md.

## Product

- `reg_ot-dftu-kpoints` green
- Public docs on `use_ot` + `plus_u` + k-mesh no longer claim multi-k
  OT+DFT+U throws
- No public `*Config` / factories
