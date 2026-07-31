# Acceptance — `xl-bomd-linearized-uks`

## Correctness

| ID | Gate | Where |
| -- | ---- | ----- |
| A1 | `xl_bomd_linearized` + Γ UKS (`n_alpha`/`n_beta`) builds (no throw) | unit |
| A2 | `xl_bomd_linearized` + UKS + multi-k still throws at `build` (as does OT) | unit |
| A3 | Dual-spin shadow energy recovers KS energy at the fixed point \(D^\sigma = n^\sigma\); residual = \(\sum_\sigma \mathrm{Tr}(H^\sigma (D^\sigma - n^\sigma))\) | unit |
| A4 | E2E O₂-triplet NVE (Γ UKS, PADE, MT — identical twin of `xl-bomd-uks` + linearized tag): cold seed `iters > 0` and converged; warm steps `iters == 0` and converged; relative total-energy drift ≤ 0.10 over 1 warm + 4 production steps; ⟨S²⟩ in the triplet window [1.5, 2.5] after seed and production; finite PE/KE; 0.9 < r_OO < 1.6 Å | `reg_xl-bomd-linearized-uks` |

## Perf

| Field | Value |
| ----- | ----- |
| `perf_exempt` | **true** when no CP2K binary for fair AIMD wall (2026-07-31 inventory: none) |
| Else | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) |

## Dual gate

Green only when A1–A4 hold and (`within_10x` or `perf_exempt`).
