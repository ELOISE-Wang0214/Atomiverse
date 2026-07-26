# Acceptance — `xl-bomd-linearized-forces`

## Correctness

| ID | Gate | Where |
| -- | ---- | ----- |
| A1 | `xl_bomd_linearized` without `xl_bomd` throws at `build` | unit |
| A2 | `xl_bomd_linearized` + UKS / multi-k / OT / CPMD / SGCP throws at `build` | unit |
| A3 | Shadow energy recovers KS energy when \(D = n\): \(\mathcal{E}[n,n] = E_{\mathrm{KS}}[n]\) (bitwise / 1e-12) | unit |
| A4 | After cold seed, warm force steps report `iters == 0` (zero iterative SCF) | unit (H₂ tiny) or reg |
| A5 | E2E H₂O NVE: relative total-energy drift ≤ 0.10 over 1 warm + 4 production steps; finite PE/KE; O–H < 1.8 Å; warm steps `iters == 0` | `reg_xl-bomd-linearized-forces` |

## Perf

| Field | Value |
| ----- | ----- |
| `perf_exempt` | **true** when no CP2K binary for fair AIMD wall |
| Else | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) |

## Dual gate

Green only when A1–A5 hold and (`within_10x` or `perf_exempt`).
