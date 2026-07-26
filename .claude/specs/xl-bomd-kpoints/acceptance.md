# Acceptance — xl-bomd-kpoints

## Correctness

| ID | Check | Where |
| -- | ----- | ----- |
| A1 | Table-I coefficients + real fixed-point algebra still green | `test_cpu_gpw_xl_bomd` |
| A2 | `XlBomdCplx` fixed point on constant complex Hermitian density; shape mismatch throws | unit |
| A3 | Multi-k RKS `GpwDft::run` accepts `guess_P_k` of size `n_kpts`; shape/count mismatch throws | unit |
| A4 | Multi-k path + real-only `guess_P` throws (no silent cold start) | unit |
| A5 | `GPWEvaluator::build` allows `xl_bomd` + non-Γ `k_mesh` (RKS); throws for multi-k UKS + XL, OT+XL, out-of-table K | unit |
| A6 | e2e H₂ multi-k NVE: relative energy drift \(\delta\le 0.05\), SCF converged every sampled step with XL \(P_k\) guesses, PE bound, finite coords, soft bond window, `n_kpts >= 2` | `reg_xl-bomd-kpoints` |

## Perf

| Field | Criterion |
| ----- | --------- |
| `within_10x` | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) when a fair offline CP2K AIMD wall exists |
| `perf_exempt` | allowed if no CP2K binary on host (same policy as `xl-bomd-water` / `xl-bomd-uks`) |

## Dual gate

- Correctness: A1–A6 green
- Speed: `within_10x` **or** `perf_exempt=true` with reason

## Done when

`acceptance_met=true` on worker report; unit + regression targets green;
commit on `wf/xl-bomd-kpoints`. Spec deleted at merge by harness.
