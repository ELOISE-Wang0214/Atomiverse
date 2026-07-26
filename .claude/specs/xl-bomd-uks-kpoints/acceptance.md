# Acceptance — xl-bomd-uks-kpoints

## Correctness

| ID | Check | Where |
| -- | ----- | ----- |
| A1 | Table-I + real/`XlBomdCplx` fixed-point algebra still green | `test_cpu_gpw_xl_bomd` |
| A2 | Multi-k RKS `guess_P_k` plumbing + real-on-k throw unchanged | unit |
| A3 | Multi-k UKS dual-spin density guess: both channels accepted; partial/mismatched count throws | unit |
| A4 | Converged multi-k UKS stores `P_k` / `P_k_beta` size `n_kpts`; warm dual guess reaches cold fixed-point energy | unit |
| A5 | `GPWEvaluator::build` **accepts** multi-k UKS + `xl_bomd`; still rejects OT+XL, bad K | unit |
| A6 | e2e multi-k UKS H₂ NVE (\(n_\alpha=n_\beta=1\)): \(\delta\le 0.05\), SCF every step with dual XL \(P_k^\sigma\), `n_kpts>=2`, PE bound, ⟨S²⟩≈0, finite coords, soft bond window | `reg_xl-bomd-uks-kpoints` |

## Perf

| Field | Criterion |
| ----- | --------- |
| `within_10x` | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) when a fair offline CP2K AIMD wall exists |
| `perf_exempt` | allowed if no CP2K binary on host (same policy as `xl-bomd-kpoints` / `xl-bomd-uks`) |

## Dual gate

- Correctness: A1–A6 green
- Speed: `within_10x` **or** `perf_exempt=true` with reason

## Done when

`acceptance_met=true` on worker report; unit + regression targets green;
commit on `wf/xl-bomd-uks-kpoints`. Spec deleted at merge by harness.
