# Acceptance — xl-bomd-uks

## Correctness

| ID | Check | Where |
| -- | ----- | ----- |
| A1 | Table-I coefficients + fixed-point algebra still green (RKS unit path) | `test_cpu_gpw_xl_bomd` |
| A2 | `GpwDft::run` accepts UKS density guess when **both** `guess_P` (α) and `guess_P_beta` (β) are non-null; shape mismatch throws | unit |
| A3 | UKS with only one of `{guess_P, guess_P_beta}` throws (no silent cold start) | unit |
| A4 | RKS with non-null `guess_P_beta` throws | unit |
| A5 | OT + any density guess still throws | unit |
| A6 | `GPWEvaluator::build` allows `xl_bomd` + `n_alpha/n_beta`; still throws for `use_ot` / out-of-table K / `cpmd` / `second_gen_cp` | unit |
| A7 | e2e O₂ NVE: relative energy drift \(\delta\le 0.05\), SCF converged every sampled step with XL spin guesses, PE bound, finite coords, soft O–O ceiling | `reg_xl-bomd-uks` |

## Perf

| Field | Criterion |
| ----- | --------- |
| `within_10x` | `T_atv / T_oracle ≤ 10` when a fair offline CP2K AIMD wall exists |
| `perf_exempt` | allowed if no CP2K binary on host (same policy as `xl-bomd-water`) |

## Dual gate

- Correctness: A1–A7 green
- Speed: `within_10x` **or** `perf_exempt=true` with reason

## Done when

`acceptance_met=true` on worker report; regression + unit targets green;
commit on `wf/xl-bomd-uks`. Spec deleted at merge by harness.
