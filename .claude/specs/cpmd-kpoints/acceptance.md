# Acceptance — cpmd-kpoints

## Correctness

| ID | Check | Where |
| -- | ----- | ----- |
| A1 | Existing real `Cpmd` RATTLE / thermostat unit tests still green | `test_cpu_gpw_cpmd` |
| A2 | `CpmdCplx` seed + advance: \(C^\dagger S C = I\), Herm\((C^\dagger S\dot C)\)≈0; empty seed throws | unit |
| A3 | Multi-k RKS fixed-density (`max_iter==0` + `guess_P_k`) returns `iters==0`, finite energy, `H_ks_k`/`S_k` size == `n_kpts` | unit |
| A4 | Multi-k fixed-density without `guess_P_k` throws | unit |
| A5 | `GPWEvaluator::build` allows `cpmd` + non-Γ `k_mesh` (RKS); throws multi-k UKS+CPMD, OT+CPMD, XL+CPMD | unit |
| A6 | e2e H₂ multi-k CP NVE: relative \(E_{\mathrm{cons}}\) drift \(\delta\le 0.10\), production SCF iters == 0, `n_kpts >= 2`, PE bound, finite coords, soft bond window | `reg_cpmd-kpoints` |

## Perf

| Field | Criterion |
| ----- | --------- |
| `within_10x` | \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) when a fair offline CP2K CP AIMD wall exists |
| `perf_exempt` | allowed if no CP2K binary on host (same policy as `cpmd-water` / `xl-bomd-kpoints`) |

## Dual gate

- Correctness: A1–A6 green
- Speed: `within_10x` **or** `perf_exempt=true` with reason

## Done when

`acceptance_met=true` on worker report; unit + regression targets green;
commit on `wf/cpmd-kpoints`. Spec deleted at merge by harness.
