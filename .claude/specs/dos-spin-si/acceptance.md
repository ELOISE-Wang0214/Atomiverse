# Acceptance — dos-spin-si

## Correctness gates (regression `reg_dos-spin-si`)

| # | Gate | Criterion |
| - | ---- | --------- |
| 1 | SCF | `scf_converged` |
| 2 | Spin | `nspin == 2`, `|magnetization| ≤ 1e-6` (M=0 fixed occ) |
| 3 | Spectrum | DOS grid non-empty; energy/dos length match |
| 4 | \(n_e\) | \(\lvert n_{\mathrm{elec\_to\_vbm}} - 8\rvert \le 1.0\) |
| 5 | Gap | `gap_ry > 0` (gapped Si twin) |
| 6 | α/β | per-spin DOS present when product exposes them; total integral consistent |

## Unit gates (`test_cpu_pw_bands_dos`)

- `gaussian_dos` with `spin_degeneracy = 1` integrates to weight (not 2×).
- `DosResult` default remains restricted (`nspin == 1`).

## Perf

- `within_10x` if CP2K wall measured; else `perf_exempt = true` with reason.

## Done when

- Dual gate green on `wf/dos-spin-si`.
- Spec deleted on merge to `dev` (harness).
