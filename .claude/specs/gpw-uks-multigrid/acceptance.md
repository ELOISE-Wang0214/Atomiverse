# Acceptance — `gpw-uks-multigrid`

## Correctness

| Check | Criterion |
|-------|-----------|
| No throw | `GpwDft` UKS + `ngrids>1` + `rel_cutoff>0` runs (no “not wired”) |
| Unit H₂ UKS multigrid | converges; finite \(E\); high-spin \(\langle S^2\rangle \approx 2\) |
| Unit closed-shell | UKS \(n_\alpha=n_\beta=1\) multigrid \(E\) ≡ RKS multigrid \(E\) (≤ 1e-6 Ha) |
| Regression O₂ | \(\lvert E_{\mathrm{ATV}}-E_{\mathrm{CP2K}}\rvert \le 1\times10^{-3}\) Ha (same pin class as `gpw-uks-o2-gamma`; multigrid ATV path) |
| Frame units | `frame.meta.get(kEnergy)` [eV] matches Ha·`kHartreeInEv` |

## Performance

| Check | Criterion |
|-------|-----------|
| Wall ratio | \(T_{\mathrm{atv}}/T_{\mathrm{cp2k}} \le 10\) on fair twin, **or** `perf_exempt` with reason in REFERENCE |

## Build / test (scoped)

```bash
cmake --build build -j$(nproc) --target atomiverse test_cpu_gpw_uks reg_gpw-uks-multigrid
ctest --test-dir build -R '^(test_cpu_gpw_uks)$' --output-on-failure
cmake --build build --target run_reg_gpw-uks-multigrid
```

## Done when

All checks pass on `wf/gpw-uks-multigrid`; `acceptance_met=true`; green dual gate.
