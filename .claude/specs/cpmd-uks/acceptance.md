# Acceptance — `cpmd-uks`

## Correctness (dual gate, correctness half)

| # | Gate | How |
| - | ---- | --- |
| 1 | `GPWEvaluator` with `cpmd=true` + `n_alpha`/`n_beta` builds (no UKS throw) | unit |
| 2 | Production CPMD steps report `last_result().iters == 0` | regression |
| 3 | \(E_{\mathrm{cons}}=T_n+T_e+E_{\mathrm{KS}}\) relative drift ≤ 0.10 over short NVE | regression |
| 4 | ⟨S²⟩ stays in triplet window [1.5, 2.5] for O₂ | regression |
| 5 | O–O bond stays in [0.9, 1.6] Å | regression |
| 6 | PE after warm-up < −100 eV; finite energies/coords | regression |
| 7 | Dual RATTLE residual unit tests still green; occ=1 force path unit-tested | unit |
| 8 | RKS CPMD path unchanged (existing regs still green if rebuilt) | no regression of API |

## Perf (dual gate, wall half)

| # | Gate |
| - | ---- |
| P1 | **perf_exempt** — NVE conservation envelope (same class as `cpmd-water`, `xl-bomd-uks`); no CP2K CPMD UKS wall twin required for green |

## Build / test commands (scoped)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
  -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON -DATV_BUILD_GPW=ON \
  -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS
cmake --build build -j$(nproc) --target atomiverse test_gpw_cpmd reg_cpmd-uks
ctest --test-dir build -R '^(test_gpw_cpmd)$' --output-on-failure
cmake --build build -j$(nproc) --target run_reg_cpmd-uks
```

## Done when

- All gates above pass on `wf/cpmd-uks`
- `regressions/dynamics/cpmd-uks/{reg_cpmd-uks.cpp,REFERENCE.md}` present
- Spec deleted only on merge close (harness)
