# Acceptance: `dftb3-water-mio`

## Dual gate

| Gate | Command (scoped) | Pass |
| ---- | ---------------- | ---- |
| Unit | `ctest --test-dir build -R '^(test_cpu_dftb3)$' --output-on-failure` | all tests PASS |
| Regression | `cmake --build build --target run_reg_dftb3-water-mio` | exit 0, prints PASS |

## Build (SE-only, once per tree)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
  -DATV_BUILD_SEMIEMPIRICAL=ON -DATV_BUILD_EVALUATE=ON \
  -DATV_BUILD_MD=OFF -DATV_BUILD_MC=OFF -DATV_BUILD_MINIMIZER=OFF \
  -DATV_BUILD_VIB=OFF -DATV_BUILD_PATH=OFF \
  -DATV_BUILD_POTENTIALS=OFF -DATV_BUILD_NEIGHBOR=OFF \
  -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
  -DATV_BUILD_KSPACE=OFF -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
  -DBLA_VENDOR=OpenBLAS
cmake --build build -j$(nproc) --target atomiverse test_cpu_dftb3 reg_dftb3-water-mio
```

## Checklist

- [ ] `DFTB3` public type, members only (no Config / factory)
- [ ] Hubbard \(U^{\mathrm{d}}\) on static element rows
- [ ] \(E^{(3)}\) and third-order shift in SCC loop
- [ ] `REFERENCE.md`: `oracle: dftb+`, `oracle_version`, eV
- [ ] Unit + regression green
- [ ] Build artifacts wiped after gate

## Status

`acceptance_met`: **true** (2026-07-27) — `test_cpu_dftb3` PASS;
`run_reg_dftb3-water-mio` PASS with pin \(E=-84.561719022489\) eV.
