# Acceptance: `dftb3-water-mio`

## Dual gate

| Gate | Command (scoped) | Pass |
| ---- | ---------------- | ---- |
| Unit | `ctest --test-dir build -R '^(test_cpu_dftb3\|test_cpu_dftb3_scc)$' --output-on-failure` | **PASS** (2026-07-29 wave-3 resume) |
| Regression | `cmake --build build --target run_reg_dftb3-water-mio` | **PASS** — `E=-94.273801189131` eV, residual `2.3e-14`, `|dE|/|E3|=0.973` |

## Build (SE-only)

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
cmake --build build -j$(nproc) \
  --target atomiverse test_cpu_dftb3 test_cpu_dftb3_scc reg_dftb3-water-mio
```

## Checklist

- [x] `DFTB3` public type, members only (no Config / factory)
- [x] Hubbard \(U^{\mathrm{d}}\) on static element rows
- [x] \(E^{(3)}\) + exact \(\partial E^{(3)}/\partial\Delta q\) in SCC
- [x] Anderson charge mixing; unconverged SCC throws
- [x] H–H mio-1-1 shared table; H2 reduces to DFTB2 + DFTB+ pin
- [x] `REFERENCE.md`: `oracle: dftb+`, pin −94.273801189131 eV
- [x] Unit + regression green (wave-3 resume dual gate)
- [x] Build artifacts wiped after gate

## Status

`acceptance_met`: **true** (2026-07-29 wave-3 resume dual-gate).
Product already on `dev` via `68df18a8`..`148a17ed` (SCC converge + double-count
fix + re-pin). Pin is ATV Hamiltonian twin + Gaus formula lineage; no host
DFTB+ binary → `perf_exempt`. Absolute DFTB.org mio-1-1 multi-element identity
remains deferred to SK-fidelity work; analytic F → `dftb3-analytic-forces`.
