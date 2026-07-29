# Acceptance — `dftb2-periodic-forces`

## Correctness

1. **Analytic vs re-SCC FD (C-chain):** on the off-equilibrium 1-D C dimer cell
   (`C` at \(0\) and \(1.20\) Å, \(L_x=2.8\), \(L_y=L_z=12\); MioSketch,
   `scc_tol=1e-10`), \(|F^{\mathrm{ana}}_{0,x} - F^{\mathrm{FD}}_{0,x}|
   \le 1\times 10^{-5}\) eV/Å where FD is central re-SCC of total \(E\) at
   \(h=10^{-4}\) Å. Same for atom 1.
2. **Nonzero force:** \(|F|^2 > 10^{-8}\) at that off-eq geometry (not the
   equilibrium \(r=1.4\) pin of `dftb2-periodic-gamma`).
3. **Action-reaction:** \(|\sum_i F_i| \le 10^{-6}\) eV/Å per Cartesian.
4. **PBC matters:** energy of the dense cell differs from the same dimer in a
   20 Å cubic box by \(> 0.05\) eV (lattice images in range).
5. **Molecular limit preserved:** large-box H2O analytic force vs re-SCC FD
   still passes existing `test_cpu_dftb2` / `reg_dftb2-analytic-forces` gates
   (min-image ≡ lattice sum when \(L \gg r_{\mathrm{cut}}\)).
6. **Polar γ lattice (unit test):** water in a tight orthorhombic box with
   images inside \(r_{\mathrm{cut}}\) has analytic \(F\) matching re-SCC FD
   within \(1\times 10^{-5}\) eV/Å on at least one gated component.

## Perf

```text
perf_exempt: true
perf_exempt_reason: no offline DFTB+ binary on host; formula-level pin only
```

Optional: if DFTB+ wall available, \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\).

## Build / test (scoped)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
  -DATV_BUILD_SEMIEMPIRICAL=ON -DATV_BUILD_EVALUATE=ON \
  -DATV_BUILD_MD=OFF -DATV_BUILD_MC=OFF -DATV_BUILD_MINIMIZER=OFF \
  -DATV_BUILD_VIB=OFF -DATV_BUILD_PATH=OFF -DATV_BUILD_POTENTIALS=OFF \
  -DATV_BUILD_NEIGHBOR=OFF -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF \
  -DATV_BUILD_GPW=OFF -DATV_BUILD_KSPACE=OFF -DATV_BUILD_CUDA=OFF \
  -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
cmake --build build -j$(nproc) --target atomiverse test_cpu_dftb2 reg_dftb2-periodic-forces
ctest --test-dir build -R '^(test_cpu_dftb2)$' --output-on-failure
cmake --build build --target run_reg_dftb2-periodic-forces
```

## Artifacts

| Path | Role |
| ---- | ---- |
| `regressions/semiempirical/dftb2-periodic-forces/` | e2e gate + REFERENCE |
| `tools/tb_oracle/inputs/dftbplus/dftb2-periodic-forces/` | offline twin / note |
| `tests/cpu/semiempirical/test_dftb2.cpp` | unit FD gates |
| `.claude/specs/dftb2-periodic-forces/` | this spec (delete at done) |

## Done when

- `acceptance_met = true`
- `reg_dftb2-periodic-forces` PASS
- `test_cpu_dftb2` PASS (incl. new PBC force tests)
- commit on `wf/dftb2-periodic-forces`
