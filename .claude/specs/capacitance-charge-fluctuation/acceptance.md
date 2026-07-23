# Acceptance — `capacitance-charge-fluctuation`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/capacitance-charge-fluctuation/REFERENCE.md`:
   - `oracle: literature` (MetalWalls / Scalfi–Salanne / Limmer charge-fluctuation
     \(C\); no LAMMPS twin for this analysis; no MW binary on host).
   - `oracle_version:` citation string for the formula source.
   - Documented \(T\) [K], optional \(\Delta\Psi\) [V], synthetic \(Q\) sample
     series [e], closed-form \(\langle Q\rangle\), \(\mathrm{Var}(Q)\),
     \(C_{\mathrm{diff}}\), \(C_{\mathrm{int}}\).
   - Absolute tolerances with rationale (exact rational moments + \(k_{\mathrm{B}}\)
     CODATA constant).
   - Perf note: `perf_exempt` reason (no oracle twin binary for fair wall).

2. **Public API** `cpu::CapacitanceChargeFluctuation`:
   - Constructor: electrode type id list (non-empty), temperature \(T\) [K]
     (`T > 0`), optional \(\Delta\Psi\) [V] (default 0 → integral \(C\) not
     defined / NaN or skipped).
   - No public `*Config`, no `make_*` / `build_*`.
   - `Observer` kind; default bind `md::Phase::end_of_step`.
   - Each `launch`: sum owned \(Q\) for electrode type 0 (MPI reduce), update
     running moments; publish members + meta.
   - Does not modify charges or coordinates.
   - Members: `n_samples()`, `mean_charge()` [e], `charge_variance()` [e²],
     `differential_capacitance()` [e/V], `integral_capacitance()` [e/V]
     (when \(\Delta\Psi \ne 0\)).

3. **Unit test** `tests/cpu/observe/test_capacitance_charge_fluctuation.cpp`:
   - Synthetic Frame, known charge rewrites across launches.
   - Asserts \(C_{\mathrm{diff}}\) and \(C_{\mathrm{int}}\) vs closed form
     (tol ≤ 1e-12 relative/absolute as appropriate).
   - Asserts empty types / non-positive \(T\) throw.
   - No full MD e2e in `tests/`.

4. **E2E regression** `regressions/md/capacitance-charge-fluctuation/`:
   - Out-of-library consumer (public headers + link `atomiverse` only).
   - Asserts observables vs offline pin within stated tol.
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU.
   - CMake: `atv_add_regression(md capacitance-charge-fluctuation)` under
     `ATV_BUILD_MD`.

5. **MD-only build PASS**:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
     -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
   ```
   Unit test + `reg_capacitance-charge-fluctuation` PASS.

## Performance (dual gate part 2)

6. Fair wall vs same oracle binary if available: \(r \le 10\).
7. Else `perf_exempt=true` with reason (analysis micro-kernel; no MetalWalls /
   LAMMPS twin for charge-fluctuation \(C\) on this host). Green allowed under
   exempt.

## Process

8. Spec written **before** production code on `wf/capacitance-charge-fluctuation`.
9. Commit(s) on that branch only; no push; never commit on `dev`.
10. No public `*Config` bag; no factory free functions; no new Driver category;
    no charge solve inside `KSpace*`.
