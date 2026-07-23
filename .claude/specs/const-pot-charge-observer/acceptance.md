# Acceptance — `const-pot-charge-observer`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/const-pot-charge-observer/REFERENCE.md`:
   - `oracle: lammps`
   - `oracle_version:` string matching the binary used (e.g. LAMMPS 30 Mar 2026
     development build).
   - Geometry / charges for the parallel-plate toy (type-1 left wall, type-2
     right wall, \(N_{\mathrm{wall}}\) atoms each, per-atom \(q\)).
   - Pin values \(Q_L^{\mathrm{pin}}\), \(Q_R^{\mathrm{pin}}\) from LAMMPS
     `compute reduce sum` on group charge (input twin under
     `tools/metalwalls_oracle/inputs/lammps/`).
   - Absolute tolerance with rationale (machine-eps scale for exact sum of
     fixed charges).
   - Perf note: measured ratio or `perf_exempt` reason.

2. **Public API** `cpu::ElectrodeChargeObserver`:
   - Constructor takes electrode type id list (no public `*Config`, no
     `make_*` / `build_*`).
   - `Observer` kind; default bind `md::Phase::end_of_step`.
   - `total_charge(e)` returns sum of `kAtom_charge` over owned atoms with
     `kAtom_type == electrode_types[e]` (after optional MPI sum).
   - Does not modify charges or coordinates.
   - Writes `frame.meta.get(kElectrodeCharge)` [e] for the first ≤3 electrodes.

3. **Unit test** `tests/cpu/observe/test_electrode_charge_observer.cpp`:
   - Synthetic Frame, fixed charges, two electrode types.
   - Asserts exact \(Q\) match (tol ≤ 1e-12).
   - Asserts non-electrode atoms do not contribute.
   - Single-function / method focus; no full MD e2e in `tests/`.

4. **E2E regression** `regressions/md/const-pot-charge-observer/`:
   - Out-of-library consumer (public headers + link `atomiverse` only).
   - Asserts \(Q_L, Q_R\) vs offline pin within stated tol.
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU.
   - CMake: `atv_add_regression(md const-pot-charge-observer)` under
     `ATV_BUILD_MD`.

5. **MD-only build PASS**:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
     -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
   ```
   Unit test + `reg_const-pot-charge-observer` PASS.

## Performance (dual gate part 2)

6. If a fair wall twin vs the same LAMMPS binary is recorded:
   \(r = T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\).
7. Else `perf_exempt=true` with reason (e.g. micro-observer O(N) sum; LAMMPS
   wall dominated by setup for `run 0` — not a fair force-loop twin).
   Green allowed under exempt.

## Process

8. Spec written **before** production code on `wf/const-pot-charge-observer`.
9. Commit(s) on that branch only; no push; never commit on `dev`.
10. No public `*Config` bag; no factory free functions; no new Driver category.
