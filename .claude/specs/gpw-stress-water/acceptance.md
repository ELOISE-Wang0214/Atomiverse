# Acceptance — `gpw-stress-water`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/gpw/gpw-stress-water/REFERENCE.md` with:
   - `oracle: cp2k` (or documented self-consistency fallback)
   - `oracle_version: <string>` when CP2K used
   - Geometry, basis/PP, XC, cutoff, cell, BC
   - Voigt stress table [Ha/Bohr³] and optional \(E_\mathrm{tot}\) [Ha]
   - Absolute tolerances stated
   - **No** live CP2K at ctest time

2. **Public API** (no Config bag / no factory):
   - `gpw_stress(atoms, method, h)` declared on public `gpw_solver.h`
   - Returns Voigt 6-vector [Ha/Bohr³]; orthorhombic v1: off-diag = 0
   - `GPWEvaluator::compute_stress` public member; writes `kStress` [eV/Å³]
   - Optional `last_stress_ha_bohr3()` for pin path

3. **Unit test** `tests/cpu/gpw/test_cpu_gpw_stress.cpp`:
   - Finite diagonal stress components on a cheap closed-shell system
   - Off-diagonal components abs ≤ 1e-12 (v1 zeros) or soft symmetry bound
   - No full water CP2K e2e in ctest

4. **E2E regression** `regressions/gpw/gpw-stress-water/reg_gpw-stress-water.cpp`:
   - Public headers only; links `atomiverse`
   - Same water geometry as forces twin
   - `GPWEvaluator` with `compute_stress=true`, Evaluate run
   - Assert max_k |σ_atv[k] − σ_cp2k[k]| ≤ **tol** for diagonal k∈{0,1,2}
   - Soft check off-diag |σ_cp2k| and ATV==0
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU preferred
     (document if strain-FD needs longer; still ≤ 180 s hard cap)

5. **Default tol** (product twin, may tighten if pin allows):
   - Diagonal: **5×10⁻⁶ Ha/Bohr³** max abs component error vs CP2K
     **NUMERICAL** stress (same-method FD twin). Do **not** pin against
     ANALYTICAL on nonperiodic MT (gross mismatch documented in REFERENCE).
   - Energy soft: |ΔE| ≤ **1×10⁻³ Ha** vs pin (same as forces twin)

6. **CMake**: `atv_add_regression(gpw gpw-stress-water)` under `ATV_BUILD_GPW`;
   unit test registered when GPW+tests ON.

7. **Build + PASS**:
   ```bash
   bash -lc '
     module load buildenv-gcc/2026.03-mpich
     export OPENBLAS_ROOT="$OPENBLASDIR"
     export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
     export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
     cmake -B build -DCMAKE_BUILD_TYPE=Release \
       -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON \
       -DATV_BUILD_GPW=ON -DATV_BUILD_REGRESSIONS=ON \
       -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS
     cmake --build build --target reg_gpw-stress-water test_cpu_gpw_stress -j$(nproc)
     ctest --test-dir build -R gpw_stress --output-on-failure
     ./build/regressions/reg_gpw-stress-water
   '
   ```

## Performance (dual gate part 2)

8. Fair wall vs the **same** CP2K binary on the twin input; record
   \(r = T_{\mathrm{atv}}/T_{\mathrm{cp2k}}\) in REFERENCE.md.
   Require \(r \le 10\) **or** `perf_exempt=true` with reason if no fair twin
   can be run offline.

## Process

9. Spec written **before** production code on `wf/gpw-stress-water`.
10. Commit(s) on `wf/gpw-stress-water` only; never `dev`; never push.
11. No public Config bags / factories introduced.
12. No live third-party process from ctest / regression binary.
