# Acceptance — `constant-force-pull-ion`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/constant-force-pull-ion/REFERENCE.md` with:
   - `oracle: lammps`
   - `oracle_version: <string>` (LAMMPS binary used to produce the pin)
   - Geometry, type selection, \(\mathbf{F}^{\mathrm{add}}\) [eV/Å], box [Å], mass, dt, steps
   - Per-atom added forces and bias energy \(U\) [eV] at `run 0`
   - Short NVE final unwrapped \(z\) of the pulled ion and work proxy \(W=F_z\Delta z\)
   - Absolute tolerances stated
   - **No** live LAMMPS at ctest time

2. **Public component** `ATV::cpu::md::FixAddForce`:
   - Constructor `(Fx, Fy, Fz, type_id)`; public members `Fx, Fy, Fz, type_id`
   - Inherits `cpu::Bias` (not a new Driver category)
   - No public `*Config` bag; no `make_*` / `build_*` factory
   - `launch`: selected type gets \(F \mathrel{+}= F^{\mathrm{add}}\),
     \(U \mathrel{+}= -F^{\mathrm{add}}\cdot r^{\mathrm{u}}\) into Frame
   - Default MD binding at `md::Phase::post_force` (LAMMPS POST_FORCE mask)

3. **Unit test** `tests/cpu/dynamics/test_fix_addforce.cpp`:
   - Force added only to matching type
   - Energy \(U = -\sum F\cdot r^{\mathrm{u}}\) closed form
   - Image-unwrapped energy when `kAtom_image` non-zero (ortho box)

4. **E2E regression** `regressions/md/constant-force-pull-ion/reg_constant-force-pull-ion.cpp`:
   - Public headers only; links `atomiverse`
   - Builds the REFERENCE geometry on a `cpu::Frame`
   - **Gate A:** `cpu::Evaluate` + `FixAddForce` → forces + \(U\) vs pin
   - **Gate B:** short NVE MD with `FixAddForce` → final \(z\) of pulled ion and
     work proxy \(W\) vs pin
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU

5. **CMake**: `atv_add_regression(md constant-force-pull-ion)` under `ATV_BUILD_MD`;
   unit test registered; `fix_addforce.cpp` in `ATOMIVERSE_SRCS` when MD on.

6. **MD-only build + PASS**:
   ```bash
   bash -lc '
     module load buildenv-gcc/2026.03-mpich
     export OPENBLAS_ROOT="$OPENBLASDIR"
     export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
     export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
     cmake -B build -DCMAKE_BUILD_TYPE=Release \
       -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
       -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
       -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
       -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
     cmake --build build --target reg_constant-force-pull-ion test_cpu_fix_addforce -j$(nproc)
     ./build/regressions/reg_constant-force-pull-ion
     ctest --test-dir build -R test_cpu_fix_addforce --output-on-failure
   '
   ```
   Forbidden: HF/PW/GPW/CUDA/Torch **ON**.

## Performance (dual gate part 2)

7. Fair wall vs the **same** LAMMPS binary on a comparable pure-addforce workload
   (pair zero / inert pair, NVE, constant addforce on a type group); record
   \(r = T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) in REFERENCE.md.
   Require \(r \le 10\) **or** `perf_exempt=true` with reason if no fair twin
   can be run offline.

## Process

8. Spec written **before** production code on `wf/constant-force-pull-ion`.
9. Commit(s) on `wf/constant-force-pull-ion` only; never `dev`; never push.
10. No public Config bags / factories introduced.
