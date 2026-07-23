# Acceptance — `external-efield-md`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/external-efield-md/REFERENCE.md` with:
   - `oracle: lammps`
   - `oracle_version: <string>` (LAMMPS binary used to produce the pin)
   - Geometry, charges, \(\mathbf{E}\) [V/Å], box [Å]
   - Per-atom forces \(F_i = q_i E\) and total field energy \(U\) [eV]
   - Absolute tolerances stated
   - **No** live LAMMPS at ctest time

2. **Public component** `ATV::cpu::md::FixEfield`:
   - Constructor `(Ex, Ey, Ez)` in V/Å; public members `Ex, Ey, Ez`
   - Inherits `cpu::Bias` (not a new Driver category)
   - No public `*Config` bag; no `make_*` / `build_*` factory
   - `launch`: \(F \mathrel{+}= qE\), \(U \mathrel{+}= -q E\cdot r^{\mathrm{u}}\) into Frame
   - Default MD binding at `md::Phase::post_force` (LAMMPS fix efield mask)

3. **Unit test** `tests/cpu/dynamics/test_fix_efield.cpp`:
   - Single ion / multi-ion \(F=qE\) and energy vs closed form
   - Image-unwrapped energy when `kAtom_image` non-zero (ortho box)

4. **E2E regression** `regressions/md/external-efield-md/reg_external-efield-md.cpp`:
   - Public headers only; links `atomiverse`
   - Builds the REFERENCE geometry on a `cpu::Frame`
   - Runs via `cpu::Evaluate` + `on(evaluate::Phase::evaluate).add(FixEfield{…})`
   - Asserts per-atom forces and total energy vs offline pin within tol
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU

5. **CMake**: `atv_add_regression(md external-efield-md)` under `ATV_BUILD_MD`.

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
     cmake --build build --target reg_external-efield-md -j$(nproc)
     ./build/regressions/reg_external-efield-md
   '
   ```
   Forbidden: HF/PW/GPW/CUDA/Torch **ON**.

## Performance (dual gate part 2)

7. Fair wall vs the **same** LAMMPS binary on a comparable pure-efield (pair zero)
   workload; record \(r = T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) in REFERENCE.md.
   Require \(r \le 10\) **or** `perf_exempt=true` with reason if no fair twin
   can be run offline.

## Process

8. Spec written **before** production code on `wf/external-efield-md`.
9. Commit(s) on `wf/external-efield-md` only; never `dev`; never push.
10. No public Config bags / factories introduced.
