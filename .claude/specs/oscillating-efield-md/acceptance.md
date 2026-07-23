# Acceptance — `oscillating-efield-md`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/oscillating-efield-md/REFERENCE.md` with:
   - `oracle: lammps`
   - `oracle_version: <string>` (LAMMPS binary used to produce the pin)
   - Geometry, charges, \(\mathbf{E}_0\) [V/Å], period \(T\), \(\Delta t\), box [Å]
   - Multi-step table: at least steps covering max / zero-crossing / min of
     \(\cos(2\pi t/T)\) with per-atom \(F_z=q E_z(t)\) and \(U(t)\) [eV]
   - Absolute tolerances stated
   - **No** live LAMMPS at ctest / regression runtime

2. **Public component** `ATV::cpu::md::FixOscillatingEfield`:
   - Constructor `(Ex0, Ey0, Ez0, period, dt, phase=0)` —
     amplitudes [V/Å], period/dt [fs], phase [rad]
   - Public members; **no** public `*Config` bag; no `make_*` / `build_*`
   - Inherits `cpu::Bias` (not a new Driver category)
   - `launch`: \(\mathbf{E}(t)=\mathbf{E}_0\cos(2\pi\,n\Delta t/T+\varphi)\),
     \(F\mathrel{+}=qE\), \(U\mathrel{+}=-q E\cdot r^{\mathrm{u}}\) into Frame
   - Default MD binding at `md::Phase::post_force`
   - Helpers `field_x/y/z(step)` for the instantaneous field

3. **Unit test** `tests/cpu/dynamics/test_fix_oscillating_efield.cpp`:
   - Schedule: \(E(0)=E_0\), \(E(T/4\Delta t)=0\), \(E(T/2\Delta t)=-E_0\) (within tol)
   - Multi-ion \(F=qE(t)\) and energy vs closed form at a non-trivial step

4. **E2E regression** `regressions/md/oscillating-efield-md/reg_oscillating-efield-md.cpp`:
   - Public headers only; links `atomiverse`
   - REFERENCE geometry on `cpu::Frame`
   - Evaluates `FixOscillatingEfield` at pin steps via `Component::launch`
     (`ctx.step` = MD step index)
   - Asserts per-atom forces **and** total field energy vs offline pin (dual gate)
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU

5. **CMake**: `atv_add_regression(md oscillating-efield-md)` under `ATV_BUILD_MD`;
   unit test registered; source under `ATV_BUILD_MD` atomiverse src list.

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
     cmake --build build --target reg_oscillating-efield-md test_cpu_fix_oscillating_efield -j$(nproc)
     ./build/regressions/reg_oscillating-efield-md
     ctest --test-dir build -R test_cpu_fix_oscillating_efield --output-on-failure
   '
   ```
   Forbidden: HF/PW/GPW/CUDA/Torch **ON**.

## Performance (dual gate part 2)

7. Fair wall vs the **same** LAMMPS binary on a comparable pure-efield
   oscillating workload (pair zero / inert, variable efield + nve); record
   \(r = T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) in REFERENCE.md.
   Require \(r \le 10\) **or** `perf_exempt=true` with reason if no fair twin
   can be run offline.

## Process

8. Spec written **before** production code on `wf/oscillating-efield-md`.
9. Commit(s) on `wf/oscillating-efield-md` only; never `dev`; never push.
10. No public Config bags / factories introduced.
