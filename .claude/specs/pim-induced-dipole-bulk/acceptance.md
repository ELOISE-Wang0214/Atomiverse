# Acceptance — `pim-induced-dipole-bulk`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/pim-induced-dipole-bulk/REFERENCE.md` with:
   - `oracle: metalwalls` (MW-only classical SCF PIM; not LAMMPS Drude)
   - `oracle_version: <string>` documenting formula source + that no MW
     binary was available (or MW version if offline run used)
   - Geometry, charges, \(\alpha\), Thole \(a\), cutoff, box
   - Hard-coded \(E_{\mathrm{pol}}\) [eV] and per-atom forces [eV/Å]
   - Absolute tolerances stated
   - **No** live MetalWalls / LAMMPS at ctest time

2. **Public component** `ATV::cpu::PimInducedDipole`:
   - Header `include/atomiverse/cpu/potentials/pim_induced_dipole.h`
   - Constructor `(cutoff, alpha_by_type, thole_a=0.39, pol_tol=1e-10,
     pol_max_iter=100)`; optional public `tt_b` for Tang–Toennies
   - Inherits `cpu::Potential` (Component on MD/Evaluate — **not** a new
     Driver category; charge solve **not** inside `KSpace*`)
   - No public `*Config` bag; no `make_*` / `build_*` factory
   - Self-consistent \(\mu=\alpha E^{\mathrm{tot}}\) with Thole-damped
     \(\mathbf{T}\); \(E_{\mathrm{pol}}=-\frac12\sum\mu\cdot E^{\mathrm{perm}}\)
     [eV]; Hellmann–Feynman forces (q–μ + μ–μ)
   - Frame I/O only (charges/types/positions → energy/forces)

3. **Unit test** `tests/cpu/potentials/test_pim_induced_dipole.cpp`
   (ctest `test_cpu_pim_induced_dipole`):
   - Single site: \(\mu=\alpha E\) under uniform permanent field proxy
     (two-charge limit or injected field path)
   - Two-site SCC energy vs closed-form / high-tol reference
   - Thole \(a>0\) changes \(E_{\mathrm{pol}}\) vs undamped
   - TT (`tt_b>0`) path runs and is finite
   - Analytic forces match central FD of \(E_{\mathrm{pol}}\) within tol

4. **E2E regression** `regressions/md/pim-induced-dipole-bulk/`:
   - Public headers + link `atomiverse` only
   - Builds REFERENCE geometry on `cpu::Frame`
   - Runs via `cpu::Evaluate` + `PimInducedDipole`
   - Asserts \(E_{\mathrm{pol}}\) and forces vs offline pin within tol
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU

5. **CMake**: `atv_add_regression(md pim-induced-dipole-bulk)` under
   `ATV_BUILD_MD`; unit test registered in `tests/CMakeLists.txt`.

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
     cmake --build build --target reg_pim-induced-dipole-bulk \
       test_cpu_pim_induced_dipole -j$(nproc)
     ./build/regressions/reg_pim-induced-dipole-bulk
     ./build/tests/test_cpu_pim_induced_dipole
   '
   ```
   Forbidden: HF/PW/GPW/CUDA/Torch **ON**.

## Performance (dual gate part 2)

7. Fair wall vs the **same** MetalWalls binary on a comparable bulk PIM
   workload if available; record
   \(r=T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) in REFERENCE.md.
   Require \(r\le 10\) **or** `perf_exempt=true` with reason (expected:
   no MW binary on agent host → exempt).

## Process

8. Spec written **before** production code on `wf/pim-induced-dipole-bulk`.
9. Commit(s) on `wf/pim-induced-dipole-bulk` only; never `dev`; never push.
10. No public Config bags / factories introduced.
11. Does **not** implement electrode-joint polarization (E5).
