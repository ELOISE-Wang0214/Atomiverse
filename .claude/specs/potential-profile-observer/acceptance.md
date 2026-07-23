# Acceptance — `potential-profile-observer`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/md/potential-profile-observer/REFERENCE.md`:
   - `oracle: literature` (planar Poisson / parallel-plate capacitor;
     MetalWalls potential-profile analysis lineage; no LAMMPS twin for this
     named observer; no MW binary on host).
   - `oracle_version:` citation string (Poisson planar continuum + MW JOSS
     analysis tools).
   - Documented geometry: \(L_x=L_y=10\) Å, \(A=100\) Å², \(z\in[0,20)\) Å,
     \(n_{\mathrm{bins}}=20\), left sheet \(Q=+1\,e\) at bin-5 center,
     right sheet \(Q=-1\,e\) at bin-14 center.
   - Pin values: mid-gap \(E_z\), right-edge \(\Psi\), selected bin
     \(\Psi_i\), continuum \(E=4\pi k_{\mathrm{e}}\sigma\).
   - Absolute tolerances with rationale.
   - Perf note: `perf_exempt` reason.

2. **Public API** `cpu::PotentialProfile`:
   - Constructor: `n_bins` (\(>0\)), `z_min`, `z_max` (\(z_{\max}>z_{\min}\)).
   - No public `*Config`, no `make_*` / `build_*`.
   - `Observer` kind; default bind `md::Phase::end_of_step`.
   - Each `launch`: bin owned charges → \(\rho\) → planar Poisson \(E,\Psi\);
     publish members + `kPotentialDrop` meta.
   - Does not modify charges or coordinates.
   - Members: `n_bins()`, `dz()`, `rho()`, `field()`, `psi()`, `z_center(i)`.

3. **Unit test** `tests/cpu/observe/test_potential_profile.cpp`:
   - Synthetic parallel-plate Frame; assert mid-gap \(E\) and key \(\Psi\)
     vs closed form (tol ≤ 1e-12 absolute on doubles matching the discrete
     scheme).
   - Asserts invalid ctor args throw.
   - No full MD e2e in `tests/`.

4. **E2E regression** `regressions/md/potential-profile-observer/`:
   - Out-of-library consumer (public headers + link `atomiverse` only).
   - Asserts observables vs offline pin within stated tol.
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU.
   - CMake: `atv_add_regression(md potential-profile-observer)` under
     `ATV_BUILD_MD`.

5. **MD-only build PASS**:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_MD=ON -DATV_BUILD_KSPACE=ON \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_HF=OFF -DATV_BUILD_PW=OFF -DATV_BUILD_GPW=OFF \
     -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF -DBLA_VENDOR=OpenBLAS
   ```
   Unit test + `reg_potential-profile-observer` PASS.

## Performance (dual gate part 2)

6. Fair wall vs same oracle binary if available: \(r \le 10\).
7. Else `perf_exempt=true` with reason (planar Poisson bin micro-kernel;
   no MetalWalls / LAMMPS twin for this named observer on this host).
   Green allowed under exempt.

## Process

8. Spec written **before** production code on `wf/potential-profile-observer`.
9. Commit(s) on that branch only; no push; never commit on `dev`.
10. No public `*Config` bag; no factory free functions; no new Driver category;
    no charge solve inside `KSpace*`.
