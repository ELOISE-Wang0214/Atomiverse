# Acceptance — `fd-phonon-supercell`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin** at `regressions/pw/fd-phonon-supercell/REFERENCE.md` with:
   - Oracle note: CP2K `VIBRATIONAL_ANALYSIS` if run; else SW literature + acoustic sum rule
   - Geometry, supercell factors, force engine, \(h\) [Å]
   - Acoustic tol and optical band stated
   - **No** live CP2K at ctest time

2. **Public API**:
   - `cpu::vib::expand_supercell(const Frame&, na, nb, nc) → Frame`
   - `cpu::vib::Phonon` with ctor `(na, nb, nc)` and public members `na, nb, nc`
   - `cpu::vib::PhononResult` with `frequencies` [cm⁻¹], `n_atoms`, `na/nb/nc`, `hessian`
   - No public `*Config` bag; no `make_*` / `build_*` factory
   - Force evaluator bound on `vib::Phase::force` (composition)

3. **Unit test** `tests/cpu/vibration/test_phonon_supercell.cpp`:
   - `expand_supercell` atom count, cell scaling, position images
   - Acoustic near-zero on SW (or zero-force) tiny crystal when forces are consistent
   - `Phonon` rejects missing force kind at `build()`

4. **E2E regression** `regressions/pw/fd-phonon-supercell/reg_fd-phonon-supercell.cpp`:
   - Public headers only; links `atomiverse`
   - Diamond Si + SW force twin per REFERENCE
   - Asserts mode count, acoustic tol, optical band
   - Exit non-zero with got/want/tol on failure; ≤ 60 s single-rank CPU

5. **CMake**: sources under `ATV_BUILD_VIB`; `atv_add_regression(pw fd-phonon-supercell)` when MD+VIB available for SW twin.

6. **Build + PASS** with project toolchain (`buildenv-gcc/2026.03-mpich`, OpenBLAS, …) and
   `-DATV_BUILD_TESTS=ON -DATV_BUILD_PW=ON -DATV_BUILD_GPW=ON -DATV_BUILD_REGRESSIONS=ON -DATV_BUILD_CUDA=OFF`.

## Performance (dual gate part 2)

7. Fair wall vs CP2K on a comparable VIBRATIONAL_ANALYSIS twin if binary available offline;
   \(r = T_\mathrm{atv}/T_\mathrm{cp2k} \le 10\). Else `perf_exempt=true` with reason in REFERENCE.md.

## Process

8. Spec written **before** production code on `wf/fd-phonon-supercell`.
9. Commit(s) on `wf/fd-phonon-supercell` only; never `dev`; never push.
10. No public Config bags / factories introduced.
