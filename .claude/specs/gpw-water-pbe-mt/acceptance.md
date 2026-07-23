# Acceptance — `gpw-water-pbe-mt`

Every item must be verified before `acceptance_met=true` / `status=green`.

## Correctness (dual gate part 1)

1. **Offline pin exists** at
   `regressions/gpw/gpw-water-pbe-mt/REFERENCE.md` with:
   - CP2K \(E_\mathrm{tot} = -17.21233714458516\) Ha (getting-started
     `h2o.inp` family: PBE, DZVP-MOLOPT-GTH, GTH-PBE-q6/q1, CUTOFF 400 Ry,
     PERIODIC NONE / PSOLVER MT, 10 Å cubic cell, geometry below).
   - Explicit provenance (not a golden of Atomiverse itself).
   - Stated absolute tolerance in Hartree and rationale.
2. **Geometry** (Å), identical to `examples/09_gpw.cpp` / CP2K tutorial:

   | Atom | x | y | z |
   | ---- | --: | --: | --: |
   | O | 5.0000 | 5.0000 | 5.0000 |
   | H | 5.7586 | 5.0000 | 5.5043 |
   | H | 4.2414 | 5.0000 | 5.5043 |

   Box: cubic \(L = 10\) Å.
3. **Public API regression** `reg_gpw-water-pbe-mt`:
   - Includes only public `atomiverse/…` headers; links `atomiverse`.
   - Constructs `cpu::gpw::Kind` with
     `GaussianBasis(Element, "DZVP-MOLOPT-GTH")` and
     `GthPotential(Element, "GTH-PBE-q6"|"GTH-PBE-q1")` — lookup constructors,
     no factories / public Config bags.
   - Runs `cpu::gpw::GPWEvaluator({O,H}, 400, 200, 1e-7, 0.3, XcKind::PBE,
     /*nonperiodic=*/true)` via `cpu::single_point(frame, eval)`.
   - Asserts `last_result().converged == true`.
   - Asserts
     \(|E_\mathrm{tot}^\mathrm{ATV} - (-17.21233714458516)| \le 1\times 10^{-3}\) Ha
     (known residual ~0.42 mHa on this stack; 1 mHa e2e band).
   - Exit non-zero with `got / want / tol` on failure; runs ≤ 60 s single CPU.
4. **CMake**: case registered under `ATV_BUILD_GPW` as
   `atv_add_regression(gpw gpw-water-pbe-mt)`.
5. **Build + run PASS** with:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release \
     -DATV_BUILD_TESTS=ON -DATV_BUILD_HF=ON -DATV_BUILD_PW=ON \
     -DATV_BUILD_GPW=ON -DATV_BUILD_REGRESSIONS=ON \
     -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS
   cmake --build build --target reg_gpw-water-pbe-mt -j$(nproc)
   ./build/regressions/reg_gpw-water-pbe-mt
   ```

## Performance (dual gate part 2)

6. If a CP2K binary is available offline: record
   \(r = T_\mathrm{atv}/T_\mathrm{cp2k}\) in REFERENCE.md; require \(r \le 10\).
7. If no CP2K binary: set `perf_exempt=true` in REFERENCE.md with reason
   (e.g. no `cp2k` module/binary on the build host). Green is allowed under
   exempt.

## Process

8. Spec written **before** production regression on `wf/gpw-water-pbe-mt`.
9. Commit(s) on `wf/gpw-water-pbe-mt` only; no push; never commit on `dev`.
10. No public `*Config` bag and no `make_*`/`build_*` factory introduced.
