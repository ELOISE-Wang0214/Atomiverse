# Spec: `gpw-water-pbe-mt`

| Field | Value |
| ----- | ----- |
| Slug | `gpw-water-pbe-mt` |
| Branch | `wf/gpw-water-pbe-mt` |
| Domain | GPW (Γ-point closed-shell RKS) |
| Observable | total energy \(E_\mathrm{tot}\) [Hartree] |
| Geometry | isolated H₂O (CP2K getting-started `h2o.inp`) |
| Basis / PP | DZVP-MOLOPT-GTH + GTH-PBE-q6 (O) / GTH-PBE-q1 (H) |
| XC | PBE (GGA) |
| Electrostatics | Martyna–Tuckerman isolated Poisson (`nonperiodic = true`) |
| Cutoff | 400 Ry |
| Box | 10 Å cubic |
| Public API | `cpu::gpw::GPWEvaluator` + `cpu::single_point` (no public `*Config` bag) |
| Regression | `regressions/gpw/gpw-water-pbe-mt/` |
| Queue row | `docs/developer/cp2k-method-queue.md` `## in_progress` `[~] gpw-water-pbe-mt` |

## Goal

Ship the **first** offline-CP2K-pinned end-to-end GPW regression: isolated water
total energy via the public GPW evaluator, matching the CP2K Quickstep
getting-started twin within the stated tolerance. Product narrative is tutorial
ch09; this slug creates `regressions/gpw/` (none exists on parent `dev` yet).

## Non-goals

- Analytic forces / stress / UKS / k-points / D3 / BOMD (other queue slugs).
- Growing public `*Config` bags or `make_*` / `build_*` factories.
- Live CP2K (or any third-party) at ctest time.
- MPI multi-rank scaling.
- Engine redesign: reuse existing `GPWEvaluator` / `GpwDft` path.

## Tasks

1. **Offline CP2K pin** in `regressions/gpw/gpw-water-pbe-mt/REFERENCE.md`:
   cite geometry, basis/PP, cutoff, XC, MT Poisson, \(E_\mathrm{tot}\) Ha, and
   provenance (CP2K getting-started `h2o.inp` pin already used by
   `examples/09_gpw.cpp` / ch09). Optional twin input under
   `tools/cp2k_oracle/inputs/` (not required at ctest).
2. **E2E regression** `regressions/gpw/gpw-water-pbe-mt/reg_gpw-water-pbe-mt.cpp`:
   public API only; Frame geometry; `GPWEvaluator` (PBE, MT, 400 Ry);
   `cpu::single_point`; assert SCF converged + \(E_\mathrm{tot}\) vs pin.
3. **Wire CMake**: `atv_add_regression(gpw gpw-water-pbe-mt)` behind
   `ATV_BUILD_GPW` (implies HF+PW+FFTW).
4. **Build + run** regression until PASS (`ATV_BUILD_REGRESSIONS=ON`,
   `ATV_BUILD_GPW=ON`, `ATV_BUILD_PW=ON`, `ATV_BUILD_HF=ON`).
5. **Perf gate**: fair wall vs CP2K if binary available; else
   `perf_exempt=true` with reason recorded in REFERENCE.md.
6. Commit on `wf/gpw-water-pbe-mt` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md` (binding).
