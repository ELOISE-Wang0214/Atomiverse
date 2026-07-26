---
slug: esmm-cpmd-01-electronic-home
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: []
unblocks: [esmm-cpmd-02, esmm-cpmd-04, esmm-cpmd-05]
grilled: true
---

# Spec: `esmm-cpmd-01-electronic-home` — shared electronic-propagation layer

**Slug:** `esmm-cpmd-01-electronic-home`
**Domain:** electronic-state propagators (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)
**Kind:** pure refactor — zero behavior change; the whole suite stays green bit-for-bit.

## Summary

The electronic-state propagators are today parked inside method domains:
`Cpmd`/`CpmdCplx` and `XlBomd`/`XlBomdCplx` in `cpu/gpw/`, `ASPCPredictor`/
`ASPCPredictorCplx` in `cpu/scf/` (already consumed cross-domain by GPW — the
smell being fixed, not the precedent to extend). Wave 3 needs all five usable
from `MeanFieldMethod` without domain-crossing includes. Create the root infra
cluster **`cpu/electronic/`** (namespace `ATV::cpu::electronic`; "dynamics =
nuclear-DOF motion, electronic = electronic-DOF motion") and move them there.

Blocking dependency to untangle first: `src/cpu/gpw/cpmd.cpp` calls
`gpw::internal::ot_gradient` / `ot_orthonormalize` (+`_k`) from
`gpw/internal/ot.h`, whose `ot.cpp` compiles **only under `ATV_BUILD_GPW`** —
an `electronic/` cluster including `gpw/internal` would invert layering and
break HF-only builds. Extract the generic S-metric orbital primitives into
`cpu/linalg/orbital_tangent.{h,cpp}`, renamed to speak function not origin.

## Design

### New: `cpu/linalg/orbital_tangent.h` (namespace `ATV::cpu::linalg`)

```cpp
/// @brief Projected band-energy gradient g = (I − S C Cᵀ)(2·occ·H C) [Hartree].
RealMatrix band_gradient(const RealMatrix& H, const RealMatrix& S,
                         const RealMatrix& C, double occ = 2.0);
/// @brief S-tangent projection (I − C Cᵀ S) d.
RealMatrix s_tangent(const RealMatrix& C, const RealMatrix& S, const RealMatrix& d);
/// @brief Löwdin re-orthonormalization C ← C (CᵀSC)^{-1/2}.
RealMatrix lowdin_orthonormalize(const RealMatrix& C, const RealMatrix& S);
// + CplxMatrix twins band_gradient_k / s_tangent_k / lowdin_orthonormalize_k
```

Bodies move verbatim from `src/cpu/gpw/internal/ot.cpp` (`ot_gradient`,
`ot_tangent`, `ot_orthonormalize` + `_k` twins). `gpw/internal/ot.h` keeps the
OT-specific parameterization (`ot_parameterize`), preconditioners
(`OtPrecond`, `OtPrecondSingleInverse` + `K`), and the fixed-H minimizers; it
includes the linalg header and its ~10 call sites in `gpw_solver.cpp` /
`ot.cpp` / `tests/cpu/gpw/test_cpu_gpw_ot.cpp` update mechanically. **No
back-compat aliases.**

### New: `cpu/electronic/` (namespace `ATV::cpu::electronic`)

| File | Moved from | Classes |
| ---- | ---------- | ------- |
| `electronic/cpmd.h` + `src/cpu/electronic/cpmd.cpp` | `gpw/cpmd.{h → deleted}`, `src/cpu/gpw/cpmd.cpp` | `Cpmd`, `CpmdCplx` |
| `electronic/xl_bomd.h` + `src/cpu/electronic/xl_bomd.cpp` | `gpw/xl_bomd.h`, `src/cpu/gpw/xl_bomd.cpp` | `XlBomd`, `XlBomdCplx` |
| `electronic/aspc_predictor.h` + `src/cpu/electronic/aspc_predictor.cpp` | `scf/aspc_predictor.h`, `src/cpu/scf/aspc_predictor.cpp` | `ASPCPredictor`, `ASPCPredictorCplx` |

`cpmd.cpp` now includes `cpu/linalg/orbital_tangent.h` (not `gpw/internal/ot.h`).
Namespace moves to `electronic::` — a class physically in `cpu/electronic/`
carrying `gpw::` would be a permanent lie. Old headers are **deleted**, no
forwarding shims (house rule: no backward compatibility anywhere).

### Include-site migration (verified blast radius, 12 files)

`gpw_evaluator.h` (+`.cpp`) member types → `electronic::Cpmd` etc.;
`src/cpu/scf/mean_field_method.cpp`; `mean_field_method.h:828` forward decl →
qualified `electronic::ASPCPredictor`; `src/cpu/electronic/{cpmd,xl_bomd,
aspc_predictor}.cpp` themselves; `tests/cpu/gpw/test_gpw_cpmd.cpp`,
`test_gpw_xl_bomd.cpp`; `tests/cpu/scf/test_scf_evaluator.cpp`,
`test_accelerators.cpp`; `tests/cpu/gpw/test_cpu_gpw_ot.cpp` (renamed
primitives).

### CMake

- `src/cpu/linalg/orbital_tangent.cpp` → unconditional CPU sources.
- The three propagator `.cpp` move from the `ATV_BUILD_GPW` block
  (`CMakeLists.txt:645,648`) / HF block (`:577`) into the **`ATV_BUILD_HF`**
  block unconditionally — they depend only on `cpu/linalg`, and GPW already
  requires HF, so GPW keeps seeing them.

### Drive-by (named debt #6)

`src/cpu/gpw/gpw_evaluator.cpp:28` local `constexpr double kAng2Bohr = …` →
`units::kAngstromInBohr` (the shared constant every other Å→Bohr site uses).

## Files to create or modify

- NEW `include/atomiverse/cpu/linalg/orbital_tangent.h`, `src/cpu/linalg/orbital_tangent.cpp`
- NEW `include/atomiverse/cpu/electronic/{cpmd.h,xl_bomd.h,aspc_predictor.h}`,
  `src/cpu/electronic/{cpmd.cpp,xl_bomd.cpp,aspc_predictor.cpp}`
- DELETE `include/atomiverse/cpu/gpw/{cpmd.h,xl_bomd.h}`,
  `include/atomiverse/cpu/scf/aspc_predictor.h`,
  `src/cpu/gpw/{cpmd.cpp,xl_bomd.cpp}`, `src/cpu/scf/aspc_predictor.cpp`
- `include/atomiverse/cpu/gpw/internal/ot.h`, `src/cpu/gpw/internal/ot.cpp`,
  `src/cpu/gpw/gpw_solver.cpp` (renamed primitive call sites)
- `include/atomiverse/cpu/gpw/gpw_evaluator.h`, `src/cpu/gpw/gpw_evaluator.cpp`
  (types + kAng2Bohr)
- `include/atomiverse/cpu/scf/mean_field_method.h`, `src/cpu/scf/mean_field_method.cpp`
- `tests/cpu/gpw/test_gpw_cpmd.cpp`, `test_gpw_xl_bomd.cpp`, `test_cpu_gpw_ot.cpp`,
  `tests/cpu/scf/test_scf_evaluator.cpp`, `test_accelerators.cpp`
- `CMakeLists.txt`

## Tasks

- [ ] Extract `orbital_tangent.{h,cpp}` from `ot.cpp`; update `ot.h`/`ot.cpp`/
      `gpw_solver.cpp`/`test_cpu_gpw_ot.cpp` call sites; build green
- [ ] Move the three propagator pairs into `cpu/electronic/` with namespace
      `electronic::`; delete old headers/sources; migrate the 12 include sites
- [ ] Move propagator `.cpp`s to the `ATV_BUILD_HF` CMake block;
      `orbital_tangent.cpp` unconditional
- [ ] Verify an HF-only configure+build (`-DATV_BUILD_HF=ON -DATV_BUILD_GPW=OFF`)
      compiles the `electronic/` sources
- [ ] Replace `kAng2Bohr` in `gpw_evaluator.cpp` with `units::kAngstromInBohr`
- [ ] Scoped build + full ctest + cpmd/xl/sgcp regressions — bit-identical green

## Testing strategy

Pure refactor: **no new regression** (a move is not a feature). Gates are the
existing suites: `test_gpw_cpmd`, `test_gpw_xl_bomd`, `test_cpu_gpw_ot`,
`test_scf_evaluator`, `test_accelerators` plus the 15 cpmd/sgcp/xl-bomd
regressions in `regressions/dynamics/` — all must pass with unchanged numerics.
New coverage: none needed beyond compile-proof of the HF-only configuration.

## Out of scope

- Any behavior change in the propagators (advance/seed/thermostat untouched)
- New public members on any method (esmm-cpmd-02/04/05)
- CUDA mirroring (`cuda/electronic/` when CUDA grows AIMD — not now)
- Namespace alignment of other gpw/scf internals

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- The propagators must not include anything from `gpw/` or `scf/` after the
  move — `cpu/linalg` + std only.
- Grep-verify zero remaining references: `grep -rn "gpw/cpmd.h\|gpw/xl_bomd.h\|scf/aspc_predictor.h\|gpw::Cpmd\|gpw::XlBomd\|scf::ASPCPredictor" include/ src/ tests/ regressions/`.
