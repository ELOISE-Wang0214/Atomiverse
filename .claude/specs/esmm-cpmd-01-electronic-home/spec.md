---
slug: esmm-cpmd-01-electronic-home
status: code-complete
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
- **Architect deltas (2026-07-26 approve-with-notes):**
  - 7 `regressions/dynamics/*/REFERENCE.md` carry old qualified names
    (cpmd-water, cpmd-electronic-thermostat, cpmd-kpoints, cpmd-rattle-water,
    cpmd-uks, second-gen-cpmd-kpoints, second-gen-cpmd-water) — AC-002's grep
    scans `regressions/`, so they are in scope.
  - 8 docs pages reference old names/paths: `docs/tutorial/methods/{cpmd-water,
    cpmd-electronic-thermostat,cpmd-rattle-water,second-gen-cpmd-kpoints,
    xl-bomd-linearized-forces,xl-bomd-uks-kpoints,xl-bomd-water}.md`,
    `docs/roadmap/molecular-dynamics.md`.
  - Tests-mirror-source: pure propagator unit tests split out to NEW
    `tests/cpu/electronic/{test_cpmd,test_xl_bomd,test_aspc_predictor}.cpp`
    (+ `tests/CMakeLists.txt` registrations); GPW/SCF integration guards stay.
  - CLAUDE.md + AGENTS.md backend-directory tables gain the `electronic/` row.

## Tasks

- [x] Extract `orbital_tangent.{h,cpp}` from `ot.cpp`; update `ot.h`/`ot.cpp`/
      `gpw_solver.cpp`/`test_cpu_gpw_ot.cpp` call sites; build green
      (also relocated `eigh_hermitian` → `cpu/linalg/eigensolve.h` — the GPW-gated
      kpoints home would have re-inverted layering; transitional `using` in
      kpoints.h to be removed in tasks 2/5)
- [x] Move the three propagator pairs into `cpu/electronic/` with namespace
      `electronic::`; delete old headers/sources; migrate the 12 include sites
      (also deduped cpmd.cpp's private c_mul/c_dagger_mul → linalg::)
- [x] Move propagator `.cpp`s to the `ATV_BUILD_HF` CMake block;
      `orbital_tangent.cpp` unconditional
- [x] Update the 7 REFERENCE.md + 8 docs pages to `electronic::` names/paths;
      add the `electronic/` row to CLAUDE.md + AGENTS.md dir tables
      (+ grep-miss fixes in the CLAUDE/AGENTS scf-utilities line; roadmap
      multi-k CPMD status corrected to shipped-RKS)
- [x] Split pure propagator unit tests to `tests/cpu/electronic/` (mirror rule);
      GPW/SCF integration guards stay in place (41 moved / 59 stay; counts
      conserved; full ctest 375/375)
- [x] Verify an HF-only configure+build (`-DATV_BUILD_HF=ON -DATV_BUILD_GPW=OFF`)
      compiles the `electronic/` sources (scratch build green; electronic+
      orbital_tangent objects present; gpw objects = dispersion family only,
      the pre-existing named debt)
- [x] Replace `kAng2Bohr` in `gpw_evaluator.cpp` with `units::kAngstromInBohr`
      (NOTE: real numeric change — old local constant was truncated, rel Δ
      7.2e-8; "bit-identical" amended to "pins pass with no tolerance edits")
- [x] Scoped build + full ctest + cpmd/xl/sgcp regressions — pins pass with
      no tolerance edits (full ctest 375/375; 15/15 spec regressions PASS +
      bonus reg_gpw-bomd-water PASS; kAng2Bohr precision fix absorbed by
      envelope pins with zero tolerance edits)

- [x] Hygiene: /mol:simplify ran clean (dead transitional using deleted; 11
      clang-format line fixes on diff-touched lines; zheevd_ linkage warning
      fixed ABI-identically; full ctest green after batch)

> Architect (2026-07-26) also NAMED pre-existing debt, routed separately (not
> this spec): 🚨 scf ⇄ gpw cycle via `src/cpu/scf/dftd3.cpp:13-14` including
> `gpw/internal/dftd3{,_params}.h` (shared Grimme engine parked in a foreign
> domain's internal/ — same anti-pattern this spec fixes, running backwards);
> HF-only build check (2026-07-26) shows the debt is the WHOLE dispersion
> family: dftd3/dftd3_params/hirshfeld_volumes/mbd_dispersion/rvv10/
> ts_dispersion/vdw_df/xdm_dispersion .cpp all live in gpw/internal yet
> compile under the ATV_BUILD_HF block — a future `cpu/dispersion/` (or
> similar) extraction is the clean fix, mirroring this spec's move;
> 🟢 ASPCPredictor's raw `xt::xtensor<double,2>` vs `linalg::RealMatrix`
> layout inconsistency (no silent change in a zero-behavior move).

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
