# Spec: `gpw-forces-water` — GPW analytic forces (water)

**Slug:** `gpw-forces-water`  
**Domain:** `gpw` (Gaussian-and-plane-waves DFT)  
**Queue:** CP2K capability loop — incomplete `[~]`  
**Branch:** `wf/gpw-forces-water`  
**Status:** open  
**Resume from:** queue `[~]`; no prior open spec (write-first)

## Problem

Hellmann–Feynman + Pulay analytic GPW forces already ship in-engine
(`gpw_forces_analytic`, `GPWEvaluator::compute_forces`). There is **no**
out-of-library e2e regression under `regressions/gpw/` that pins water
cart forces against an offline CP2K `FORCE_EVAL` twin. The deleted unit
suite (`test_cpu_gpw_forces.cpp`) only checked analytic-vs-FD self-
consistency, never CP2K parity. `local_pp_on_grid` must stay energy/force
consistent (historical ZnO 0.53 Ha/Bohr bug).

## Goal

Deliver a **CP2K-pinned** water analytic-force regression on the public
`cpu::gpw::GPWEvaluator` path (Frame I/O, eV on surface, pin in Ha/Bohr),
plus a fast unit test that analytic forces match the FD oracle on a tiny
system. Enables `gpw-stress-water` / `gpw-bomd-water` downstream.

## Out of scope

- Cell stress (`gpw-stress-water`)
- BOMD trajectory pin (`gpw-bomd-water`)
- Open-shell GPW UKS (`gpw-uks-o2-gamma`)
- k-point GPW
- Re-architecting force kernels (only fix if pin fails)
- Public `*Config` bags / `make_*` factories
- Live CP2K from ctest / regression binary
- CUDA GPW

## Public API (minimal — already mostly shipped)

```cpp
cpu::gpw::Kind O{molrs::Element::O,
                 cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
                 cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
cpu::gpw::Kind H{molrs::Element::H,
                 cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
                 cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

cpu::gpw::GPWEvaluator eval(
    /*kinds*/ {O, H},
    /*cutoff_ry*/ 400.0,
    /*max_iter*/ 200,
    /*tol*/ 1e-7,
    /*mix*/ 0.3,
    /*xc*/ cpu::gpw::XcKind::PBE,
    /*nonperiodic*/ true,
    /*compute_forces*/ true);   // analytic forces → frame.span(kAtom_f*) [eV/Å]

cpu::Evaluate e;
e.on(evaluate::Phase::evaluate).add_ref(eval);
e.build();
e.run(frame);
// energy: frame.meta.get(kEnergy) [eV]
// forces: kAtom_fx/y/z [eV/Å]; convert to Ha/Bohr for pin compare
```

Rules:

- Parameters are **constructor args / members** of `GPWEvaluator` / `GpwDft`
  (iron law: no public `*Config`).
- Regression uses **public headers only** (no `internal/gpw_forces.h`).
- Units on Frame: eV, eV/Å. Pin values recorded and asserted in **Ha/Bohr**
  (and optional total energy Ha) after conversion via `core/units.h`.

## Geometry / method (binding twin)

| Field | Value |
| ----- | ----- |
| System | isolated H₂O (CP2K getting-started coords) |
| Cell | cubic 10 Å, `PERIODIC NONE` / Atomiverse `nonperiodic=true` (MT) |
| Coordinates (Å) | O (5.0000, 5.0000, 5.0000); H (5.7586, 5.0000, 5.5043); H (4.2414, 5.0000, 5.5043) |
| Basis | DZVP-MOLOPT-GTH |
| PP | GTH-PBE-q6 (O), GTH-PBE-q1 (H) |
| XC | PBE |
| Cutoff | 400 Ry (product path; regression may lower if ≤60 s still meets tol) |
| Observable | cart forces [Ha/Bohr]; optional \(E_\mathrm{tot}\) [Ha] |
| `local_pp_on_grid` | default **true** (product path); document if pin uses alternate |

## Tasks

1. **Spec** (this file) + INDEX Live row; commit early on `wf/gpw-forces-water`.
2. **Offline CP2K pin:** run twin `ENERGY_FORCE` via
   `tools/cp2k_oracle/inputs/` (record version, input hash path, forces).
   Write `regressions/gpw/gpw-forces-water/REFERENCE.md`.
3. **Unit test** `tests/cpu/gpw/test_cpu_gpw_forces.cpp` (or slim equivalent):
   analytic vs FD on a tiny closed-shell system (H₂ or distorted H₂O),
   gate max\|F_a−F_fd\| ≤ 5e-5 Ha/Bohr (or documented tighter). Uses
   public/`GpwDft` + FD API; **no** full CP2K e2e in ctest.
4. **Regression** `regressions/gpw/gpw-forces-water/reg_gpw-forces-water.cpp`:
   public `GPWEvaluator{…, compute_forces=true}` + Evaluate (or single_point
   pattern with forces); assert max component error vs CP2K pin within tol.
   Also: finite converged energy; translational invariance
   \(\|\sum_A \mathbf{F}_A\|\) small (sanity).
5. Wire `atv_add_regression(gpw gpw-forces-water)` under `ATV_BUILD_GPW`.
6. **Perf gate:** wall \(T_\mathrm{atv}/T_\mathrm{cp2k}\) if CP2K binary used
   for pin; else `perf_exempt` with reason in REFERENCE.
7. Commit on `wf/gpw-forces-water` only (never `dev`, never push).

## Acceptance criteria

| # | Criterion | How verified |
| - | --------- | ------------ |
| A1 | Spec + Live INDEX row exist for `gpw-forces-water` | Files present |
| A2 | No public `*Config` / `make_*` added for forces | Code review / factory scrub |
| A3 | Unit test: analytic vs FD max error ≤ **5×10⁻⁵ Ha/Bohr** on a tiny system | `ctest -R gpw_forces` (or named test) pass |
| A4 | Regression `reg_gpw-forces-water` builds with `ATV_BUILD_GPW` + `ATV_BUILD_REGRESSIONS` | CMake + link |
| A5 | SCF converges (or finite energy with `converged` documented); \(E\) finite | Regression exit 0 |
| A6 | max_a,k \|F_atv − F_cp2k\| ≤ **1×10⁻³ Ha/Bohr** (default product tol; may tighten if pin allows) | Regression vs REFERENCE pin |
| A7 | \(\|\sum_A \mathbf{F}_A\|\) ≤ **5×10⁻⁴ Ha/Bohr** (translational residual; CP2K twin residual ~1.5e-4) | Regression |
| A8 | `REFERENCE.md` cites CP2K version + input settings + force table **or** explicit unavailability + FD/self-consistency pin with literature method cite (VandeVondele et al. CPC 167, 103 (2005)) | File present |
| A9 | Dual gate: correctness (A5–A7) **and** (`within_10x` **or** `perf_exempt=true` with reason) | REFERENCE + report |
| A10 | No live third-party process from ctest / regression binary | Design / grep |
| A11 | `local_pp_on_grid` default path is what the regression exercises (document if false) | REFERENCE |

## Offline pin strategy

1. Prefer **CP2K 2024+ / 2026.1** Quickstep `RUN_TYPE ENERGY_FORCE`, PBE,
   DZVP-MOLOPT-GTH, GTH-PBE, MT Poisson, 10 Å cell, coords above, CUTOFF 400.
2. Record Atomic Forces from FORCE_EVAL print in Ha/Bohr into REFERENCE.
3. If CP2K binary unavailable on host: pin via **analytic vs FD** of the same
   `GpwDft` energy (self-consistency) + energy anchor `-17.21233714` Ha from
   CP2K getting-started water (already cited in `examples/09_gpw.cpp` /
   ch09); set `perf_exempt=true`. Do **not** invent absolute force digits
   from memory.
4. Fair wall: same geometry; single-rank; write ratio to REFERENCE.

## Files (expected touch list)

- `.claude/specs/gpw-forces-water/spec.md` (this)
- `.claude/specs/INDEX.md` (Live row)
- `tools/cp2k_oracle/inputs/h2o_forces_mt.inp` (optional twin input)
- `tests/cpu/gpw/test_cpu_gpw_forces.cpp` + `tests/CMakeLists.txt`
- `regressions/gpw/gpw-forces-water/{reg_gpw-forces-water.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt`
- Force kernel sources **only if** pin fails (`src/cpu/gpw/internal/gpw_forces.cpp`, …)

## Done when

All acceptance rows A1–A11 hold; `acceptance_met=true`; branch committed;
report JSON with `status=green` (or `blocked` only for speed after one
optimize pass / unfixable host dep).
