# Spec: `gpw-stress-water` — GPW cell stress (water)

**Slug:** `gpw-stress-water`  
**Domain:** `gpw` (Gaussian-and-plane-waves DFT)  
**Queue:** CP2K capability loop — incomplete `[~]`  
**Branch:** `wf/gpw-stress-water`  
**Status:** open  
**Resume from:** queue `[~]`; known_spec path empty on staff tree — write-first  
**Deps:** energy pin `gpw-water-pbe-mt` ✓; forces pin `gpw-forces-water` ✓; PW strain-FD stress ✓

## Problem

Plane-wave cell stress ships as strain-FD (`cpu::pw::pw_stress` → Frame
`kStress`). GPW still has **no** stress tensor on the product path. Variable-cell
relax (`gpw-variable-cell-ice`) is blocked on this foothold. Roadmap
`docs/roadmap/dft.md`: “GPW stress remains open”.

## Goal

Deliver **GPW cell stress** (v1 = **strain-FD of SCF energy**, same strategy as
PW v1) for closed-shell Γ-point GPW, with a **CP2K-pinned** water e2e regression
and a unit test. Enables `gpw-variable-cell-ice` productization.

## Out of scope

- Analytic Nielsen–Martin / Pulay stress decomposition (follow-on)
- Triclinic / off-diagonal shear strain (GPW cell is **orthorhombic** only;
  off-diagonal Voigt components are **0** in v1)
- BOMD trajectory pin (`gpw-bomd-water`)
- UKS stress, k-points, D3 stress
- Public `*Config` bags / `make_*` factories
- Live CP2K from ctest / regression binary
- CUDA GPW

## Public API (binding shape)

```cpp
// Free function — Voigt σ [Hartree/Bohr³], order xx,yy,zz,xy,xz,yz
// (matches gpw_forces_* internal units; orthorhombic: xy=xz=yz=0)
std::array<double, 6>
gpw_stress(const std::vector<GpwAtom>& atoms, const GpwDft& method,
           double h = 2e-3);

// Product path: member on GPWEvaluator (no Config bag)
cpu::gpw::GPWEvaluator eval(
    /*kinds*/ {O, H}, /*cutoff_ry*/ 400.0, /*max_iter*/ 200, /*tol*/ 1e-7,
    /*mix*/ 0.3, /*xc*/ cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
    /*compute_forces*/ false);
eval.compute_stress = true;   // public member

cpu::Evaluate driver;
driver.on(evaluate::Phase::evaluate).add_ref(eval);
driver.build();
driver.run(frame);
// frame.meta.get(kEnergy)  [eV]
// frame.meta.get(kStress)  [eV/Å³] Voigt 6
// last_stress_ha_bohr3()   [Ha/Bohr³] for pin compare
```

Rules:

- Parameters are **constructor args / members** of `GPWEvaluator` / `GpwDft`
  (iron law: no public `*Config`).
- Regression uses **public headers only** (no `internal/gpw_stress.h` in e2e).
- Units on Frame: energy eV, stress **eV/Å³**. Pin values recorded and
  asserted in **Ha/Bohr³** after conversion via `core/units.h` factors
  (`kHartreeInEv`, `kBohrInAngstrom`).
- v1 algorithm: central FD of \(E\) under diagonal strain \(\varepsilon_{aa}\) at
  fixed fractional coordinates (scale \(L_a\) and \(r_a\) by \(1\pm h\)).
  \(\sigma_{aa} = \Omega^{-1} \partial E/\partial\varepsilon_{aa}\).

## Geometry / method (binding twin)

Same water twin as `gpw-forces-water` / `gpw-water-pbe-mt`:

| Field | Value |
| ----- | ----- |
| System | isolated H₂O (CP2K getting-started coords) |
| Cell | cubic 10 Å, `PERIODIC NONE` / Atomiverse `nonperiodic=true` (MT) |
| Coordinates (Å) | O (5.0000, 5.0000, 5.0000); H (5.7586, 5.0000, 5.5043); H (4.2414, 5.0000, 5.5043) |
| Basis | DZVP-MOLOPT-GTH |
| PP | GTH-PBE-q6 (O), GTH-PBE-q1 (H) |
| XC | PBE |
| Cutoff | 400 Ry |
| Observable | cell stress Voigt [Ha/Bohr³]; optional \(E_\mathrm{tot}\) [Ha] |
| `local_pp_on_grid` | default **true** |

Twin input: `tools/cp2k_oracle/inputs/h2o_stress_mt.inp` (`STRESS_TENSOR ANALYTICAL`
or NUMERICAL; print stress).

## Tasks

1. **Spec** (this file) + acceptance + INDEX Live row; commit early on
   `wf/gpw-stress-water`.
2. **Offline CP2K pin:** run twin via `tools/cp2k_oracle/`; write
   `regressions/gpw/gpw-stress-water/REFERENCE.md` with stress table + version.
3. **Implement** `gpw_stress` (public declaration on `gpw_solver.h`; impl in
   `src/cpu/gpw/internal/gpw_stress.cpp`) + `GPWEvaluator::compute_stress` writing
   `kStress` [eV/Å³].
4. **Unit test** `tests/cpu/gpw/test_cpu_gpw_stress.cpp`: finite diagonal
   components; off-diag ≈ 0; FD self-consistency soft gate on a cheap system.
5. **Regression** `regressions/gpw/gpw-stress-water/reg_gpw-stress-water.cpp`:
   public `GPWEvaluator` + `compute_stress`; assert max |σ_atv − σ_cp2k| within tol.
6. Wire `atv_add_regression(gpw gpw-stress-water)` under `ATV_BUILD_GPW`.
7. **Perf gate:** wall \(T_\mathrm{atv}/T_\mathrm{cp2k}\) if CP2K binary used;
   else `perf_exempt` with reason.
8. Commit on `wf/gpw-stress-water` only (never `dev`, never push).

## Acceptance criteria

See `acceptance.md`.

## Offline pin strategy

1. Prefer **CP2K 2026.1** Quickstep, PBE, DZVP-MOLOPT-GTH, GTH-PBE, MT Poisson,
   10 Å cell, coords above, CUTOFF 400, `STRESS_TENSOR` print.
2. Record Voigt stress (xx,yy,zz,xy,xz,yz) in Ha/Bohr³ into REFERENCE.
3. If CP2K binary unavailable: document unavailability; pin ATV strain-FD vs
   independent second-h FD of the same energy (self-consistency) + energy
   anchor from `gpw-water-pbe-mt`; set `perf_exempt=true`. Do **not** invent
   absolute stress digits from memory.
4. Fair wall: same geometry; single-rank; write ratio to REFERENCE.

## Files (expected touch list)

- `.claude/specs/gpw-stress-water/{spec.md,acceptance.md}`
- `.claude/specs/INDEX.md` (Live row)
- `tools/cp2k_oracle/inputs/h2o_stress_mt.inp`
- `include/atomiverse/cpu/gpw/gpw_solver.h` (declare `gpw_stress`)
- `include/atomiverse/cpu/gpw/gpw_evaluator.h` (+ `compute_stress`, last stress)
- `include/atomiverse/cpu/gpw/internal/gpw_stress.h`
- `src/cpu/gpw/internal/gpw_stress.cpp`
- `src/cpu/gpw/gpw_evaluator.cpp`
- `CMakeLists.txt` (add source)
- `tests/cpu/gpw/test_cpu_gpw_stress.cpp` + `tests/CMakeLists.txt`
- `regressions/gpw/gpw-stress-water/{reg_gpw-stress-water.cpp,REFERENCE.md}`
- `regressions/CMakeLists.txt`

## Done when

All acceptance rows hold; `acceptance_met=true`; branch committed; report JSON
with `status=green` (or `blocked` only for speed after one optimize pass /
unfixable host dep).
