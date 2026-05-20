---
title: Phase 1 — Productize the CPU SCF + gradient + opt stack
status: in-progress
created: 2026-04-14
revised: 2026-05-20
---

# Phase 1 — Productization

## Status: 35/38 landed, 3 remaining

The Phase 1 foundation is complete and in use: multi-basis support
(STO-3G / 6-31G / 6-31G* / def2-SVP via `tools/gen_basis_data.py` +
`basis_data/`), d-shell Rys ERIs end-to-end, the 12-molecule
`bench_qc_suite`, the minimizer/SCF recorder + observer components,
and the locked `default_config` / `baseline_config` presets. All
STO-3G-only entry points and back-compat shims have been removed —
every consumer uses the general-basis `build_integrals` /
`accumulate_two_electron_gradient` path (see commit
`refactor(scf): remove legacy/back-compat SCF surface`).

Three genuine follow-up items remain. They are deferred refinements,
not blockers — the production stack runs today at RHF/RKS up to
def2-SVP.

## Remaining tasks

### B2c — f-shell (l=3) Rys ERIs

`include/atomiverse/cpu/scf/rys_eri.h` has `kRysMaxL = 2` and
`kRysMaxRoots = 5`. Extending to l=3 unlocks def2-TZVP and
heavier-element bases.

- [ ] Raise `kRysMaxL` to 3; audit `kRysMaxRoots` for the worst-case
      f-quartet **and** its derivative-raised (`+1`) temporaries.
- [ ] Size the Hermite / Boys / E-coefficient tables in
      `native_integrals.cpp` from `kMaxL` so the derivative path that
      raises `ai` by `+1` (f → g-like temp) does not overflow.
- [ ] FD self-consistency test on `(f,s|s,s)`, `(f,p|s,s)`,
      `(f,d|s,s)`, `(f,f|f,f)` quartets and their derivatives.

### KS/d-shell analytic-gradient accuracy

`bench_qc_suite` reports max-force ~1 eV/Bohr at the PySCF-optimized
water geometry where PySCF reports ~1e-4 Ha/Bohr. SCF energy matches
PySCF to ~80 μHa on the same geometry, so the defect is in the KS
gradient assembly for d-shells, not the SCF. Armijo line search
currently rescues convergence but takes many backtracking iterations.

- [ ] Validate the XC piece of `gradient.cpp` against finite
      differences at 6-31G* and def2-SVP.
- [ ] Add a d-shell KS FD test alongside `test_def2svp_gradient`.

### Denser SG-1/SG-2 integration grid

`bench_qc_suite` leaves a ~5 mHa residue between ATV and PySCF on
def2-SVP because the SG-0-ish `(n_radial=50, n_angular=110)` grid in
`default_config` is coarser than PySCF's default. The
`bench_qc_suite` energy threshold is held at 1e-2 Ha until this lands.

- [ ] Implement a denser SG-1 / SG-2 radial+angular scheme.
- [ ] Tighten the `bench_qc_suite` energy threshold toward the spec's
      1 μHa target once the grid residue is closed.

## Binding rules (apply here)

- All "outputs" are `cpu::ComponentBase` subclasses bound to a
  `minimizer::Phase` or `scf::Phase`. No `Workflow` / `Runner` /
  `Report` facade classes.
- Basis-set data flows through plain value types. No runtime string
  lookups inside hot loops.
- PySCF is ground truth for integrals, SCF energies, orbital
  energies, and geometries.

## Acceptance (remaining)

- **f-shell:** FD self-consistency holds for all f-containing
  quartets and their derivatives; def2-TZVP `build_integrals`
  succeeds.
- **KS/d-shell gradient:** analytic XC gradient matches FD to
  ≤ 1e-2 eV/Å at 6-31G* and def2-SVP on water.
- **Grid:** ATV − PySCF def2-SVP single-point residue ≤ 1 μHa on the
  H2O/NH3/CH4 triplet; `bench_qc_suite` threshold tightened to match.

## Risks

- The f-shell Rys extension is the highest-risk item. Gate: FD
  self-consistency on all f-quartets must pass before def2-TZVP is
  advertised as supported.
- PySCF reference generation must stay version-pinned
  (`pyscf==2.x.y`) — re-running on a new minor release causes
  spurious test regressions.
