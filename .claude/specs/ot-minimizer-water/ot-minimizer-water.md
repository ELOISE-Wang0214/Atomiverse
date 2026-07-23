---
slug: ot-minimizer-water
status: approved
created: 2026-07-23
track: cp2k-method-loop
chapter: gpw
grilled: false
---

# ot-minimizer-water

## Summary

Productize the existing GPW **Orbital Transformation (OT)** SCF minimizer as a
first-class public switch (`use_ot` on `cpu::gpw::GpwDft` and
`cpu::gpw::GPWEvaluator`), and pin **OT total energy** for the CP2K
getting-started isolated water (PBE / DZVP-MOLOPT-GTH / GTH-PBE / MT Poisson)
against an offline **CP2K OT** twin. OT and DIIS must agree on the same
geometry within a stated energy tolerance (same self-consistent fixed point).

## Domain basis

- OT parameterization and FULL_SINGLE_INVERSE-style preconditioner already in
  `cpu/gpw/internal/ot.h` (Kühne et al., JCP 152, 194103 (2020) §VII C;
  VandeVondele & Hutter, JCP 118, 4365 (2003)).
- GPW SCF already branches on `GpwDft::use_ot` in `gpw_solver.cpp`.
- Public surface still forces OT only via diagnostic env `ATV_GPW_OT` in
  `GPWEvaluator::build` — product path must be a **member**, not an env hack.
- Literature/tutorial twin: `examples/09_gpw.cpp` / CP2K h2o getting-started
  (MT water, −17.21233714 Ha class).

## Design

### Public API (no `*Config` bags, no factories)

```cpp
cpu::gpw::GPWEvaluator eval({O, H}, /*cutoff_ry*/400.0, /*max_iter*/200,
                            /*tol*/1e-7, /*mix*/0.3, XcKind::PBE,
                            /*nonperiodic*/true);
eval.use_ot = true;   // product OT path
cpu::single_point(frame, eval);
// last_result().e_total [Hartree]; frame.meta kEnergy [eV]
```

- `GpwDft::use_ot` already exists as a member — keep it.
- Add matching **public member** `GPWEvaluator::use_ot` (default `false`) and
  wire it in `build()` into the owned `GpwDft`.
- Retain `ATV_GPW_OT` env as an optional diagnostic force-on (OR with member);
  do not require env for product code or the regression.
- No new public Config type; no `make_ot_*` factories.

### Physics / system (regression pin)

| Field | Value |
| ----- | ----- |
| Geometry | CP2K h2o tutorial H₂O [Å] in 10 Å cubic cell (same as `examples/09_gpw.cpp`) |
| Method | closed-shell Γ-point GPW, PBE, nonperiodic (MT) |
| Basis / PP | DZVP-MOLOPT-GTH + GTH-PBE-q6 (O) / GTH-PBE-q1 (H) |
| Cutoff | 400 Ry (density grid) |
| Observable | \(E_{\mathrm{tot}}\) [Hartree]; also `iters` reported (not gated hard) |
| OT vs DIIS | both paths converge; \(\|E_{\mathrm{OT}}-E_{\mathrm{DIIS}}\|\) ≤ tol_ot_diis |

### Offline CP2K oracle

- Binary: CP2K 2026.1 (Arrhenius module) or equivalent cited version.
- Twin input: OT minimizer (`&OT` / `MINIMIZER DIIS` or CP2K default OT) vs
  traditional diagonalization for the same FORCE_EVAL — record both energies
  in `REFERENCE.md`. Primary pin is **CP2K OT energy**.
- ctest **never** launches CP2K.

## Files to create or modify

- `.claude/specs/ot-minimizer-water/*` (this spec + acceptance)
- `include/atomiverse/cpu/gpw/gpw_evaluator.h` — `use_ot` member + ctor doc
- `src/cpu/gpw/gpw_evaluator.cpp` — wire `use_ot` (keep env OR)
- `tests/cpu/gpw/test_cpu_gpw_ot.cpp` (new) — unit: OT fixed-H, OT vs DIIS H2,
  evaluator member toggles OT path
- `tests/CMakeLists.txt` — register test under `ATV_BUILD_GPW`
- `regressions/gpw/ot-minimizer-water/reg_ot-minimizer-water.cpp`
- `regressions/gpw/ot-minimizer-water/REFERENCE.md`
- `tools/cp2k_oracle/inputs/ot-minimizer-water.inp` (optional twin)
- `regressions/CMakeLists.txt` — `atv_add_regression(gpw ot-minimizer-water)`
  under `ATV_BUILD_GPW`

## Tasks

- [x] Spec + acceptance committed on `wf/ot-minimizer-water`
- [x] Expose `GPWEvaluator::use_ot`; wire into `GpwDft` at build
- [x] Unit tests: internal OT fixed-H vs diag subspace; SCF OT≈DIIS on mini-H2;
      evaluator `use_ot` member default/ctor
- [x] Offline CP2K OT pin → `REFERENCE.md` (250 Ry; 400 Ry class cross-check)
- [x] E2E regression: public API OT energy vs CP2K pin; OT vs DIIS parity
- [x] Build with `ATV_BUILD_GPW=ON`; unit + regression pass
- [x] Perf gate: fair OT wall ratio ≈9.0 (OMP=4) ≤ 10

## Testing strategy

| Tier | What |
| ---- | ---- |
| `tests/cpu/gpw/test_cpu_gpw_ot.cpp` | Single-function units only; no full CP2K e2e; ≤ seconds |
| `regressions/gpw/ot-minimizer-water` | One textbook e2e vs offline CP2K pin; ≤ 60 s |

## Out of scope

- Sparse / linear-scaling OT preconditioner
- Open-shell OT / UKS OT
- k-point OT
- Driver-phase SCF rewrite of GPW
- Molecular post-HF (pyscf loop)
- MPI multi-rank scaling

## Acceptance criteria (binding)

See `ot-minimizer-water.acceptance.md`. `acceptance_met=true` only when every
criterion is verified.

---

**Librarian note:** `cpu/gpw/`; confidence high — OT engine already ships.
