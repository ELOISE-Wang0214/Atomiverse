# xl-bomd-water — acceptance

Each criterion is observable + tolerance. Verified by the scoped gate
(unit tests `test_cpu_gpw_xl_bomd`, existing `test_cpu_gpw_solver` /
`test_cpu_gpw_bomd` untouched-behavior guards, and the regression
`run_reg_xl-bomd-water`).

- [ ] **A1 — coefficient pin.** `XlBomd(5)` exposes exactly the Niklasson
  JCP 130, 214109 (2009) Table-I row: `κ = 1.82`, `α = 0.018`,
  `c = {−6, 14, −8, −3, 4, −1}` (bitwise `EXPECT_DOUBLE_EQ`). For every
  `K ∈ [3, 9]`: `Σ_k c_k == 0` and `Σ_k k·c_k == 0` (exact integer sums).
- [ ] **A2 — order validation.** `XlBomd(2)` and `XlBomd(10)` throw
  `std::invalid_argument`; `advance` with an empty or shape-mismatched
  matrix throws.
- [ ] **A3 — seed + fixed point.** First `advance(D)` makes
  `ready() == true` and `aux() == D` elementwise; ≥ 20 further advances
  with the same constant `D` keep `max|aux − D| ≤ 1e-12` (the recursion's
  fixed point — Σc_k = 0 makes the dissipation vanish on a constant
  history).
- [ ] **A4 — bounded response.** 200 advances driven by
  `D_n = D0 + 0.01·sin(0.3 n)·B` (fixed symmetric `B`, `‖B‖∞ = 1`) keep
  `max|aux_n − D_n| ≤ 0.1` (no secular growth / instability at K = 5).
- [ ] **A5 — guess_P plumbing.** On the tiny in-test H₂ `GpwDft` system
  (8 Bohr box, 120 Ry), `run(atoms, nullptr, &P_conv)` with the cold run's
  converged density converges with `|e_total − e_cold| ≤ 1e-10` Ha and
  `iters ≤ iters_cold`. Non-null `guess_P` with UKS or `use_ot` configured
  throws `std::invalid_argument`.
- [ ] **A6 — evaluator guards.** `GPWEvaluator` with `xl_bomd = true` and
  (`n_alpha + n_beta > 0` or `use_ot = true`) throws at `build()`.
- [ ] **A7 — e2e XL-BOMD NVE envelope (regression).** H₂O
  (CP2K getting-started geometry, 10 Å cubic MT cell, PBE,
  DZVP-MOLOPT-GTH, GTH-PBE, 200 Ry, SCF tol 1e-5, dt 0.5 fs), `xl_bomd =
  true`, `xl_bomd_k = 5`: 1 warm-up + 4 production steps give relative
  total-energy drift `max_t |E(t) − E0| / max(1, |E0|) ≤ 0.05`; PE
  < −100 eV after warm-up; KE finite, ≥ 0, non-trivial; all coordinates
  finite; both r(O–H) < 1.8 Å at the end; the final SCF is converged
  (`last_result().converged == true`, iterations < max_iter).
- [ ] **A8 — regression budget.** `reg_xl-bomd-water` wall ≤ 60 s
  (login-class, OMP as configured).
- [ ] **A9 — no-regression guard.** Existing scoped tests
  (`test_cpu_gpw_solver`, `test_cpu_gpw_bomd`, `test_cpu_gpw_forces`) still
  pass — default (non-XL) behavior unchanged (new parameters default to
  the old code path).
- [ ] **A10 — perf gate.** Fair CP2K AIMD twin wall if a CP2K binary is
  available on the worker host; else `perf_exempt = true` with the reason
  recorded in `REFERENCE.md` (2026-07-25 host inventory: no cp2k binary /
  module — expected exempt).
