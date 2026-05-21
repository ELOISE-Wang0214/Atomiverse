---
title: B3LYP / PBE0 analytic nuclear gradient and hybrid geometry optimization
status: approved
created: 2026-05-21
---

# B3LYP / PBE0 analytic nuclear gradient and hybrid geometry optimization

## Summary

This is the third and final link of the hybrid-DFT chain. It implements the analytic nuclear gradient for the B3LYP and PBE0 global-hybrid functionals so that hybrid DFT supports forces, not only single-point energies. The work scales the exchange (K) channel of the existing two-electron gradient by the functional's exact-exchange fraction `alpha_hf` (0.20 for B3LYP, 0.25 for PBE0) and routes the B88/LYP semi-local exchange-correlation contribution through the existing GGA XC-gradient machinery. Once landed, a user can hand a `b3lyp_config` or `pbe0_config` (produced by link 2) to the existing `cpu::minimizer::Minimize` + `SCFEvaluator` + `LBFGS` stack and run a hybrid geometry optimization with no further code — that composition is verified by a checked-in test. After this link the full hybrid-DFT capability (energy + forces + geometry optimization) is complete.

## Domain basis

The hybrid nuclear gradient differs from the pure-DFT gradient only in the exact-exchange admixture. Verified this session against PySCF `grad/rks.py` (`Gradients.get_veff`):

```
vj, vk = ks_grad.get_jk(mol, dm)
vk *= hyb               # hyb == alpha_hf
vxc += vj - vk * .5
```

- **K-channel multiplier.** The two-electron nuclear gradient is `dE_2e/dR = dJ/dR - alpha_hf * dK/dR`, where the closed-shell exchange spin factor (the `* .5` above, folded into Atomiverse's existing `factor_K` as the `-0.25` restricted prefactor) is unchanged. The exact-exchange fraction `alpha_hf` is the *only* new multiplier: it is `0.0` for pure DFT (K skipped), `1.0` for pure HF, and exactly `alpha_hf` for a global hybrid — `0.20` for B3LYP, `0.25` for PBE0. `omega = 0` for both functionals (global, not range-separated), so there is no range-separated K-derivative.
- **Semi-local XC gradient is untouched by `alpha_hf`.** The derivative of the semi-local `E_xc^sl` (`accumulate_xc_gradient`) uses the semi-local functional *exactly* as the pure-DFT path. For B3LYP the XC gradient must evaluate the B88 exchange + LYP correlation GGA contributions — the `vrho`/`vsigma` chain rule through `grad-rho`; for PBE0 the existing PBE XC-gradient path is reused with the functional's 0.75 semi-local exchange coefficient. The GGA XC-gradient assembly is functional-agnostic given `vsigma`, so it must already generalize once link 1's B88/LYP `vrho`+`vsigma` are routed to `evaluate_xc_batch`.
- **Grid-response term unchanged.** The Becke fuzzy-cell weight-derivative and grid-coordinate terms (merged spec `ks-xc-grid-response-gradient`, present in `accumulate_xc_gradient`) apply unchanged. They make `sum_A F_A = 0` exact, so the translational-invariance projection (`SCFConfig::xc_grad_project_translation`) stays off.

References:
- Becke, "Density-functional thermochemistry. III. The role of exact exchange", J. Chem. Phys. 98, 5648 (1993). DOI: 10.1063/1.464913 — B3LYP.
- Adamo & Barone, "Toward reliable density functional methods without adjustable parameters: The PBE0 model", J. Chem. Phys. 110, 6158 (1999). DOI: 10.1063/1.478522 — PBE0.
- Johnson, Gill & Pople, "The performance of a family of density functional methods", J. Chem. Phys. 98, 5612 (1993). DOI: 10.1063/1.464906 — DFT nuclear-gradient AO-derivative convention.
- PySCF `pyscf/grad/rks.py` — `get_veff` hybrid `vk *= hyb` reference implementation.

Performance note (non-binding): the hybrid K-derivative re-walks every shell quartet, roughly doubling the two-electron-gradient cost relative to pure DFT for `alpha_hf > 0`. Optimizing the hybrid K-derivative is explicitly out of scope.

## Design

This link is a *consumer* of link 1 and link 2. The dependency surface it assumes is in place:

- Link 2 extends `XCType` (in `xc_functional.h`) with `B3LYP` and `PBE0`, and provides a hybrid-fraction accessor — `hybrid_exchange_fraction(XCType)` returning `0.0` for `LSDA`/`PBE`, `0.20` for `B3LYP`, `0.25` for `PBE0` — plus an `is_gga(XCType)` predicate that is `true` for `PBE`, `B3LYP`, `PBE0`. Link 2 also provides `b3lyp_config(...)` / `pbe0_config(...)` builders.
- Link 1 makes `evaluate_xc_batch` produce correct `vrho`+`vsigma` for the B88/LYP semi-local pieces underlying B3LYP.

Entities touched in this link:

- `accumulate_two_electron_gradient` (`integral_derivatives.cpp`). The K-channel currently fires only when `!is_ks` (HF, multiplier `1.0`); for KS it is skipped. This link replaces the `!is_ks` gate with a numeric `alpha_hf = hybrid_exchange_fraction(config.xc_type)` taken for KS and `1.0` for HF, then multiplies every `factor_K` term by it. The existing `-0.25` (restricted) / `-0.5` (unrestricted) closed-shell exchange prefactors and the `0.5·[P(i,k)P(j,l)+P(j,k)P(i,l)]` symmetrization are kept verbatim — `alpha_hf` is a pure scalar multiplier on top. Both the libcint branch and the native Rys block-deriv branch carry the identical change. When `alpha_hf == 0.0` the K loop must early-skip exactly as today (no behavior change for `LSDA`/`PBE`).
- `accumulate_xc_gradient` (`gradient.cpp`). The `is_gga` flag is currently `config.xc_type == XCType::PBE`; it must become `is_gga(config.xc_type)` so B3LYP and PBE0 route through the GGA `vsigma` path. No new GGA-gradient math: the `gga_contraction` lambda and Term1/Term2 chain rule are functional-agnostic given the `vsigma_*` arrays returned by `evaluate_xc_batch`. Verify B88/LYP `vsigma` flows through unchanged; if any `XCType::PBE`-specific branch remains beyond the `is_gga` flag, generalize it.
- Geometry optimization is composition-only — no new production code. `cpu::minimizer::Minimize` + `SCFEvaluator{b3lyp_config}` + `LBFGS` already compose; this link only adds a test exercising that path.

Ownership/lifecycle: no new types, no new headers, no new allocations in any `launch()` path. `scf_gradient` and `SCFEvaluator` are unchanged at the signature level — the hybrid behavior is entirely data-driven through `SCFConfig::xc_type`.

## Files to create or modify

- `src/cpu/scf/integral_derivatives.cpp` — scale the K-channel `factor_K` by `alpha_hf` in both the libcint and native Rys branches of `accumulate_two_electron_gradient`; replace the `!is_ks` K-gate with `alpha_hf > 0`.
- `src/cpu/scf/gradient.cpp` — replace `is_gga = (config.xc_type == XCType::PBE)` with `is_gga(config.xc_type)` in `accumulate_xc_gradient`; confirm no other `XCType::PBE`-specific branch.
- `src/cpu/scf/xc_functional.cpp` — only if link 1 has not already wired the B88/LYP GGA-gradient `vrho`+`vsigma` potential pieces into `evaluate_xc_batch`; otherwise no change.
- `tests/cpu/scf/test_hybrid_gradient.cpp` (new) — analytic-vs-PySCF, analytic-vs-finite-difference, and translational-invariance tests for B3LYP and PBE0 at 6-31G*, plus a B3LYP/6-31G* geometry-optimization composition test, all against checked-in PySCF reference data.
- `tests/CMakeLists.txt` — register `test_hybrid_gradient` in the CPU SCF test group (guarded by `ATV_BUILD_RHF`).

## Tasks

- [ ] Write failing tests for the hybrid two-electron K-channel scaling (tests/cpu/scf/test_hybrid_gradient.cpp) — B3LYP and PBE0 analytic force vs checked-in PySCF `grad.RKS` reference on H2O/CH4/NH3 at 6-31G*
- [ ] Write failing analytic-vs-finite-difference and translational-invariance tests for B3LYP/PBE0 (tests/cpu/scf/test_hybrid_gradient.cpp), moving-grid central difference h=1e-3 Bohr, SCF tol 1e-9
- [ ] Write failing B3LYP/6-31G* geometry-optimization composition test (tests/cpu/scf/test_hybrid_gradient.cpp) using cpu::minimizer::Minimize + SCFEvaluator + LBFGS
- [ ] Scale the K-channel by alpha_hf in accumulate_two_electron_gradient (src/cpu/scf/integral_derivatives.cpp) — libcint and native Rys branches, replacing the !is_ks K-gate with hybrid_exchange_fraction(config.xc_type)
- [ ] Generalize the GGA route in accumulate_xc_gradient (src/cpu/scf/gradient.cpp) to use is_gga(config.xc_type) so B3LYP/PBE0 reach the vsigma chain rule
- [ ] Verify B88/LYP vrho+vsigma flow through evaluate_xc_batch for the B3LYP XC gradient (src/cpu/scf/xc_functional.cpp); wire only if link 1 left a gap
- [ ] Add Doxygen per doc.style with units to every changed function and the new test fixtures (src/cpu/scf/integral_derivatives.cpp, src/cpu/scf/gradient.cpp)
- [ ] Register test_hybrid_gradient in tests/CMakeLists.txt under ATV_BUILD_RHF
- [ ] Run full check + test suite

## Testing strategy

- **Happy path — analytic vs PySCF reference.** B3LYP analytic force matches `pyscf.grad.RKS(mf).kernel()` with `xc='b3lyp'` to `max|delta| <= 1e-5` Ha/Bohr on H2O, CH4, NH3 at 6-31G*; PBE0 likewise with `xc='pbe0'`. Reference forces are genuine PySCF runs, checked into the test file as inline arrays (same convention as `test_def2svp_gradient.cpp`), with the PySCF version and exact level of theory recorded in a comment.
- **Self-consistency — analytic vs finite difference.** B3LYP and PBE0 analytic force matches a moving-grid central-difference oracle (`scf_finite_difference_gradient`, grid regenerated per displaced geometry, `h = 1e-3` Bohr, SCF tol `1e-9`) to `<= 1e-5` Ha/Bohr. This guards against any PySCF convention mismatch and is the independent oracle.
- **Edge cases.** Translational invariance: `sum_A F_A = 0` per Cartesian axis to `<= 1e-8` Ha/Bohr for both functionals (exact because the merged Becke grid-response term is in place). Regression: `alpha_hf == 0` for `LSDA`/`PBE` must leave the existing pure-DFT and HF gradients bit-unchanged — the full ctest suite is the guard.
- **Domain validation.** Hybrid geometry optimization: a B3LYP/6-31G* optimization of water (or another small molecule) converges, and the optimized geometry RMS is within `1e-3` Å and the final energy within `1e-6` Ha of the PySCF optimized reference. This is the end-to-end proof that energy + force + optimizer compose.
- **Suite integrity.** Full `ctest` stays green — no regression to HF, LSDA, or PBE gradients.

## Out of scope

- Range-separated hybrids (`omega != 0`, e.g. CAM-B3LYP, wB97X) — global hybrids only.
- Analytic Hessian / second derivatives.
- Performance optimization of the hybrid K-derivative — the K-channel re-walk is accepted as-is (non-binding performance note above).
- CUDA backend hybrid gradient — CPU SCF only; CUDA SCF migration is tracked separately in `cuda-driver-taxonomy.md`.
