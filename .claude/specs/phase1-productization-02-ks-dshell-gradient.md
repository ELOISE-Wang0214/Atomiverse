---
title: KS d-shell analytic-gradient accuracy fix
status: approved
created: 2026-05-20
---

# KS d-shell analytic-gradient accuracy fix

## Summary
The Kohn-Sham (DFT) analytic nuclear gradient is wrong by roughly four orders of
magnitude for basis sets that include d-shell functions. At a PySCF-optimized
water geometry where PySCF reports a near-zero max force (~1e-4 Ha/Bohr), the
Atomiverse KS gradient reports ~1 eV/Bohr (~0.037 Ha/Bohr), even though the SCF
*energy* on the identical geometry matches PySCF to ~80 microHartree. The defect
is isolated to the exchange-correlation term of the gradient assembly in
`src/cpu/scf/gradient.cpp` for d-shell (l=2) basis functions; s/p gradients and
the HF (exact-exchange) gradient are already correct. This spec diagnoses the
XC-gradient assembly for d-shells, fixes it, and adds a finite-difference
regression test so the geometry PySCF calls optimized is recognized as
near-stationary. The user-visible outcome is that geometry optimization at
PBE/6-31G* and PBE/def2-SVP converges with correct forces instead of relying on
Armijo backtracking to mask the defect.

## Domain basis

The total KS nuclear gradient decomposes into five terms:

```
dE/dR_A = dE_nn/dR_A                              (nuclear repulsion)
        + sum_munu P_munu dH_munu/dR_A            (one-electron core: dT/dR + dV_ne/dR)
        + 1/2 sum P_munu P_lamsig d(munu|lamsig)/dR_A   (two-electron Coulomb dJ/dR)
        - sum_munu W_munu dS_munu/dR_A            (Pulay / overlap-derivative term)
        + dE_xc/dR_A                              (exchange-correlation)
```

with `W` the energy-weighted density matrix, `W_munu = sum_i eps_i c_mui c_nui`.

The XC gradient `dE_xc/dR_A` has three contributions:

- (a) basis-function-derivative contribution: `sum_g w_g (df_xc/dR_A)|explicit`
- (b) grid-weight-derivative contribution: `sum_g (dw_g/dR_A) f_xc(r_g)`
- (c) grid-point-coordinate chain-rule contribution: `sum_g w_g (df_xc/dr_g)(dr_g/dR_A)`

For a GGA (e.g. PBE), term (a) expands via the chain rule through both rho and
grad-rho: the XC-potential pieces `df_xc/drho` ("vrho") and
`df_xc/d|grad rho|^2` ("vsigma") contracted with derivatives of rho and
grad-rho, which carry derivatives (and, for grad-rho, second derivatives) of the
Cartesian Gaussian basis functions. References: Johnson, Gill, Pople,
*J. Chem. Phys.* 98, 5612 (1993), DOI:10.1063/1.464906; Hirata et al.,
arXiv:2204.13769 Eqs. 6 & 9.

Most likely bug class for d-shells: the Cartesian-polynomial-prefactor
derivative. A Cartesian d-function is `chi = x^i y^j z^k exp(-alpha r^2)`,
`i+j+k=2`. With `x = (r - R_A)_x`, the nuclear derivative is

```
d chi/dX_A = -i x^(i-1) y^j z^k g  +  2 alpha x^(i+1) y^j z^k g
```

The `-i x^(i-1)` polynomial-lowering term vanishes for s (i=0), is a trivial
constant for p, but for d-functions it is non-trivial **and component-dependent
within one shell** (i=2 for xx, i=1 for xy, i=0 for yy). A hand-written gradient
correct for s/p will silently drop or mis-index this term for d. The
energy-right / gradient-wrong / d-isolated signature is exactly consistent with
this. Second-most-likely: the GGA chain rule through grad-rho needs **second**
derivatives of d basis functions (the gradient of a d-function contains f-like
and p-like Cartesian terms) — the d-shell second-derivative table.

An energy match with a grossly wrong gradient points to a missing/wrong
basis-function-derivative term (a), not a grid-weight (b) omission — grid-weight
omission produces only sub-mHa/Bohr residual forces on a converged grid.

FD validation protocol: central differences `F_FD = -[E(R+h) - E(R-h)]/(2h)`,
step `h = 1e-3 Bohr` (h_opt ~ (eps_E / 2C)^(1/3) ~ 1e-3 for SCF energies
converged to ~1e-9 Ha). Use the **same fixed (non-rotating) grid** for the + and
- displacements or grid-weight noise dominates. Tighten SCF convergence to
<= 1e-9 Ha. Expected analytic-vs-FD agreement for a correct GGA gradient is
<= 1e-5 Ha/Bohr per component (1e-6 with a dense grid). Worse than ~1e-4 Ha/Bohr
indicates a genuine bug.

## Design

The defect site is `accumulate_xc_gradient()` in `src/cpu/scf/gradient.cpp`,
which assembles `dE_xc/dR_A` from per-grid-point quantities. It evaluates basis
functions through the `eval_basis_*` helpers in
`include/atomiverse/cpu/scf/basis_eval.h` (value, gradient, Hessian), then
contracts them against `vrho`/`vsigma` and the density matrix. The LDA path
forms `drho/dR_A` from `dPhi`; the GGA path additionally forms `dgrad_rho/dR_A`,
whose "Term2" uses the basis-function Hessian as the nuclear derivative of
`dPhi`.

The nuclear derivative `d/dR_A` of any quantity built from a basis function on
atom A equals the negative of its spatial derivative `-d/dr`. The s/p paths get
this right because their polynomial prefactors are constant or linear; the d
path is where a component-dependent polynomial-lowering term enters. Diagnosis
proceeds term by term, isolating which of the XC sub-contributions (a)/(b)/(c)
is wrong by comparing the d-shell case against FD with a fixed grid.

Two candidate fix surfaces:

1. **Assembly bug in `gradient.cpp`** — `accumulate_xc_gradient` mis-indexes or
   drops the d-shell polynomial-lowering / Hessian term. This is the primary
   suspect; the fix is local to `gradient.cpp`.
2. **Helper bug in `basis_eval.h`** — `poly_grad` / `poly_hess` /
   `eval_basis_all` produce a wrong d-shell derivative table. Less likely
   because `eval_basis` is also exercised in V_xc construction (energy is
   correct), but the gradient path uses `eval_basis_all`/`eval_basis_val_grad`
   which the energy path may not. The existing `test_basis_eval_dshell.cpp`
   already FD-checks `poly_{grad,hess}`; if it passes, the bug is in (1).

Entities touched: `accumulate_xc_gradient` (gradient assembly), possibly the
`detail::poly_grad`/`poly_hess`/`eval_basis_all` helpers. No new public symbols;
no lifecycle or ownership changes. `scf_finite_difference_gradient` is reused
as-is as the FD oracle. The diagnosis is recorded in the corrected code's
Doxygen comment so the failure mode is documented for future basis extensions.

## Files to create or modify

- `src/cpu/scf/gradient.cpp` — fix the d-shell XC-gradient assembly in
  `accumulate_xc_gradient`; update the explanatory comment block with the
  diagnosed term.
- `include/atomiverse/cpu/scf/basis_eval.h` — only if diagnosis shows the
  d-shell `poly_grad`/`poly_hess`/`eval_basis_all` table is wrong; otherwise
  unchanged.
- `tests/cpu/scf/test_ks_dshell_gradient.cpp` (new) — d-shell KS finite-
  difference regression test (PBE/6-31G*, PBE/def2-SVP on water; PySCF-optimized
  geometry near-stationary check).
- `tests/CMakeLists.txt` — register `test_cpu_scf_ks_dshell_gradient` under the
  `ATV_BUILD_RHF` block.

## Tasks

- [ ] Write failing d-shell KS-gradient FD test reproducing the ~1 eV/Bohr defect (tests/cpu/scf/test_ks_dshell_gradient.cpp)
- [ ] Register test_cpu_scf_ks_dshell_gradient in tests/CMakeLists.txt under the ATV_BUILD_RHF block
- [ ] Diagnose the d-shell XC-gradient term in accumulate_xc_gradient via per-term FD comparison (src/cpu/scf/gradient.cpp)
- [ ] Implement the d-shell XC-gradient fix in accumulate_xc_gradient (src/cpu/scf/gradient.cpp)
- [ ] Correct the d-shell poly_grad/poly_hess/eval_basis_all table in basis_eval.h if diagnosis implicates it, else leave unchanged
- [ ] Add docstring per doxygen documenting the diagnosed d-shell term and its units in gradient.cpp
- [ ] Verify analytic KS gradient matches FD to <= 1e-2 eV/A on water at PBE/6-31G* and PBE/def2-SVP
- [ ] Verify max-force at the PySCF-optimized water geometry drops to <= 1e-3 Ha/Bohr at PBE/def2-SVP
- [ ] Run full check + test suite

## Testing strategy

Happy path:

- `test_ks_dshell_gradient.cpp` builds water at 6-31G* and def2-SVP, runs PBE
  RKS to SCF tolerance <= 1e-9 Ha, computes the analytic KS gradient via
  `scf_gradient`, and compares every component against `scf_finite_difference_gradient`
  to <= 1e-2 eV/Å (the project's existing d-shell gradient tolerance, matching
  the `test_def2svp_gradient.cpp` d-shell tests).

Edge cases:

- Translational invariance: `sum_A dE/dR_A = 0` per axis to <= 1e-7 (carried over
  from `test_def2svp_gradient.cpp`).
- C2v symmetry of the water gradient: H1 and H2 force components mirror-equal.
- 6-31G* control vs def2-SVP: both d-bearing bases must pass, confirming the fix
  is not basis-data-specific.

Domain validation (`science.required`):

- Tight FD agreement: analytic vs central-difference FD <= 1e-5 Ha/Bohr per
  component on a fixed (non-rotating) dense grid with SCF converged to 1e-9 Ha,
  step h = 1e-3 Bohr, per the FD protocol in Domain basis. Achievable if the fix
  is complete; if grid noise prevents 1e-5 the criterion records the achieved
  bound and the 1e-2 eV/Å criterion remains the binding bar.
- Near-stationary recognition: at the PySCF-optimized water geometry the
  max-force component drops from ~1 eV/Bohr to <= 1e-3 Ha/Bohr at PBE/def2-SVP.
- Full `ctest` suite stays green — in particular the s/p gradient tests in
  `test_gradient.cpp`, the d-shell tests in `test_def2svp_gradient.cpp`, and the
  HF-gradient-vs-libcint test `test_gradient_vs_libcint.cpp` show no regression.

## Out of scope

- HF (exact-exchange) gradient — already correct per `test_gradient_vs_libcint.cpp`;
  not touched.
- f-shell (l=3) gradients — depend on the f-shell Rys ERI work in chain link 01
  (`phase1-productization-01-*`); this spec covers d-shells only.
- The denser SG-1/SG-2 integration grid — a separate Phase 1 sub-spec; this fix
  uses the existing grid and requires only that the same grid is used for the +
  and - FD displacements.
- Performance optimization of the GGA per-point Hessian loop — correctness only;
  any speedup is a follow-up.
