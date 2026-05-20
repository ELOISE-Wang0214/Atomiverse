---
title: KS XC grid-response nuclear gradient (Becke weight + grid-coordinate terms)
status: approved
created: 2026-05-20
---

# KS XC grid-response nuclear gradient (Becke weight + grid-coordinate terms)

## Summary
The CPU Kohn-Sham analytic nuclear gradient currently computes the exchange-correlation
contribution `dE_xc/dR_A` on a grid held fixed in space and then projects out the residual
net-translational quadrature error (~4e-4 Ha/Bohr) to restore translational invariance. This
makes *relative* forces exact but leaves the *absolute* per-atom XC gradient resting on a
projection shortcut rather than first principles. This spec implements the genuine
grid-response term — the contribution from the nuclear-coordinate dependence of the Becke
fuzzy-cell quadrature weights and from grid points riding rigidly with their parent atoms —
so the absolute gradient is correct from first principles and `sum_A dE_xc/dR_A` vanishes
*without* any projection. The implementation is a verbatim port of PySCF's
`grids_response_becke`; once verified, the translational-invariance projection at
`gradient.cpp:447-469` is removed, with an optional debug gate left behind. Scope is CPU
RKS/UKS with LDA and GGA functionals only.

## Domain basis
The total KS XC nuclear gradient decomposes into three additive pieces:

```
dE_xc/dR_A = (AO-derivative term)          # already implemented, accumulate_xc_gradient
           + (weight-derivative term)      # NEW — Becke partition response
           + (grid-coordinate term)        # NEW — grid rides with parent atom
```

**(1) Weight-derivative term.** A grid weight factorises as `W_g = w_radial * w_angular *
w_n(R)`, where only the Becke fuzzy-cell partition weight `w_n` carries nuclear-position
dependence. Its contribution is

```
dE_xc/dR_A|weight = sum_g (dw_g/dR_A) * eps_xc(rho(r_g)) * rho(r_g)
```

The `* rho` factor is mandatory here: Atomiverse stores `exc` as XC energy **per electron**
(`xc_functional.h:22`; total energy is assembled as `E_xc = sum_g w * exc * (rho_a + rho_b)`
in `xc_grid_engine.h:177`). This matches PySCF's convention `einsum('r,r,nxr->nx', exc, rho,
weight1)`, so the port carries the `rho` factor verbatim. The authoritative reference is
`grids_response_becke` in `pyscf/grad/rks.py`; its core (verbatim):

```
def get_du(ia, ib):                       # JCP 98, 5612 (1993), eq. (B10)
    uab = atm_coords[ia] - atm_coords[ib]
    duab  = 1./atm_dist[ia,ib] * grid_norm_vec[ia]
    duab -= uab[:,None]/atm_dist[ia,ib]**3 * (grid_dist[ia]-grid_dist[ib])
    return duab
# for each atom pair (ia, ib), ib < ia:
g     = 1/atm_dist[ia,ib] * (grid_dist[ia]-grid_dist[ib])
g_adj = fadjust(ia,ib,g)                  # size adjustment nu = mu + a(1-mu^2)
smoothed, t_uab = get_smooth_function_derivative(g_adj)
t_uab *= 0.5 * gadjust(ia,ib,g)           # gadjust = 1 - 2*a*mu
s_uab = .5*(1 - smoothed); s_uba = .5*(1 + smoothed)
pbecke[ia] *= s_uab; pbecke[ib] *= s_uba
pt_uab = -t_uab/s_uab; pt_uba = t_uab/s_uba
# scatter pt_uab/pt_uba * get_du(...) into dpbecke[natm,natm,3,ngrid], including
# the extra block dpbecke[atom_id,ia/ib] -= pt_uab*(ua-ub)/R_ab for ia,ib != atom_id
dpbecke[:,ia] *= pbecke[ia]               # finally
z  = 1/pbecke.sum(0)                      # normalized-weight derivative
w1 = dpbecke[:,ia]*z - pbecke[ia]*z**2 * dpbecke.sum(1)
w1 *= vol
```

For the `original_becke` switching scheme, the smoothing function and its derivative (from
the PySCF C kernel `grid_basis.c`, verbatim):

```
# iterated polynomial f(x) = (3x - x^3)/2 applied 3x to nu = mu + a(1-mu^2)
get_smooth_function_derivative(mu):
    p0 = mu; p1 = (3-p0^2)*p0*.5; p2 = (3-p1^2)*p1*.5; p3 = (3-p2^2)*p2*.5
    t_uab = 27/8 * (1-p2^2)*(1-p1^2)*(1-p0^2)   # = dp3/dmu without the -1/2 factor
    return p3, t_uab
s(mu) = .5*(1 - p3);   ds/dmu = -.5 * t_uab
```

Size adjustment: `nu = mu + a_ij*(1-mu^2)`, with `a_ij` from the same atomic-radii table the
grid builder uses, clamped to `[-1/2, 1/2]`; `fadjust(i,j,g) = g + a*(1-g^2)`;
`gadjust(i,j,g) = 1 - 2*a*g`. The confocal-elliptical coordinate is
`mu_ij = (r_i - r_j)/R_ij`, with `r_i` the grid-point-to-atom-i distance; the chain rule
`dmu/dR_A` runs through `r_i` (`dr_i/dR_A = -unit_vec * delta_iA`) and through `R_ij`. The
self term `dw_A/dR_A` is computed indirectly via translational invariance
`dw_A/dR_A = -sum_{B!=A} dw_A/dR_B` because `mu -> 0` on an atom causes numerical blow-up;
this PySCF convention is ported faithfully.

**(2) Grid-coordinate term.** Each grid point rides rigidly with its parent atom P:
`dr_g/dR_P = I`, `dr_g/dR_A = 0` for `A != P`. PySCF does not differentiate this spatially;
it reuses the AO-derivative-contracted XC potential matrix (`vtmp`, already assembled for the
standard XC gradient) contracted with the density matrix and attributed to the parent atom:

```
excsum[parent_atom] += einsum('xij,ji->x', vtmp, dm) * 2
```

This term is nearly free once the AO-gradient-contracted matrix exists.

**(3) Convention risk.** See item (1): the `* rho` factor is load-bearing and a mismatch is
a silent factor-of-rho bug invisible to any single-geometry energy check. Atomiverse's `exc`
being per-electron has been confirmed against `xc_functional.h` and `xc_grid_engine.h`; the
port keeps the `rho` factor.

References (PySCF source is the binding reference; the JGP appendix B8/B10 are noted by
PySCF itself as incomplete, so they are cited as background only):
- PySCF source: `pyscf/grad/rks.py` — `grids_response_becke`, `get_vxc_full_response`.
- PySCF source: `pyscf/lib/dft/grid_basis.c` — `VXCgen_grid`, `VXCbecke_weight_derivative`.
- A. D. Becke, J. Chem. Phys. 88, 2547 (1988), DOI 10.1063/1.454033 — background.
- B. G. Johnson, P. M. W. Gill, J. A. Pople, J. Chem. Phys. 98, 5612 (1993),
  DOI 10.1063/1.464906 — background (appendix B8/B10 incomplete per PySCF).

## Design
- **New grid helper (`grid.h` / `grid.cpp`).** Add `becke_weight_derivative(atoms, natm,
  points, dw_out)` — given the molecule and the existing flat grid points, it returns
  `dw_g/dR_A` as a `[natm, ngrid, 3]` (or equivalent flat) buffer, mirroring
  `grids_response_becke`. It internally recomputes the Becke geometry (per-point per-atom
  distances, unit vectors, `atm_dist`) and the size-adjusted `a_ij` from the same atomic-radii
  table `becke_partition_weights` uses, so consistency with the production weights is
  structural, not coincidental. The helper also exposes the per-point parent-atom id (needed
  by the grid-coordinate term). `gradient.cpp` stays a pure assembly site: it does not
  re-derive Becke geometry.
- **Weight-derivative accumulation (`gradient.cpp`).** Inside `accumulate_xc_gradient`
  (`gradient.cpp:144`), after `exc_arr` and `rho_a/rho_b` are evaluated, contract
  `dw_g/dR_A` with `eps_xc(g) * rho(g)` to add `dE_xc/dR_A|weight` into `grad_xc`.
- **Grid-coordinate accumulation (`gradient.cpp`).** Reuse the AO-gradient-contracted XC
  potential already assembled in `accumulate_xc_gradient` (the `WT_a/WT_b` /
  per-point GGA contraction in section 5), attribute its contracted value to each grid
  point's parent atom, and add it into `grad_xc`. If the existing assembly does not expose a
  per-grid-point AO-contracted quantity in a reusable form, factor it out so both the
  standard AO-derivative term and this term consume it.
- **Projection removal (`gradient.cpp:447-469`).** Once both new terms are verified, replace
  the unconditional translational-invariance projection with a no-op by default, gated behind
  a new `SCFConfig` debug flag (`xc_grad_project_translation`, default `false`). The genuine
  net translational force must come out `<= 1e-8` Ha/Bohr per axis without it.
- **Self-term convention.** The grid helper computes `dw_A/dR_A` as
  `-sum_{B!=A} dw_A/dR_B` exactly as PySCF does; no direct `mu -> 0` differentiation.
- **Spin handling.** The weight-response term uses total `rho = rho_a + rho_b` and the
  per-electron `exc`; it is spin-agnostic and identical for RKS and UKS (the existing
  `xc_spin_scale` logic at `gradient.cpp:213-222` governs only the density fed to the
  functional, not the weight response).

## Files to create or modify
- `include/atomiverse/cpu/scf/grid.h` — declare `becke_weight_derivative` (and a parent-atom
  accessor if not already derivable); Doxygen with units.
- `src/cpu/scf/grid.cpp` — implement `becke_weight_derivative` as a verbatim port of
  `grids_response_becke` plus `get_smooth_function_derivative`.
- `src/cpu/scf/gradient.cpp` — add weight-derivative and grid-coordinate accumulation in
  `accumulate_xc_gradient`; gate the projection at lines 447-469 behind the new debug flag.
- `include/atomiverse/cpu/scf/scf_solve.h` — add the `xc_grad_project_translation` debug
  flag to `SCFConfig` (default `false`).
- `tests/cpu/scf/test_ks_grid_response_gradient.cpp` (new) — translational-invariance,
  moving-grid FD, and projection-removal tests.
- `tests/cpu/scf/test_ks_dshell_gradient.cpp` — adjust only if a fixed-grid test there
  implicitly depended on the projection; otherwise leave untouched (regression target).
- `tests/CMakeLists.txt` — register the new test target.

## Tasks
- [ ] Write failing tests for `becke_weight_derivative` correctness (tests/cpu/scf/test_ks_grid_response_gradient.cpp) — per-point `dw_g/dR_A` against a central-difference of `becke_partition_weights_for_test`, and the partition-sum identity `sum_A dw_g/dR_A == 0`
- [ ] Implement `becke_weight_derivative` and `get_smooth_function_derivative` in src/cpu/scf/grid.cpp, declared in include/atomiverse/cpu/scf/grid.h — verbatim port of PySCF `grids_response_becke`
- [ ] Add Doxygen to `becke_weight_derivative` per doxygen style with all parameter units (Bohr, Hartree/Bohr) and a note that `dw_A/dR_A` follows PySCF's translational-invariance convention
- [ ] Write failing tests for the assembled XC grid-response gradient: zero net translational force on water PBE/def2-SVP, and moving-grid central-difference FD parity (tests/cpu/scf/test_ks_grid_response_gradient.cpp)
- [ ] Implement weight-derivative and grid-coordinate accumulation in `accumulate_xc_gradient` in src/cpu/scf/gradient.cpp, contracting `dw_g/dR_A` with `exc(g)*rho(g)` and attributing the AO-contracted grid-coordinate term to each parent atom
- [ ] Add `xc_grad_project_translation` debug flag (default `false`) to `SCFConfig` in include/atomiverse/cpu/scf/scf_solve.h and gate the projection block at gradient.cpp:447-469 behind it
- [ ] Register tests/cpu/scf/test_ks_grid_response_gradient.cpp in tests/CMakeLists.txt
- [ ] Verify against PySCF reference: moving-grid FD parity for PBE/def2-SVP and PBE/6-31G* on water to <= 1e-5 Ha/Bohr per component
- [ ] Run full check + test suite

## Testing strategy
Happy path:
- Translational invariance — with the projection disabled (`xc_grad_project_translation =
  false`), `||sum_A dE_xc/dR_A||` per axis on water PBE/def2-SVP is `<= 1e-8` Ha/Bohr. This
  needs no external reference data and is the core first-principles diagnostic.
- Moving-grid FD parity — the analytic KS gradient matches a central-difference FD that
  *regenerates the grid at each displaced geometry* (`h = 1e-3` Bohr, SCF tol `1e-9`), to
  `<= 1e-5` Ha/Bohr per component, for water PBE/def2-SVP and PBE/6-31G*. This is distinct
  from the existing fixed-grid FD test in `test_ks_dshell_gradient.cpp`; the moving-grid FD
  is precisely what the projection shortcut cannot pass.

Edge cases:
- `becke_weight_derivative` per-point unit check — `dw_g/dR_A` matches a central difference
  of `becke_partition_weights_for_test` to tight tolerance, and `sum_A dw_g/dR_A == 0` per
  grid point (rigid-shift identity), exercised on a multi-atom geometry.
- Both LDA (`xc_type` LSDA) and GGA (`xc_type` PBE) functional paths covered by the
  moving-grid FD test; RKS and UKS both exercised on water.

Regression / domain validation:
- The existing fixed-grid tests in `test_ks_dshell_gradient.cpp` still pass (with the
  projection re-enabled via the debug flag where a test relies on the old fixed-grid
  convention).
- The full `ctest` suite stays green.
- Domain validation — the moving-grid FD parity itself is the physics check: it confirms the
  ported Becke weight response reproduces the true `dE_xc/dR` against an independent numerical
  oracle, guarding against the factor-of-rho convention bug and against the ~30 eV/Bohr
  garbage-error class that sank the prior attempt.

## Out of scope
- f-shell gradients (the grid-response port is independent of basis angular momentum, but
  f-shell AO derivatives remain out of scope per the predecessor spec).
- CUDA backend — CPU RKS/UKS only.
- meta-GGA tau-dependent grid response.
- NLC / VV10 grid response.
- SG-1 / SG-2 grid construction work (already complete).
- The `treutler` switching-function variant for the weight derivative — only
  `original_becke` is ported, matching the production grid builder; revisited only if the
  grid builder's switching scheme changes.
