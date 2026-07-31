---
slug: gpw-vdw-multigrid
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: gpw / nonlocal vdW
queue: cp2k-method-queue
depends_on: []
grilled: false
---

# Spec: `gpw-vdw-multigrid` — nonlocal vdW-DF/DF2/rVV10 on the REL_CUTOFF multigrid

## Problem

`src/cpu/gpw/gpw_solver.cpp` rejects nonlocal vdW (`vdw_df` / `vdw_df2` /
`rvv10`) whenever the CP2K-style REL_CUTOFF multigrid product path is active
(`ngrids > 1 && rel_cutoff > 0`):

- the self-consistent V_nl is gated `&& !use_prod` inside `build_H_and_energy`
  (silently skipped on the product path, then …)
- … the post-SCF final E_c^nl block throws
  `"GPW <tag>: multigrid not supported (single-grid only)"` (was :5096).

CP2K's own default is `NGRIDS 4`, so every CP2K vdW input is a multigrid
input; ATV cannot twin the default until the gate opens.

## Why the gate can simply open

Nonlocal vdW is a functional of the **finest-grid** total density. On the
product path `build_H_and_energy` already prolongs every level's collocated
density to the finest grid *before* Poisson/XC, and it restricts the finest
`veff` back per level (exact adjoint) *after*. `nonlocal_vdw_on_grid` reads
that fine `rho` and adds `V_nl` into `veff` at exactly the point where
Hartree + XC potentials are summed — the multigrid transfer machinery is
orthogonal to it. The SCF therefore stays variational: the integrated KS
matrix remains the exact gradient of the collocated
E_H + E_xc + E_c^nl.

## Change (Γ-point closed-shell RKS only; k path stays rejected — queued `gpw-vdw-kpoints`)

1. **SC potential:** drop the `&& !use_prod` gate — apply
   `nonlocal_vdw_on_grid(grid, rho, …, potential=true)` on the fine grid in
   every SCF iteration on both single-grid and multigrid paths.
2. **Final E_c^nl:** delete the throw; collocate the final density via the
   product path (products per level + spectral prolongation) instead of the
   `cbg` gather, which is deliberately empty under multigrid. Extract the
   density-collocation branch (products vs cbg) into one lambda — it now has
   three call sites (`build_H_and_energy`, `fold_dispersion`, final E_c^nl).
3. **Found rot (fix in the same cut, fail-fast):** Γ-point **UKS** with any
   nonlocal vdW flag silently drops E_c^nl today (no SC potential, no
   post-SCF fold — `e_disp` stays 0 while the run "succeeds"). Mirror the
   k-path rule: reject at validation with a how-to. UKS vdW is a separate
   capability (spin-resolved fold), not this slug.

No public-surface change: `GpwDft` members `vdw_df` / `vdw_df2` / `rvv10`,
`ngrids`, `rel_cutoff` already exist; they merely become legal together.

## Oracle

CP2K 2026.1 (`cp2k.psmp`), offline via `tools/cp2k_oracle` (never in ctest).
Input extends the `gpw-vdw-df-ice` twin: water dimer, 12 Å periodic cell,
PBE, DZVP-MOLOPT-SR-GTH + GTH-PBE, CUTOFF 200 → `&MGRID NGRIDS 4
REL_CUTOFF 40` + `&NON_LOCAL TYPE DRSLL`. Pin `ENERGY| Total` and
`Dispersion energy` for NGRIDS 4 (gate) and NGRIDS 1 (context row).

## Tests

- Unit (`tests/cpu/gpw/`): validation-level only — UKS + vdW rejects
  fail-fast; multigrid + vdW is accepted by `build()`. No SCF solve in
  `tests/` (e2e forbidden).
- Regression `regressions/gpw/gpw-vdw-multigrid/`: Γ RKS water dimer,
  `vdw_df=true`, `ngrids=4`, `rel_cutoff=40`:
  - internal: |E_tot(mg) − E_tot(sg)| tight (measured, 1e-5 Ha class or
    better at the regression eps_default);
  - pin: E_disp(mg) vs CP2K NGRIDS-4 DRSLL `Dispersion energy` (1e-4 Ha
    class, matching `gpw-vdw-scf-potential`; tighten to what is measured).

## Out of scope

- k-point path vdW (queued `gpw-vdw-kpoints`)
- UKS vdW (fail-fast reject only)
- vdW forces/stress × multigrid gates (force path unchanged here)
- CUDA
