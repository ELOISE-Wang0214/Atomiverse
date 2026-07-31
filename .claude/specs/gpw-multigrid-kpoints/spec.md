---
slug: gpw-multigrid-kpoints
status: approved
created: 2026-07-31
revised: 2026-07-31
domain: gpw / k-point sampling / multigrid
queue: cp2k-method-loop
deferred: false
depends_on: []
unblocks: [gpw-multigrid-kpoints, gpw-vdw-multigrid]
grilled: false
---

# Spec: `gpw-multigrid-kpoints` — Multi-k GPW with REL_CUTOFF multigrid collocation

**Slug:** `gpw-multigrid-kpoints`
**Domain:** GPW k-point SCF (CPU), CP2K CUTOFF/REL_CUTOFF multigrid
**Queue row:** `gpw-multigrid-kpoints` (GapCensus wave 14; assigned wave 1 2026-07-31)
**Oracle:** CP2K 2026.1 (`cp2k.psmp`, module `CP2K/2026.1` on Arrhenius), offline pin

## Summary

`GpwDft::build` rejects `ngrids > 1` on the k-point path
(`gpw_solver.cpp` k-gate: `reject("multigrid (ngrids > 1)")`), so multi-k
GPW always pays the single fine grid: diffuse MOLOPT products span huge
windows on the CUTOFF grid and the discretization differs from CP2K's
NGRIDS 4 default. Γ RKS and Γ/UKS multigrid product collocation already
shipped (`internal/multigrid.h`: `build_products` / `collocate_products` /
`integrate_products` + `prolong_fft` / `restrict_fft`). This spec extends the
**product machinery** to Bloch sampling so the closed-shell multi-k RKS SCF
(DIIS and OT, with smearing and DFT+U riding the same chokepoints) converges
with `ngrids > 1 && rel_cutoff > 0`.

## Domain basis

The Γ product expansion already enumerates primitive-pair products with the
**ket image-shifted** (`rb = B.center + tL`, bra home) and collocates over the
unwrapped window with wrapped scatter — which supplies the bra-image lattice
sum. Under Bloch sampling the same geometric product list serves all k; only
the density weight and the KS scatter acquire phases:

- **Density:** with χ_μ^k per `internal/kpoints.h` conventions,
  ρ(g) = Σ_t Σ_μν D_μν(t) · [periodized product (μ home, ν at ket shift t)],
  where **D(t) = Σ_k w_k Re[P^k e^{−2πi κ·t}]** is the real-space density
  matrix stack over the distinct ket image shifts t appearing in the product
  list. At Γ, D(t) = P for every t — the shipped Γ weight.
- **KS integration:** per-level potential integrals accumulate into real
  per-cell matrices **V_t(μν) = Σ_g w V(g) φ_μ(r_g) φ_ν(r_g − tL)**
  (the exact transpose of the collocation); the per-k grid KS matrix is the
  phase sum **V^k = Σ_t e^{+2πi κ·t} V_t** (Hermitian because
  V_{−t} = V_tᵀ by the lattice-shifted-window symmetry).
- Level transfer (prolong to the finest grid for Poisson + XC, adjoint
  restrict back per level) is unchanged from Γ multigrid — the density on the
  fine grid is real regardless of k.

Time-reversal folding stays exact: D and V_t use only Re[P e^{−ik·t}] /
real V_t, and P^{−k} = conj(P^k).

## Deliverables

1. `internal/multigrid.h/.cpp`: `GridProduct` records its ket image shift;
   `ProductSet` records the distinct-shift cell list; new
   `collocate_products_cells` (per-cell D stack weight) and
   `integrate_products_cells` (per-cell V_t stack out) siblings of the Γ
   functions. Γ entry points unchanged (bit-identical single-density path).
2. `gpw_solver.cpp` `run_gpw_scf_kpoints` (RKS): under
   `ngrids > 1 && rel_cutoff > 0` skip the `KpointBasisGrid` build; density
   via D-stack product collocation per level + prolong; per-k grid KS via
   per-level restrict + `integrate_products_cells` + phase sum. DIIS, OT,
   smearing, DFT+U, dispersion folds ride the same chokepoints
   (`build_k_energy` / `rebuild_all_H` / Hirshfeld fold).
3. `GpwDft::build` k-gate: accept `ngrids > 1` for RKS; multi-k **UKS**
   multigrid stays rejected fail-fast (scope: this row is the closed-shell
   observable; the UKS twin is a census residual).
4. Regression `regressions/gpw/gpw-multigrid-kpoints/` + `REFERENCE.md` with
   the offline CP2K 2026.1 pin (twin input
   `tools/cp2k_oracle/inputs/gpw-multigrid-kpoints/si_kpoints_mg.inp`).
5. Unit tests (`tests/cpu/gpw/test_cpu_gpw_multigrid.cpp`): collocate/integrate
   adjointness of the `_cells` pair on a toy grid; Γ-phase equivalence
   (D-stack of identical matrices ≡ Γ product path); D/V_t phases vs the
   independent `KpointBasisGrid` machinery; time-reversal D-stack folding.

## Acceptance

See `acceptance.md`.

## Out of scope

- Multi-k **UKS** multigrid (gate keeps throwing with a how-to).
- Analytic multi-k forces/stress consistency with the multigrid functional
  (force paths keep their own single-grid assembly, as Γ multigrid ships).
- Nonlocal vdW (vdw_df/vdw_df2/rvv10) under k-points — found silently
  dropped on the k-path during this spec; now rejected fail-fast at build
  (named debt, sibling row `gpw-vdw-multigrid` holds the multigrid part).
