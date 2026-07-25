---
slug: gpw-kpoints-si
status: approved
created: 2026-07-25
revised: 2026-07-25
domain: gpw / periodic DFT
queue: cp2k-method-queue
deferred: false
depends_on: []   # Γ energy/forces pins all shipped (waves 1-2)
unblocks: [gpw-kpoints-si]
---

# Spec: `gpw-kpoints-si` — k-point GPW RKS (Quickstep k-points, Si)

**Slug:** `gpw-kpoints-si`
**Domain:** GPW Kohn–Sham DFT under PBC (CPU)
**Queue row:** `k-point GPW RKS | gpw | Si | DZVP-MOLOPT-GTH + GTH-PBE | e_tot Ha`
**Twin:** CP2K Quickstep `&KPOINTS SCHEME MONKHORST-PACK` (offline pin only)

## Summary

The GPW solver (`cpu::gpw::GpwDft`) is Γ-only: real Bloch-summed S/T/V_pp
lattice sums, real generalized eigensolve, real density collocation. This
spec adds **Monkhorst–Pack k-point sampling for closed-shell RKS**: complex
Bloch-phased one-electron / nonlocal-PP matrices per k, per-k complex
Hermitian generalized diagonalization, k-weighted density accumulation on the
real-space grid, and a k-summed Pulay DIIS SCF — the CP2K Quickstep k-point
loop shape (real-space operators → per-k Fourier phase → complex eigensolve →
k-summed density).

Target system for the regression: **Si diamond in the 4-atom
body-centred-tetragonal cell** (orthorhombic-compatible: the solver's cell is
`double cell[3]`), DZVP-MOLOPT-GTH + GTH-PBE-q4, PBE XC — all already in the
in-tree tables (`gpw_tables.h`, shipped by `gpw-pbesol-si`).

## Domain formulation

Bloch AO basis per k (fractional κ, `k = 2π diag(1/L) κ`):

    χ_μ^k(r) = Σ_T e^{i k·T} φ_μ(r − T),   T = (t₁L₁, t₂L₂, t₃L₃)

- **Analytic parts** (one image sweep fills every k; per-image real integral ×
  phase): `S_k(μν) = Σ_T e^{ik·T} ⟨φ_μ|φ_ν(·−T)⟩`, `T_k` likewise
  (McMurchie–Davidson, unchanged); GTH nonlocal via phased Kleinman–Bylander
  projector overlaps `O_μ^k[A,slot] = Σ_T e^{ik·T}⟨φ_μ(·−T)|p_A⟩`,
  `V_nl^k = Σ_A O^k† h O^k`.
- **Grid parts** are k-independent operators: V_H, V_xc, V_loc^SR live on the
  real-space grid built from the k-summed density
  `ρ(r) = Σ_k w_k Σ_μν P^k_μν χ^k_μ(r) χ^{k*}_ν(r)` (real by Hermiticity of
  P^k). Matrix elements back per k: `V^k_μν = ∫ χ^{k*}_μ V_eff χ^k_ν`.
  Implementation: per-k **pre-summed** complex basis-on-grid values
  χ^k_μ(g) (per active grid point) — the image sum with phases is folded at
  build time, so the per-iteration collocate/integrate cost is O(n_pts·n_b²)
  per k with no per-image duplication (the Γ path's duplicated-triple CSR
  does not scale to diffuse MOLOPT in a 7 Bohr cell).
- **Eigensolve:** per-k canonical orthogonalization of complex Hermitian S_k
  (drop eigenvalues < `lindep_threshold`), then zheevd on X_k† H_k X_k.
  Direct LAPACK `zheevd_` binding with explicit workspace (the bundled
  LAPACK's LWORK=-1 query is broken — same rationale as `dsygvd_` in
  `linalg/eigensolve.h`).
- **SCF:** per-k KS matrices share one Pulay DIIS extrapolation
  (`B_ij = Σ_k w_k Re Tr(e_i^k† e_j^k)` — the standard k-summed DIIS);
  convergence = |ΔE| < tol AND max_k ‖[H_k,P_k]_{S_k}‖_∞ < 1e-5 (the Γ-path
  gates, k-generalized). Occupations: fixed integer (insulator), lowest
  `n_occ` bands doubly occupied at every k.
- **Monkhorst–Pack mesh** (CP2K `SCHEME MONKHORST-PACK N N N` convention):
  κ_d = (2r − N_d − 1)/(2N_d), r = 1..N_d (even N ⇒ Γ-avoiding shifted mesh,
  e.g. N=2 → ±1/4 → genuinely complex phases {±1, ±i}). Optional
  time-reversal reduction k ↔ −k (weights doubled; ρ uses Re, so TR is
  exact for RKS).
- **Energy:** E = Σ_k w_k Re[Tr(P_k T_k) + Tr(P_k V_nl^k)] + ∫ρV_loc^SR
  + (E_H − E_self) + E_xc (+ E_disp). Reported in Hartree on `GpwResult`
  exactly as the Γ path.

## Public surface (members on the concrete method — no Config bag)

`cpu::gpw::GpwDft` gains:

- `int k_mesh[3] = {1, 1, 1};` — MP mesh (CP2K MONKHORST-PACK convention).
  `{1,1,1}` keeps the existing real Γ path bit-unchanged.
- `bool k_time_reversal = true;` — fold ±k pairs (CP2K default behaviour).
- `std::vector<std::array<double, 3>> k_points;` +
  `std::vector<double> k_point_weights;` — explicit fractional k list
  (CP2K `SCHEME GENERAL` / explicit `KPOINT` lines). Non-empty overrides
  `k_mesh`. Needed by the supercell-folding acceptance gate and later band
  paths.

`cpu::gpw::GpwResult` gains (result payload, not config):

- `std::vector<std::array<double, 3>> k_frac;` — sampled fractional k.
- `std::vector<double> k_weights;`
- `std::vector<std::vector<double>> eps_k;` — band energies per k [Ha].

Scope guards (fail-fast throw, no fallback): the k path is RKS + periodic +
`local_pp_on_grid` + single grid only — UKS, OT, DFT+U, multigrid,
`nonperiodic`, analytic local-PP, and `level_shift != 0` throw
`std::invalid_argument`/`runtime_error` when combined with a non-Γ k
configuration. D3/D3(BJ)/D4 stay allowed (geometry-only additive).
`gpw_forces_fd` / `gpw_stress` propagate the new members (`copy_gpw`).

## Internals

`include/atomiverse/cpu/gpw/internal/kpoints.h` + `src/cpu/gpw/internal/kpoints.cpp`
(never on the public include path):

- `struct KPoint { double frac[3]; double weight; };`
- `std::vector<KPoint> monkhorst_pack(const int mesh[3], bool time_reversal);`
- `overlap_kinetic_k(...)` — one image sweep filling S_k, T_k for all k.
- `add_nonlocal_kb_k(...)` — phased KB assembly (mirrors `add_nonlocal_kb`).
- `KpointBasisGrid` — per-k pre-summed complex χ^k on active grid points;
  `collocate` (adds w_k ρ_k) and `integrate` (V^k from V_eff).
- `diagonalize_hermitian_generalized(H_k, S_k, lindep)` — canonical-orth +
  `zheevd_` (explicit workspace).
- k-SCF loop hosted in `gpw_solver.cpp` beside the UKS/OT branches, sharing
  the Γ setup (grid, n_ion, v_locsr, e_self, dispersion).

## Out of scope (follow-ons)

Symmetry (IBZ) reduction beyond time reversal; k-point UKS; analytic k-point
forces/stress (FD works through the public surface); band-structure paths /
DOS output; smearing/fractional occupations; multigrid + k; OT + k.

## Oracle / pin

No CP2K binary or module on the worker host (verified 2026-07-25:
`module avail cp2k` empty, no `cp2k*` on PATH) — same status as the
`gpw-pbesol-si` / `gpw-d4-water` pins. Strategy (binding, mirrors precedent):

1. **Internal exactness gates** (the decisive k-point correctness checks —
   they hold to SCF tolerance regardless of any oracle):
   supercell folding identity + time-reversal identity + Γ-limit identity
   (see acceptance.md).
2. **Frozen first-green ATV absolute band** for e_tot (setup integrity, not a
   third-party oracle).
3. Twin input `tools/cp2k_oracle/inputs/gpw-kpoints-si.inp` committed so the
   absolute pin can be replaced by a cited CP2K number the moment an offline
   CP2K becomes available.

## Tests

- `tests/cpu/gpw/test_cpu_gpw_kpoints.cpp` — unit, single-function:
  MP mesh generation (count/weights/CP2K shift/TR fold), S_k(Γ) == real
  periodic overlap, Hermiticity + S_{−k} = conj(S_k), zheevd wrapper on an
  analytic 2×2 complex Hermitian pair, KpointBasisGrid Γ-limit vs
  `collocate_density`, grid commensuration for the supercell gate cutoffs.
- Regression `regressions/gpw/gpw-kpoints-si/reg_gpw-kpoints-si.cpp`
  (public headers only, flat main, ≤ 60 s, hard-coded numbers; name follows
  the `atv_add_regression(<domain> <case>)` → `reg_<case>.cpp` convention).
