# Spec — `cpmd-uks-kpoints`

**Slug:** `cpmd-uks-kpoints`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Branch:** `wf/cpmd-uks-kpoints`  
**Method:** Multi-k UKS Car–Parrinello MD (dual-spin Bloch orbital set
\(C_k^\sigma\) propagated per spin and per sampled k)

## Problem

CPMD is product on Γ RKS (`cpmd-water` / `cpmd-rattle-water`), Γ UKS
(`cpmd-uks`, dual real α/β `Cpmd` propagators) and multi-k RKS
(`cpmd-kpoints`, one `CpmdCplx` per sampled k). The crossing is a declared
residual: `GPWEvaluator::build` throws
`"multi-k UKS + cpmd is not productized"` (`src/cpu/gpw/gpw_evaluator.cpp:242`
at spec time). Underneath, the solver has no multi-k **UKS** fixed-density
evaluation either: `run_gpw_scf_kpoints_uks` has no `max_iter == 0` branch,
and `GpwResult` has no `H_ks_k_beta` payload — so even a hand-rolled dual-spin
Bloch CP loop cannot get \(H_k^\beta\) at a supplied \(\{P_k^\alpha,
P_k^\beta\}\).

## Science

Same references as the two parents, crossed:

- Car & Parrinello, *Phys. Rev. Lett.* **55**, 2471 (1985) — fictitious-mass
  electronic dynamics.
- H. C. Andersen, *J. Comput. Phys.* **52**, 24 (1983) — RATTLE (electronic
  orthonormality constraint, per spin and k: \(C_k^{\sigma\dagger} S_k
  C_k^\sigma = I\)).
- Marx & Hutter, *Ab Initio Molecular Dynamics* (CUP) — CP conserved energy.
  For collinear spin on a k-mesh:
  \(E_{\mathrm{cons}} = T_n + T_e + E_{\mathrm{KS}}\) with
  \(T_e = \sum_k w_k \sum_{\sigma} \tfrac12 \mu \lVert \dot C_k^\sigma
  \rVert_F^2\).
- Spin-diagonal EOM (occupancy 1 per spin-orbital, mirroring the Γ-UKS
  product): \(\mu \ddot C_k^\sigma = -(I - S_k C_k^\sigma
  C_k^{\sigma\dagger})(2 H_k^\sigma C_k^\sigma) + \text{RATTLE}\).
  \(H_k^\sigma\) is the spin-resolved KS matrix at the **fixed** dual-spin
  k-density \(\{P_k^\alpha, P_k^\beta\}\) rebuilt from the propagated
  coefficients (occ = 1).
- Forces: Hellmann–Feynman at fixed density via `gpw_forces_from_scf`'s
  multi-k UKS path (`C_occ_k`/`eps_k` + `C_occ_k_beta`/`eps_k_beta`, occ 1
  per spin) with \(\varepsilon_i \approx \mathrm{Re}(C^\dagger H C)_{ii}\)
  per spin/k (same approximation the multi-k RKS CPMD branch ships).

## Implementation

1. **`GpwResult::H_ks_k_beta`** (`gpw_solver.h`): multi-k UKS β fixed-density
   KS matrices, paired with `H_ks_k` (α), size `n_kpts`; empty otherwise.
2. **`run_gpw_scf_kpoints_uks` fixed-density branch** (`gpw_solver.cpp`):
   `cfg.max_iter == 0` → require both `guess_P_k` (α) and `guess_P_k_beta`
   (β); refuse `use_ot`; copy the guesses into `Pa_k`/`Pb_k` on every rank;
   one collocate→Poisson→spin-XC→integrate energy build (hoist the existing
   `build_uks_k_energy` lambda out of the OT block — second caller); rebuild
   full \(H_k^\alpha / H_k^\beta\) on every rank from the shared spin V_eff;
   fill `H_ks_k` / `H_ks_k_beta` / `S_k` / `P_k` / `P_k_beta` / `k_frac` /
   `k_weights`; integer ⟨S²⟩; `converged = true`, `iters = 0`; fold
   dispersion. Mirror of the multi-k RKS `max_iter == 0` branch.
3. **`GPWEvaluator`** (`gpw_evaluator.{h,cpp}`):
   - members: `cpmd_k_beta_` (`std::vector<electronic::CpmdCplx>`) +
     `cpmd_guess_p_k_beta_` scratch, siblings of the RKS-k pair;
   - `build()`: delete the multi-k UKS refusal; keep `n_alpha >= 1` and the
     multi-k **thermostat** refusal (NVE electrons only on the k path —
     `cpmd-thermostat-kpoints` is a separate queue row); clear the β k-set
     wherever `cpmd_k_` is cleared;
   - `launch()` warm branch: multi-k readiness now also requires every β
     channel ready when `n_beta > 0`; standalone `overlap_k` → per-spin
     `align_to_overlap`; dual-spin density rebuild (occ 1) → fixed-density
     eval with both k-guesses; energy publish; forces populate
     `C_occ_k`/`eps_k` (α from \(H_k^\alpha\)) and
     `C_occ_k_beta`/`eps_k_beta` (β from \(H_k^\beta\)), `converged = true`,
     `gpw_forces_from_scf`; advance α on \(H_k^\alpha\) and β on
     \(H_k^\beta\) (occ 1 each); \(T_e = \sum_k w_k (T_e^{\alpha,k} +
     T_e^{\beta,k})\) published on the frame;
   - seed branch: after the first converged multi-k UKS SCF, seed one
     `CpmdCplx` per (k, σ) from `C_k`/`C_occ_k` (α, `n_alpha` columns) and
     `C_k_beta`/`C_occ_k_beta` (β, `n_beta` columns). `n_beta == 0`
     half-metal edge: no β propagators, zero β guess densities (mirrors the
     Γ-UKS edge).
4. **Unit tests** (`tests/cpu/gpw/test_gpw_cpmd.cpp`): flip
   `BuildRejectsMultiKUksCpmd` → `BuildAcceptsMultiKUksCpmd`; solver-level
   fixed-density dual-spin export test (tiny inline H2, `n_alpha=n_beta=1`,
   1×1×2): `H_ks_k`+`H_ks_k_beta`+`S_k`+`P_k`+`P_k_beta` sized `n_kpts`,
   `iters == 0`, energy ≈ cold SCF at the same density; Γ-limit identity:
   the k-path fixed-density eval at one explicit Γ point (`k_points =
   {{0,0,0}}`) reproduces the Γ dual-spin fixed-density eval (energy, and
   \(H_k^\sigma\) vs real \(H^\sigma\)) — anchors the new branch to the
   shipped `cpmd-uks` product surface.

No public-surface change beyond deleting a refusal: `CarParrinello{gpw, mu}`
with `gpw.n_alpha/n_beta` + `gpw.k_mesh` is the whole user story (no Config
bag, region-free, Frame in/out).

## Regression (`regressions/dynamics/cpmd-uks-kpoints/`)

Family crossing of `cpmd-kpoints` × `cpmd-uks`, geometry twin of
`xl-bomd-uks-kpoints`: periodic H₂ (R = 1.0 Bohr, L = 12 Bohr cube,
DZV-GTH-PADE, GTH-PADE-q1, PADE), mesh 1×1×2 unreduced, **closed-shell
fillings on the UKS path** (`n_alpha = n_beta = 1`) so both spin channels
propagate and the molecule stays bound (triplet H₂ is repulsive and cannot
pin NVE — the xl-bomd-uks-kpoints precedent). μ = 400 a.u., dt = 0.1 fs,
1 seed + 1 baseline + 4 production steps.

Gates: CP conserved-energy relative drift ≤ 0.10 (family envelope);
production `iters == 0` with `n_kpts >= 2`; `T_e > 0` once propagated;
UKS-vs-RKS closed-shell consistency — dedicated multi-k **UKS** and
multi-k **RKS** single points (`cpu::single_point`, plain `GPWEvaluator`,
same geometry/mesh/cutoff, `tol = 1e-7`) agree within 3×10⁻⁴ eV (the
`gpw-kpoints-uks` A2 identity class); finite coords; H–H window.

## Oracle / pin provenance

Queue-assigned family policy (`cpmd-kpoints` / `cpmd-uks`): **literature
NVE-drift envelope**, `perf_exempt`. Stated gate, not "no binary": CP2K
2026.1 builds and runs on this host (`gpw-bomd-kpoints` pinned against it),
but CP2K has **no fictitious-mass Car–Parrinello propagator** (its
`RUN_TYPE MD` is BOMD; its CP-like method is second-generation CPMD) — the
CP conserved quantity \(T_n + T_e + E_{\mathrm{KS}}\) with \(T_e\) of the
fictitious \(\dot C_k^\sigma\) is not a number CP2K can print for any input.
The BO surface underneath is CP2K-anchored by the shipped multi-k UKS pins
(`gpw-kpoints-uks`, `gpw-kpoint-uks-forces`) and the closed-shell RKS
crossing above.
