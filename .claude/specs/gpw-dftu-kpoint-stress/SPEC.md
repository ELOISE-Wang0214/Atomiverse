# Spec — `gpw-dftu-kpoint-stress`

**Status:** acceptance_met (pending harness merge)  
**Slug:** `gpw-dftu-kpoint-stress`  
**Domain:** gpw  
**Oracle:** CP2K (offline pins only)  
**Branch:** `wf/gpw-dftu-kpoint-stress`

## Intent

Close the residual **multi-k Dudarev DFT+U analytic cell stress** gap left after:

| Slug | Capability |
|------|------------|
| `gpw-dftu-nio` / `gpw-dftu-stress` | Γ-UKS + PLUS_U energy + analytic stress |
| `gpw-dftu-kpoints` | Multi-k UKS + PLUS_U SCF energy |
| `gpw-kpoint-stress` / `gpw-kpoint-uks-stress` | Multi-k RKS/UKS analytic stress (no U) |

Today `energy_fixed_P_k` throws when `plus_u` is non-empty. Product path must
include \(E_U[n^\sigma(\{P_k^\sigma\},\{S_k\})]\) in the Hellmann–Feynman
fixed-\(\{P_k,W_k\}\) Lagrangian under orthorhombic diagonal strain.

## Science

Dudarev / simplified DFT+U (one spin channel):

\[
n^\sigma=\sum_k w_k\,\tfrac12\mathrm{Herm}\big[(P_k^\sigma S_k)+(S_k P_k^\sigma)\big],
\quad
E_U^\sigma=\tfrac{U}{2}\big[\mathrm{Tr}\,n^\sigma-\mathrm{Tr}(n^\sigma n^\sigma)\big]
\]

For the fixed-\(P\) cell theorem, \(\{P_k^\sigma\}\) (and \(W_k^\sigma\) from
eigenvalues that already include \(V_U\)) are frozen from the reference SCF;
only geometry-dependent \(S_k\), \(T_k\), \(V_{\mathrm{nl}}^k\), grid, Poisson,
XC, and therefore \(n^\sigma(S_k)\) / \(E_U\) vary under strain \(\varepsilon\).

\[
\sigma_{aa}=-\frac{1}{\Omega}\frac{\partial E_L}{\partial\varepsilon_{aa}}
\quad\text{(CP2K STRESS\_TENSOR sign; Ha/Bohr³)}
\]

Closed-shell multi-k RKS: spin density \(P_k/2\), \(E_U=2 E_U^\sigma\).

## Acceptance

1. **No throw** on `gpw_stress(atoms, cfg, h)` when `cfg` has non-Γ `k_mesh`/`k_points`
   **and** non-empty `plus_u` (local_pp_on_grid periodic foothold).
2. **Dual-gate unit:** high-spin H₂ UKS + MP \(1\times1\times2\) (or \(2\times2\times2\))
   + PLUS_U \(L=0\), \(U=4\) eV: max\|σ_an − σ_fd\| on diag ≤ \(2\times10^{-3}\) Ha/Bohr³
   (unit smoke; product regression tighter).
3. **Regression** `regressions/gpw/gpw-dftu-kpoint-stress/`:
   - Dual-gate analytic vs multi-SCF FD (decisive).
   - Off-diagonal product ≈ 0.
   - Soft absolute energy pin (ATV first-green or CP2K when physical).
   - e_dftu non-trivial.
4. Offline CP2K twin under `tools/cp2k_oracle/inputs/gpw-dftu-kpoint-stress/`
   (NUMERICAL stress + PLUS_U + KPOINTS). Absolute stress digits only if CP2K
   Mulliken is physical; otherwise dual-gate + note (same pathology class as
   `gpw-dftu-kpoints`).
5. Public API: `GpwDft` members only — no `*Config`, no `make_*`/`build_*`.
6. Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) or `perf_exempt` + reason.

## Out of scope

- Molecular Gaussian DFT+U (PySCF loop).
- OT + multi-k + U.
- Shear (xy/xz/yz) stress (orthorhombic diagonal only).
- CUDA port.
- Solid NiO AFM multi-k + U product twin (follow-on if needed).

## Impl plan

1. Remove `plus_u` throw in `energy_fixed_P_k` (`src/cpu/gpw/internal/gpw_stress.cpp`).
2. After strained \(S_k\) build, evaluate BZ-averaged Dudarev \(E_U\) via
   `build_dftu_potential_k` (energy only; discard \(V\)) on spin densities —
   UKS: \(P_k^\alpha,P_k^\beta\); RKS: \(P_k/2\) with factor 2.
3. Add \(e_{\mathrm{dftu}}\) into the Lagrangian return.
4. Flip unit test `DftUKpointStillThrows` → dual-gate match.
5. Regression + REFERENCE + CMake `atv_add_regression`.
6. Scoped build/test; commit.

## Deps

- `gpw-dftu-kpoints` ✓ (multi-k + U energy)
- `gpw-dftu-stress` ✓ (Γ + U stress Lagrangian pattern)
- `gpw-kpoint-uks-stress` ✓ (multi-k UKS stress machinery)
