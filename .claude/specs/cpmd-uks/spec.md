# Spec — `cpmd-uks`

**Slug:** `cpmd-uks`  
**Method:** Car–Parrinello MD on Γ-UKS GPW (fictitious electronic DOF + RATTLE)  
**Domain:** dynamics / GPW  
**Branch:** `wf/cpmd-uks`  
**Status:** in progress (CP2K capability loop wave 7)

## Intent

Extend product CPMD (`Cpmd` + `GPWEvaluator::cpmd`) from closed-shell RKS Γ
to **collinear Γ-UKS**: independent occupied manifolds \(C_\alpha\), \(C_\beta\)
with occupancy-1 S-tangent forces, fixed-density dual KS evaluation, and
analytic nuclear forces at the current spin densities.

**Not** second-gen CP UKS (`second-gen-cpmd-uks`), **not** XL-BOMD UKS
(already product: `xl-bomd-uks`), **not** BOMD full SCF each step.

## Background (on tree)

| Piece | State |
| ----- | ----- |
| `Cpmd` + electronic RATTLE | ✓ RKS (`cpmd-water`, `cpmd-rattle-water`) |
| `GPWEvaluator::cpmd` | rejects `n_alpha`/`n_beta` at `build()` |
| `GpwDft` fixed-density (`max_iter==0`) | RKS only |
| `gpw_forces_at_density` | RKS only |
| Γ-UKS energy/forces | ✓ (`gpw-uks-forces`, `gpw-uks-o2-gamma`) |
| XL-BOMD dual α/β dens. | ✓ (`xl-bomd-uks`) |

## Design

1. **`Cpmd::advance(H, S, dt_fs, occ=2.0)`** — pass orbital occupation into
   `ot_gradient` (RKS default 2; UKS spin channels 1). Dual instances share μ.
2. **`GpwDft` fixed-density UKS** — when `max_iter==0` and UKS: require
   `guess_P` + `guess_P_beta`; one collocate→Poisson→spin-XC→integrate;
   populate `H_ks` (α), **`H_ks_beta`** (β), `P`/`P_beta`, `S`, `e_total`,
   `iters=0`, `converged=true`.
3. **`gpw_forces_at_density` UKS overload** — \(W_\sigma = C_\sigma(C_\sigma^\top H_\sigma C_\sigma)C_\sigma^\top\),
   \(W=W_\alpha+W_\beta\); assemble with dual densities.
4. **`GPWEvaluator`** — when UKS + `cpmd`:
   - seed both `cpmd_` (α) and `cpmd_beta_` (β) after cold SCF;
   - production: \(P_\sigma = C_\sigma C_\sigma^\top\), fixed-density eval,
     S-align both, forces at density, advance both with `occ=1`;
   - \(T_e = T_e^\alpha + T_e^\beta\).
5. **Guards** — still refuse OT, XL, second-gen CP combo; k-points remain
   out of scope (Γ only).

## Oracle / pin

- Twin geometry: O₂ triplet (same family as `xl-bomd-uks` / `gpw-uks-forces`).
- Offline CP2K input: `tools/cp2k_oracle/inputs/cpmd_uks_o2_nve.inp`
  (QS PADE/LSD, DZVP-MOLOPT-GTH, GTH-PADE, MT 12 Å, MULTIPLICITY 3).
- **ctest pin** is literature **\(E_{\mathrm{cons}}\)** drift envelope +
  structural / ⟨S²⟩ gates + production `iters==0` — same inventory class as
  `cpmd-rattle-water` / `xl-bomd-uks` (no live CP2K in ctest).
- **perf_exempt** — NVE drift envelope, no fair-wall twin on worker host.

## Acceptance

See `acceptance.md`.

## Out of scope

- Second-gen CP on UKS; electronic thermostat; multi-k CPMD; OT-CP hybrid.
- Public `*Config` / `make_*` / `build_*` factories.
