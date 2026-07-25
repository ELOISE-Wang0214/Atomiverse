# Acceptance: `wannier-centers`

## Dual gate

| Gate | Criterion |
| ---- | --------- |
| Correctness | All checks below green |
| Speed | **`perf_exempt`** — no CP2K MLWF twin wall on the pin host |

## Unit (`test_cpu_scf_wannier`)

1. `boys_pair_angle(1, 0) == 0`; `boys_pair_angle(0, 1) == π/8` (within 1e-14).
2. One occupied MO: `compute` returns 1 center, 0 sweeps of meaningful rotation
   (or `converged` with `total_spread == total_spread_initial`), no throw.
3. Synthetic two-orbital Jacobi: after localization `total_spread ≤ total_spread_initial`
   (within 1e-10 relative).
4. Not-ready / open-shell MeanFieldOutputs → `std::invalid_argument` / `runtime_error`.

## Regression (`reg_wannier-centers`)

H₂O experimental-like geometry (same as multipole dipole pin), RHF / STO-3G:

1. SCF converges; `WannierCenters::compute` returns **5** centers and spreads.
2. Every spread > 0 and finite; `total_spread < total_spread_initial` (strict).
3. **Classic Boys water pattern** (sort centers by distance to O):
   - 1 center within **0.15 Å** of O (core).
   - 2 OH-bond centers: each within **0.35 Å** of the corresponding H midpoint
     (or y-sign matched to ±H with |y| > 0.2 Å).
   - 2 lone-pair centers: |x| > 0.15 Å, opposite x signs, both closer to O than
     the bond centers on average (rabbit ears).
4. **Sum of centers** (electronic centroids) invariant under localization to
   1e-8 Å vs pre-localization sum (unit algebra).
5. Optional H₂ pin in the same binary or unit: single center at bond midpoint
   within **1e-4 Å**.
6. Wall ≤ 60 s single rank.

## Docs / hygiene

- Public header has Doxygen `@brief` / `@param` / `@return` with units.
- No public `*Config` / factory.
- Tutorial page exists; roadmap §11 notes the foothold; queue row marked done
  only after merge (worker marks on branch).

## Out of acceptance

- Numerical parity with CP2K `LOCALIZE` (offline oracle only).
- Periodic Wannier / Berry–Wannier dual.
