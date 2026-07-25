# periodic-wannier-si — acceptance

Gate = scoped build/test per loop law:
`atomiverse` + `test_cpu_pw_periodic_wannier` + `reg_periodic-wannier-si`.

## Criteria

1. **Unit — bond midpoints** (`test_cpu_pw_periodic_wannier`):
   `diamond_bond_midpoints(a)` returns the four Cartesian midpoints
   \((\pm a/8,\pm a/8,\pm a/8)\) with even number of minuses (tetrahedral set)
   within 1e-12; pure math, no SCF.
2. **Unit — SMV cubic weight**: `mv_pair_weight_cubic(|b|)` equals
   \(1/(2 b^2)\) (six neighbors, \(w_b=3/(N_b b^2)\)); throws on \(b\le 0\).
3. **Unit — phase → center (analytic)**: for a single k and six ±b neighbors
   with \(M_{nn}=\exp(-i\mathbf b\cdot\mathbf r_0)\),
   `center_from_diagonal_phases` recovers \(\mathbf r_0\) to 1e-9 (Å-scale
   test with synthetic b, w).
4. **Regression `regressions/pw/periodic-wannier-si`** (public API only,
   hard-coded literature geometry pin, ≤ 60 s):
   a. SCF converged; `n_centers == 4`.
   b. **Bond-midpoint pin** (Marzari & Vanderbilt, PRB 56, 12847 (1997) §IV;
      coarse 2×2×2 SMV mesh): after optimal assignment + min-image under the
      fcc cell, every Wannier center lies within **0.90 Å** of one tetrahedral
      Si–Si bond midpoint (bond length 2.35 Å).
   c. **Localization:** \(\Omega_\mathrm{final} < 0.25\,\Omega_\mathrm{init}\).
   d. Finite spreads (no NaN).
5. **Architecture**: no public `*Config`/factory; class + members mirrors
   `BerryPolarization` / `WannierCenters`; no extension of SCF bag debt.
6. **Perf**: fair CP2K wall twin or `perf_exempt` with reason (no CP2K module
   on build host — same precedent as `berry-polarization`).

## Status (2026-07-25)

- [x] 1–3 unit tests green (`test_cpu_pw_periodic_wannier`)
- [x] 4 regression green (`reg_periodic-wannier-si`: dmax 0.83 Å, Ω ratio 0.081)
- [x] 5 architecture: `PeriodicWannier` members, no `*Config`/factory
- [x] 6 perf: `perf_exempt=true` (no CP2K binary); ATV wall ~8 s
