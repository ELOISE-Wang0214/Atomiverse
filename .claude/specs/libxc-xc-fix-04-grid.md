---
title: Grid fix — tabulated Lebedev, Mura-Knowles radial, extended Bragg-Slater
slug: libxc-xc-fix-04-grid
status: approved
created: 2026-05-19
---

# Grid fix — tabulated Lebedev, Mura-Knowles radial, extended Bragg-Slater

## Summary

Replace the three abandoned approximations in `src/cpu/scf/grid.cpp` that the original `libxc-xc-fix.md` flagged as XC-error sources. The angular grid switches from a product (θ, φ) sweep — which is not a Lebedev rule at all — to tabulated Lebedev-Laikov octahedral points of orders 6, 50, 110, 194, 302, 590. The radial grid switches from the Becke (1-cos) mapping with a single ad-hoc radius parameter to the Mura-Knowles transform, the simpler of the two physically-correct choices and matching PySCF's default for non-pruned grids. The Bragg-Slater atomic radii table grows from the current 3-element stub to a frozen Z=1..36 lookup with hard rejection beyond Kr, satisfying the project's no-fallback policy. Outcome: XC integration grids that pass standard exactness tests independently of the later libxc parity check.

## Domain basis

**Angular quadrature**: V. I. Lebedev and D. N. Laikov, "A quadrature formula for the sphere of the 131st algebraic order of accuracy", *Doklady Mathematics* 59(3) 477-481 (1999). An N-point Lebedev rule of stated degree L exactly integrates every spherical harmonic Y_l^m with l ≤ L; the orthogonality identity `Σ_i w_i |Y_l^m(Ω_i)|² = 1` holds for l ≤ L when weights are normalised to the unit sphere.

**Radial quadrature**: M. E. Mura and P. J. Knowles, "Improved radial grids for quadrature in molecular density-functional calculations", *J. Chem. Phys.* 104, 9848 (1996), DOI 10.1063/1.471749. The transform `x ∈ (0,1) → r = -R_α ln(1 - x³)` carries the unit interval onto (0, ∞); using the N-point uniform open Newton-Cotes rule on x gives integrals of the form `∫₀^∞ f(r) r² dr = Σ_i w_i f(r_i)` with `w_i = 3 R_α x_i² r_i² / (N (1 - x_i³))`. For a Gaussian moment `∫₀^∞ r^(n+2) exp(-α r²) dr` the analytic value is `Γ((n+3)/2) / (2 α^((n+3)/2))`. Treutler-Ahlrichs (*J. Chem. Phys.* 102, 346 (1995)) is the alternative non-pruned choice and is explicitly deferred.

**Atomic radii**: J. C. Slater, *J. Chem. Phys.* 41, 3199 (1964), DOI 10.1063/1.1725697. Noble-gas entries (He, Ne, Ar, Kr) use Clementi-Raimondi-Reinhardt calculated radii, *J. Chem. Phys.* 47, 1300 (1967), DOI 10.1063/1.1712084. Tabulated for Z = 1..36 in Angstrom; consumed internally in Bohr via `kAngstromToBohr` from `core/units.h`.

Units: radii Angstrom on input, Bohr on output; Lebedev weights dimensionless and sum to 1; Mura-Knowles weights carry an r² Jacobian.

## Design

Three changes localised to `src/cpu/scf/grid.cpp`, plus growth of the already-committed `lebedev_data.h`:

- **Angular grid**: `make_angular_grid(order)` becomes a dispatch over the six tabulated Lebedev orders, returning `std::span<const LebedevPoint>` views into static arrays in `lebedev_data.h`. The current product-grid code path is deleted. Unsupported orders throw `std::invalid_argument`.
- **Radial grid**: `make_radial_grid(N, R_α)` is rewritten in terms of the Mura-Knowles transform. Quadrature nodes are uniform-open in `x_i = i / (N+1)`, i = 1..N; `r_i` and `w_i` derived from the closed-form transform.
- **Atomic radii**: `bragg_slater_radius(Z)` is replaced by a `constexpr std::array<double, 37>` indexed by Z (entry 0 unused) holding the 36 Angstrom values. Any Z outside 1..36 throws `std::invalid_argument` with a message naming Z; a deliberate hard wall — adding period-5 entries is a follow-up spec.

`LebedevPoint` (x, y, z, w) and the public surface of `grid.h` (`make_atomic_grid`, `make_molecular_grid`, Becke partition) are unchanged.

`lebedev_data.h` currently ships only orders 6 and 110 (per the header comment). Orders 50, 194, 302, 590 must be appended. The points may be regenerated with any one of: Burkardt's `sphere_lebedev_rule.c`, libxc's shipped tables, or PySCF's `grid_basis.c`; the generator script is not committed but its provenance must be noted in the header comment.

Bragg-Slater values for Z=1..36 (Angstrom): H=0.25, He=0.31, Li=1.45, Be=1.05, B=0.85, C=0.70, N=0.65, O=0.60, F=0.50, Ne=0.38, Na=1.80, Mg=1.50, Al=1.25, Si=1.10, P=1.00, S=1.00, Cl=1.00, Ar=0.71, K=2.20, Ca=1.80, Sc=1.60, Ti=1.40, V=1.35, Cr=1.40, Mn=1.40, Fe=1.40, Co=1.35, Ni=1.35, Cu=1.35, Zn=1.35, Ga=1.30, Ge=1.25, As=1.15, Se=1.15, Br=1.15, Kr=0.88.

## Files to create or modify

- `include/atomiverse/cpu/scf/lebedev_data.h` — append point tables for orders 50, 194, 302, 590; keep existing 6 and 110 untouched; update header comment.
- `src/cpu/scf/grid.cpp` — rewrite `make_angular_grid` as six-way Lebedev dispatch; rewrite `make_radial_grid` as Mura-Knowles; replace `bragg_slater_radius` with the Z=1..36 table and hard rejection; remove dead product-grid and Becke-mapping code.
- `tests/cpu/scf/test_grid.cpp` — add three test groups for spherical-harmonic orthogonality, Gaussian radial moments, and Bragg-Slater coverage/strictness.

## Tasks

- [ ] Write failing tests for spherical-harmonic orthogonality, Gaussian radial moments, and Bragg-Slater coverage (`tests/cpu/scf/test_grid.cpp`)
- [ ] Extend `include/atomiverse/cpu/scf/lebedev_data.h` with tabulated points and weights for orders 50, 194, 302, 590
- [ ] Implement Lebedev dispatch in `make_angular_grid` over orders {6, 50, 110, 194, 302, 590} in `src/cpu/scf/grid.cpp`
- [ ] Implement Mura-Knowles radial transform in `make_radial_grid` in `src/cpu/scf/grid.cpp`, removing the Becke-mapping code
- [ ] Implement the Z=1..36 Bragg-Slater table in `bragg_slater_radius` with `std::invalid_argument` outside the range
- [ ] Add Doxygen `@brief` / `@param` (with units) / `@return` to every touched public or file-scope symbol
- [ ] Run full check + test suite (`cmake --build build -j && ctest --test-dir build --output-on-failure`)

## Testing strategy

- **Lebedev exactness (happy path)**: for each of the six tabulated orders, assert `Σ_i w_i |Y_l^m(Ω_i)|² == 1` to within 1e-14 for every (l, m) with l ≤ stated_degree. Real spherical harmonics computed in-test from associated Legendre polynomials.
- **Mura-Knowles radial moments (happy path + domain validation)**: for α in {0.25, 1.0, 4.0} and n in {0, 2, 4, 6}, assert N=96 radial quadrature of `r^(n+2) exp(-α r²)` matches `Γ((n+3)/2) / (2 α^((n+3)/2))` to 1e-10 relative.
- **Bragg-Slater coverage (boundaries)**: every Z in 1..36 returns the tabulated value (converted to Bohr) bit-exactly; Z = 0, Z = 37, Z = -1 each throw `std::invalid_argument`.
- **Edge case**: `make_angular_grid` with an order not in the six-element set throws `std::invalid_argument`.
- **Smoke**: existing `make_molecular_grid` callers in the SCF test suite continue to link and run; an H2 single-point at the new grid behaves at least as well as before (no XC parity claim — sub-spec 05).

## Out of scope

- Treutler-Ahlrichs radial transform — deferred to a follow-up.
- Bragg-Slater entries for Z > 36 (period 5+); current spec rejects them by design.
- Grid pruning (Treutler / SG-1 schemes); all atoms get the full angular order at every radial shell.
- libxc-vs-native XC parity — sub-spec 05 dependency.
- Python generator script for Lebedev tables — referenced by provenance comment but not committed.
