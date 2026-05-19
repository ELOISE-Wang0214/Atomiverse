---
slug: libxc-xc-fix-04-grid
criteria:
  - id: ac-001
    summary: Lebedev grids integrate spherical harmonics exactly to stated degree
    type: scientific
    pass_when: |
      For each Lebedev order in {6, 50, 110, 194, 302, 590} with stated
      algebraic degree L in {3, 11, 17, 23, 29, 41}, the test
      `test_grid.LebedevSphericalHarmonicOrthogonality` reports
      |Σ_i w_i |Y_l^m(Ω_i)|² - 1| < 1e-14 for every (l, m) with
      0 ≤ l ≤ L and |m| ≤ l, and the test passes under
      `ctest --test-dir build --output-on-failure -R test_grid`.
    status: pending
  - id: ac-002
    summary: Mura-Knowles radial grid reproduces analytic Gaussian moments
    type: scientific
    pass_when: |
      With N = 96 radial points, for every (α, n) in
      {0.25, 1.0, 4.0} × {0, 2, 4, 6}, the test
      `test_grid.MuraKnowlesGaussianMoments` reports
      |Q_num - Γ((n+3)/2)/(2 α^((n+3)/2))| / |Q_analytic| < 1e-10,
      where Q_num is the quadrature `Σ_i w_i r_i^n exp(-α r_i²)`, and
      passes under `ctest --test-dir build --output-on-failure -R test_grid`.
    status: pending
  - id: ac-003
    summary: Bragg-Slater table covers Z=1..36 exactly and rejects Z=37
    type: code
    pass_when: |
      For every Z in 1..36, `bragg_slater_radius(Z)` returns the tabulated
      Angstrom value (H=0.25 ... Kr=0.88) multiplied by `kAngstromToBohr`
      bit-exactly; calls with Z in {-1, 0, 37, 100} each throw
      `std::invalid_argument` whose `what()` contains the offending Z.
      `test_grid.BraggSlaterCoverage` asserts both halves and passes under
      `ctest --test-dir build --output-on-failure -R test_grid`.
    status: pending
---

# Acceptance criteria

## ac-001 — Lebedev exactness (scientific)
The shipped Lebedev tables are the entire foundation of the angular grid; if they do not integrate spherical harmonics exactly to their stated degree, every downstream XC integral is biased. The 1e-14 tolerance is the realistic double-precision floor for sums of ~600 unit-magnitude squared values.

## ac-002 — Mura-Knowles exactness (scientific)
Gaussian radial moments are the canonical analytic check for radial quadrature in molecular DFT: every basis-function pair product collapses to a finite linear combination of such moments. The 1e-10 relative bar with N = 96 is conservative for Mura-Knowles on the chosen α range; failure indicates either the transform Jacobian or the node placement is wrong.

## ac-003 — Bragg-Slater strictness (code)
This criterion enforces the project's no-fallback rule at the table boundary. Bit-exact return values prevent silent precision regressions during future refactors; the throw-on-Z=37 half ensures an Rb run fails on the first call rather than producing a wrong density.
