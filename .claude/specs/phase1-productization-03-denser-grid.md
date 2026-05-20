---
title: Denser, properly-pruned DFT integration grid (SG-1 + match-PySCF)
status: approved
created: 2026-05-20
---

# Denser, properly-pruned DFT integration grid (SG-1 + match-PySCF)

## Summary
The CPU SCF DFT integration grid is currently an unpruned "SG-0-ish" grid
(`n_radial=50`, `n_angular=110` Lebedev, uniform angular order at every radial
shell). This coarseness leaves a ~5 mHa ATV-vs-PySCF energy residue on def2-SVP,
which forces `bench_qc_suite` to hold its energy acceptance threshold at a loose
1e-2 Ha. This spec replaces the grid with two named, properly-constructed grid
levels selectable through `SCFConfig` / `default_config`: (1) the **SG-1**
standard grid (Gill-Johnson-Pople 1993 — Euler-Maclaurin radial, 5-region
Lebedev pruning, Becke size-adjusted weights), which brings the residue down to a
few hundred microHartree; and (2) a **match-PySCF** high-accuracy grid
(Treutler-Ahlrichs radial + PySCF level-3 Lebedev orders + nwchem-style angular
pruning), which reaches the sub-microHartree regime. The new default and the
benchmark threshold are tightened honestly to what each mode actually achieves —
SG-1 is held to ~5e-4 Ha, match-PySCF to ~1e-6 Ha, with no aspirational
"1 microHartree with SG-1" claim.

## Domain basis
All schemes and constants below were verified this session by the `scientist`
agent and are reproduced verbatim.

**SG-1 standard grid** — Gill, Johnson, Pople, *Chem. Phys. Lett.* **209** (1993)
506, DOI 10.1016/0009-2614(93)80125-9. Parent grid EML-(50,194): Euler-Maclaurin
radial quadrature, 50 radial points, uniform for H-Ar. Lebedev angular grid
pruned into 5 radial regions with orders {6, 38, 86, 194, 86} from innermost to
outermost. Region boundaries are alpha_1*R .. alpha_4*R with atomic radius R; the
pruning alpha-quadruples are H,He = (0.25, 0.5, 1.0, 4.5); Li-Ne =
(0.1667, 0.5, 0.9, 3.5); Na-Ar = (0.1, 0.4, 0.8, 2.5). Atom-size adjustment uses
R = 1/2 x Bragg-Slater radius for every atom **except hydrogen** (no 1/2 factor
for H).

**Match-PySCF high-accuracy grid** — PySCF default `level=3`. Radial: number of
points (n_rad, n_ang) = (50, 590) for H/He, (75, 590) for Li-Ne, (80, 770) for
Z > 10. Radial scheme: Treutler-Ahlrichs M4 mapping (O. Treutler & R. Ahlrichs,
*J. Chem. Phys.* **102** (1995) 346, DOI 10.1063/1.469408). Angular pruning:
nwchem-style. Atomic partitioning: original Becke ("original_becke").

**SG-2 / SG-3** — Dasgupta & Herbert, *J. Comput. Chem.* **38** (2017) 869,
DOI 10.1002/jcc.24761. SG-2 parent (75, 302), SG-3 parent (99, 590),
double-exponential radial scheme; SG-2 targets meta-GGAs. Cited as future work
only — not implemented here.

**Becke fuzzy-cell partitioning** — Becke, *J. Chem. Phys.* **88** (1988) 2547,
DOI 10.1063/1.454033. Cell weight w_n(r) = P_n(r) / sum_m P_m(r), with
P_n = product over m != n of s(mu_nm). The confocal-elliptical coordinate is
mu_ij = (r_i - r_j) / R_ij. The switching function is
s(mu) = 1/2 [1 - f(f(f(mu)))] with f(mu) = (3/2) mu - (1/2) mu^3 iterated k = 3
times. Atomic-size adjustment is nu_ij = mu_ij + a_ij (1 - mu_ij^2), with a_ij
derived from chi = R_i / R_j via u = (chi - 1)/(chi + 1), a = u/(u^2 - 1),
clamped to |a| <= 1/2.

**Accuracy bounds (honesty constraint).** SG-1's design accuracy target is
~0.2 kcal/mol ~= 320 microHartree; it is only modestly denser than the current
SG-0-ish grid. Sub-microHartree agreement vs PySCF requires the match-PySCF grid
that reproduces PySCF's actual level-3 default — not SG-1. Acceptance thresholds
below are set per mode accordingly.

## Design
A new enum `GridLevel { SG0, SG1, MatchPySCF }` is added to `grid.h` and carried
by `SCFConfig` as a field `grid_level` (default `MatchPySCF`). `SCFConfig`
retains `n_radial` / `n_angular` only for explicit-override / legacy callers; a
named level overrides them when not `SG0`.

`grid.cpp` gains:

- `radial_grid_euler_maclaurin(int n, double R)` — Euler-Maclaurin radial
  quadrature for the SG-1 parent grid, R = atom-size-adjusted radius.
- A pruning dispatch: a free function `lebedev_order_for_shell(GridLevel level,
  int Z, double r, double R, int radial_index, int n_radial)` returning the
  Lebedev order to use at a given radial shell. For `SG1` it applies the
  5-region alpha-quadruple table; for `MatchPySCF` it applies nwchem-style
  pruning bounded by the per-element level-3 maximum order (590 or 770).
- New Lebedev tables `grid_38`, `grid_86`, `grid_770` in `lebedev_data.h`
  (orders 38, 86, 770 are absent today; SG-1 needs 38 and 86, match-PySCF
  needs 770 for Z > 10). The angular-order dispatch in `build_molecular_grid`
  is extended to recognize the new orders.
- `build_molecular_grid` is refactored so the angular order is selected
  **per radial shell** via `lebedev_order_for_shell`, instead of one uniform
  order for the whole atom. The SG-0 path (uniform order, no pruning) is kept
  as the `SG0` branch for backward compatibility and to keep the existing
  legacy `default_config` regression baselines stable.

The existing Becke partitioning (`becke_partition_weights`) is audited against
the Domain basis formulas (f iterated 3×, size-adjustment chi→u→a, |a|<=1/2
clamp); it is corrected only if a discrepancy is found, otherwise left in place
with a note in the docstring confirming "original_becke" parity.

`default_config` switches its grid to `GridLevel::MatchPySCF` and its
"SG-0-ish" doc comment is rewritten. The three internal call sites that read
`config.n_radial` / `config.n_angular` (`gradient.cpp`,
`internal/hamiltonian_builder.cpp`, `uks_solve.cpp`) are routed through the
level-aware `build_molecular_grid`. `bench_qc_suite` tightens its
`threshold_ha` from 1e-2 to the residue actually achieved by the new default.

Ownership / lifecycle: grids remain transient value objects (`MolecularGrid`)
built once per geometry inside `HamiltonianBuilder` / gradient assembly; no new
long-lived state. The Lebedev tables are `constexpr` static data, consistent
with the existing `lebedev_data.h`.

## Files to create or modify
- `include/atomiverse/cpu/scf/grid.h` — add `GridLevel` enum,
  `radial_grid_euler_maclaurin` declaration, `lebedev_order_for_shell`
  declaration, level-aware `build_molecular_grid` overload/parameter.
- `src/cpu/scf/grid.cpp` — implement Euler-Maclaurin radial, SG-1 5-region
  pruning, nwchem-style pruning, per-shell angular dispatch; audit/fix
  `becke_partition_weights`.
- `include/atomiverse/cpu/scf/lebedev_data.h` — add `grid_38`, `grid_86`,
  `grid_770` tabulated Lebedev grids.
- `include/atomiverse/cpu/scf/scf_solve.h` — add `grid_level` to `SCFConfig`;
  update `default_config` and rewrite its "SG-0-ish" doc comment.
- `include/atomiverse/cpu/scf/xc_grid_engine.h` — update only if the
  `MolecularGrid` consumer signature changes (expected: no change; touched
  only if the `build_molecular_grid` signature shift requires it).
- `src/cpu/scf/internal/hamiltonian_builder.cpp` — route grid build through
  the level-aware path.
- `src/cpu/scf/gradient.cpp` — route grid build through the level-aware path.
- `src/cpu/scf/uks_solve.cpp` — route grid build through the level-aware path.
- `tests/cpu/scf/test_grid.cpp` — extend with electron-count, Lebedev
  exactness, and SG-1/match-PySCF energy-residue tests.
- `benchmarks/qc/bench_qc_suite.cpp` — tighten `threshold_ha`.

## Tasks
- [ ] Write failing tests for Lebedev angular exactness of new orders 38, 86, 770 (tests/cpu/scf/test_grid.cpp)
- [ ] Add `grid_38`, `grid_86`, `grid_770` Lebedev tables to lebedev_data.h and extend the angular-order dispatch in build_molecular_grid
- [ ] Write failing tests for SG-1 and match-PySCF electron-count recovery on water/def2-SVP (tests/cpu/scf/test_grid.cpp)
- [ ] Implement `radial_grid_euler_maclaurin` and `lebedev_order_for_shell` (SG-1 5-region + nwchem pruning) in src/cpu/scf/grid.cpp
- [ ] Add `GridLevel` enum to grid.h and refactor `build_molecular_grid` to select angular order per radial shell
- [ ] Audit and, if needed, fix `becke_partition_weights` against the Becke 1988 size-adjusted formula in src/cpu/scf/grid.cpp
- [ ] Add `grid_level` to `SCFConfig`, route gradient.cpp / hamiltonian_builder.cpp / uks_solve.cpp through the level-aware grid build
- [ ] Update `default_config` to `GridLevel::MatchPySCF` and rewrite the "SG-0-ish" doc comment in scf_solve.h
- [ ] Write failing tests for SG-1 (<= 5e-4 Ha) and match-PySCF (<= 1e-6 Ha) energy residue on H2O/NH3/CH4 def2-SVP, then verify
- [ ] Tighten `bench_qc_suite` `threshold_ha` from 1e-2 Ha to the achieved residue and add doxygen units to all new public symbols
- [ ] Run full check + test suite

## Testing strategy
- **Happy path.** SG-1 and match-PySCF grids build without error on
  H2O/NH3/CH4 at def2-SVP; `build_molecular_grid` returns a non-empty
  `MolecularGrid` whose `points` count matches the expected per-shell pruned
  total.
- **Edge cases.** Single-atom molecule (Becke weights = 1, no pruning regression);
  hydrogen atom-size adjustment uses the no-1/2-factor radius; pure-HF path
  (`xc_type` HF) is unaffected because it builds no grid; legacy `SG0` level
  still reproduces the current (50,110) behavior bit-for-bit.
- **Domain validation** (`science.required`):
  - Lebedev exactness self-test — integrating spherical harmonics up to each
    grid's design algebraic degree returns the analytic value to machine
    precision (orders 6, 38, 50, 86, 110, 194, 302, 590, 770).
  - Electron-count integral — sum of grid weights times electron density
    recovers N_electrons on water/def2-SVP to <= 1e-4 for both SG-1 and
    match-PySCF.
  - ATV-vs-PySCF single-point energy residue on the H2O/NH3/CH4 triplet:
    SG-1 <= 5e-4 Ha, match-PySCF <= 1e-6 Ha.
  - Becke partition unity — per-grid-point cell weights sum to 1 within 1e-12.

## Out of scope
- Full SG-2 / SG-3 grid implementation (double-exponential radial, meta-GGA
  targeting). Cited as future work; only SG-1 + match-PySCF land here.
- The CUDA grid path — this spec touches CPU SCF only.
- Meta-GGA functional support (SG-2's design driver).
- Changing the radial point counts of the legacy `SG0` default beyond what is
  needed to preserve its current behavior.
