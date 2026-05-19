---
title: libxc XC Fix 05 — End-to-end UKS Validation
slug: libxc-xc-fix-05-uks-validation
status: approved
created: 2026-05-19
---

# libxc XC Fix 05 — End-to-end UKS Validation

## Summary

Closing sub-spec of the chain. With native XC matched to libxc point-by-point (spec 03), the Lebedev + Mura-Knowles grid replacing the product grid plus the Bragg-Slater table extended to Z=1..36 (spec 04), and the libxc CMake option + wrapper available (specs 01-02), prove the fix end-to-end: Atomiverse UKS-LSDA total energies agree with PySCF reference values within the per-system tolerance budget on H / He / Li atoms and on H₂ at equilibrium plus a 10-point PES scan, the native and libxc XC paths agree to 1e-12 Ha on each of those systems, and the full ctest suite passes in both `-DATV_USE_LIBXC=ON` and `-DATV_USE_LIBXC=OFF`. `python/examples/03_uks.py` is flipped from its current "KNOWN OFFSET" prose to a hard PASS/FAIL gate at the 1% relative-energy threshold.

## Domain basis

UKS-LSDA total energy: `E = E_1e + E_J + E_xc + V_nn`, with `E_xc = ∫ ρ(r) ε_xc(ρ_α(r), ρ_β(r)) dr` evaluated by quadrature on the molecular grid. The validation reference is PySCF (Sun et al., *J. Chem. Phys.* 153, 024109 (2020), DOI: 10.1063/5.0006074) running `dft.UKS` with `xc="lda,vwn"` → Slater + VWN5 (libxc id 7; Vosko, Wilk, Nusair, *Can. J. Phys.* 58, 1200 (1980)). PySCF grid: `atom_grid={s: (50, 110) for s in ...}` with `prune=None` to match Atomiverse's `(n_rad=50, n_ang=110)`.

Tolerance budget (per-system `|E_atv − E_PySCF|`):

| System              | Target  | Rationale                                    |
|---------------------|---------|----------------------------------------------|
| H atom (doublet)    | 1e-6 Ha | one-electron; grid-limited                    |
| He atom (singlet)   | 1e-5 Ha | closed-shell two-electron LSDA benchmark      |
| Li atom (doublet)   | 1e-4 Ha | multi-shell spin-polarized; larger partition  |
| H₂ at 0.74 Å        | 1e-5 Ha | equilibrium; two-atom Becke partition         |
| H₂ PES 10 points    | 1e-4 Ha | scan R ∈ {0.5, 0.6, 0.74, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5} Å |

Native vs libxc parity: with both XC evaluators driven from the same converged density on the same grid, `|E_xc^native − E_xc^libxc| < 1e-12 Ha` per system.

Units: Hartree throughout the SCF.

## Design

No new public symbols. Exercises existing infrastructure:

- `cpu::scf::uks_solve(ints, n_alpha_occ, n_beta_occ, XCType::LSDA, n_rad, n_ang, tol, max_iter)`.
- `cpu::scf::evaluate_xc(...)` — fixed by spec 03.
- `cpu::scf::evaluate_xc_libxc(...)` (under `ATV_HAS_LIBXC`) — added by spec 02.
- `cpu::scf::Grid` — fixed by spec 04.

New test fixtures encapsulate the five validation systems. PySCF reference energies are baked in as hardcoded `constexpr` constants alongside the geometry that produced them; a header comment documents the PySCF version, `mol.basis`, grid spec, `conv_tol`, and the one-line `dft.UKS` invocation so the numbers are reproducible.

Native↔libxc parity is asserted by re-evaluating `E_xc` on the converged UKS density with both XC paths in the same test binary. The block is `#ifdef ATV_HAS_LIBXC`-gated so the OFF build skips cleanly.

`03_uks.py` is restructured: the "KNOWN OFFSET" paragraph is replaced with a "validation gate" paragraph documenting the tolerance budget; the trailing summary block computes `atom_max_abs_dE`, `h2_max_abs_dE`, the 1% relative threshold check, and `sys.exit(0 if passed else 1)`.

CI side: `CMakeLists.txt` already exposes `ATV_USE_LIBXC` (spec 01). No new CMake changes — the "build matrix" regression criterion is satisfied by running ctest twice with the option flipped.

## Files to create or modify

- `tests/cpu/scf/test_uks_solve.cpp` — replace loose `EXPECT_LT(E, -2.70)` bounds with hard PySCF-reference comparisons; add the H₂ PES scan; add the native-vs-libxc parity test under `ATV_HAS_LIBXC`.
- `python/examples/03_uks.py` — replace "KNOWN OFFSET" docstring + final-line prose with a PASS/FAIL gate; add `sys.exit` with non-zero status on failure.

## Tasks

- [ ] Write failing UKS-vs-PySCF reference tests for H / He / Li / H₂@0.74 / H₂ 10-point PES scan with hardcoded reference energies and per-system tolerances (`tests/cpu/scf/test_uks_solve.cpp`)
- [ ] Write failing native-vs-libxc `E_xc` parity test on converged UKS densities for all five systems, gated on `ATV_HAS_LIBXC` (`tests/cpu/scf/test_uks_solve.cpp`)
- [ ] Tighten `kEnergyTol` and remove the loose `EXPECT_LT/GT` energy-bound assertions superseded by the new reference comparisons
- [ ] Flip `03_uks.py` from "KNOWN OFFSET" prose to a PASS/FAIL gate with per-system Ha tolerances on atoms and 1% relative threshold on H₂ converged points; exit non-zero on failure
- [ ] Run full ctest suite with `-DATV_USE_LIBXC=ON` and verify all tests pass
- [ ] Run full ctest suite with `-DATV_USE_LIBXC=OFF` and verify all tests pass with the libxc-parity block skipped

## Testing strategy

- **Per-system reference energy** (`test_uks_solve`): for each of {H, He, Li, H₂@0.74, H₂ scan}, assert `|E_atv − E_PySCF_const|` is below the system's target. The 10-point H₂ scan is a single parameterized test that fails on the first violating point.
- **Native↔libxc `E_xc` parity** (`test_uks_solve`, `#ifdef ATV_HAS_LIBXC`): converge UKS once per system, recompute `E_xc` from the same `P_α`, `P_β` on the same grid via both `evaluate_xc` and `evaluate_xc_libxc`, assert `|ΔE_xc| < 1e-12 Ha`.
- **Cross-method sanity** (kept from current file): UKS ≠ UHF for He; LSDA ≠ PBE for He; density-matrix shapes; integrated electron counts; α/β symmetry on closed-shell He.
- **Build-matrix regression**: `ctest --test-dir build --output-on-failure` passes under both `-DATV_USE_LIBXC=ON` (parity block exercised) and `-DATV_USE_LIBXC=OFF` (parity block compiles out).
- **Python example gate**: `python python/examples/03_uks.py` exits status 0 when the tolerance budget is met and non-zero otherwise.

## Out of scope

- meta-GGA / hybrid / range-separated — only LSDA here.
- AIMD energy-conservation under SCFEvaluator — separate downstream concern.
- CUDA SCF — legacy, tracked separately.
- Tighter than 1e-12 native↔libxc parity — 1e-12 is the iron-rule bar.
- Pruned / SG-1 grid scoring — fixed `(50, 110)` per spec 04.
- Replacing the bake-in `constexpr` reference with a JSON / Zarr file.
