---
title: libxc XC Fix 05 — End-to-end UKS Validation
slug: libxc-xc-fix-05-uks-validation
status: approved
created: 2026-05-19
revised: 2026-05-19
depends_on:
  - libxc-xc-fix-03-native-parity
  - libxc-xc-fix-04-grid
---

# libxc XC Fix 05 — End-to-end UKS Validation

## Summary

Closing sub-spec of the chain. With native XC matched to libxc to 1e-12 at the primitive level (spec 03) and the Lebedev + Mura-Knowles grid plus the extended Bragg-Slater table (spec 04), prove the fix end-to-end against PySCF reference numbers. Atomiverse UKS-LSDA total energies must agree with hardcoded PySCF `constexpr` references within a per-system Hartree tolerance budget on H / He / Li atoms, H₂ at 0.74 Å, and a 10-point H₂ PES scan, with `python/examples/03_uks.py` flipped from "KNOWN OFFSET" prose to a PASS/FAIL gate at the 1% relative-energy threshold (non-zero exit on failure).

**libxc is intentionally not referenced anywhere in this sub-spec.** That was the whole point of the chain: with primitive parity locked at 1e-12 in spec 03, every integrated quantity is mathematically equal to a libxc-driven integrated quantity to within numerical noise; a redundant runtime native↔libxc agreement criterion would falsely imply libxc remains in the production picture. libxc is test-only scaffolding for primitive-level parity (spec 03), never a runtime dependency. End-to-end validation is therefore against PySCF — the upstream that drives our scientific bar.

## Domain basis

UKS-LSDA total energy: `E = E_1e + E_J + E_xc + V_nn`, with
`E_xc = ∫ ρ(r) ε_xc(ρ_α(r), ρ_β(r)) dr` evaluated by quadrature on the molecular grid. Validation reference: PySCF (Sun et al., *J. Chem. Phys.* 153, 024109 (2020), DOI 10.1063/5.0006074) running `dft.UKS` with `xc="lda,vwn"` → Slater + VWN5 (Vosko, Wilk, Nusair 1980). PySCF grid: `mf.grids.atom_grid = {s: (50, 110) for s in ...}` with `prune=None`. `mf.conv_tol = 1e-12`.

Tolerance budget (per-system `|E_atv − E_PySCF|`):

| System              | Target  | Rationale                                                    |
|---------------------|---------|--------------------------------------------------------------|
| H atom (doublet)    | 1e-6 Ha | one-electron; grid-limited                                    |
| He atom (singlet)   | 1e-5 Ha | closed-shell two-electron LSDA benchmark                      |
| Li atom (doublet)   | 1e-4 Ha | multi-shell spin-polarized; larger Becke partition            |
| H₂ at 0.74 Å        | 1e-5 Ha | equilibrium; two-atom Becke partition                          |
| H₂ PES 10 points    | 1e-4 Ha | scan R ∈ {0.5, 0.6, 0.74, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5} Å |

**Why no native↔libxc runtime criterion**: spec 03 locks `ε_xc^native − ε_xc^libxc = 0` at the pointwise primitive level to 1e-12 over the full physically reachable (ρ_α, ρ_β) domain. Because `E_xc`, `V_xc`, the SCF iteration, the converged density, and `E_tot` are all deterministic functions of the primitive composed of finite arithmetic on the same grid + same integrals, primitive parity implies end-to-end integrated-quantity parity to within accumulated FP noise. A separate runtime "native vs libxc on converged density" criterion would be redundant and would falsely imply libxc is part of the production picture.

Units: Hartree throughout.

## Design

No new public symbols. Exercises existing infrastructure:

- `cpu::scf::uks_solve(ints, n_alpha_occ, n_beta_occ, XCType::LSDA, n_rad, n_ang, tol, max_iter)`.
- `cpu::scf::evaluate_xc(...)` — fixed by spec 03 to primitive-parity with libxc.
- `cpu::scf::Grid` — fixed by spec 04.

New test fixtures encapsulate the five validation systems. PySCF reference energies are baked in as hardcoded `constexpr double` constants in `tests/cpu/scf/test_uks_solve.cpp`, alongside the geometry that produced them. The test-file header comment documents the PySCF version, `mol.basis="sto-3g"`, the `dft.UKS(xc="lda,vwn")` invocation, grid spec `atom_grid=(50, 110), prune=None`, and `conv_tol=1e-12`.

The 10-point H₂ PES scan uses a value-parameterized GTest sweep (`TEST_P` over the ten R values) so a single violating R is reported by name.

`03_uks.py` is restructured: the "KNOWN OFFSET" paragraph becomes a "validation gate" paragraph documenting the tolerance budget; the trailing summary block computes `atom_max_abs_dE`, `h2_max_abs_dE`, the 1% relative threshold check on converged points, and `sys.exit(0 if passed else 1)`.

CI side: ctest runs once in its default configuration. `ATV_BUILD_XC_ORACLE` (test-only scaffolding for spec 03) is orthogonal to this sub-spec.

## Files to create or modify

- `tests/cpu/scf/test_uks_solve.cpp` — replace loose `EXPECT_LT(E, ...)` bounds with `EXPECT_NEAR(E, kPyscfRefXxx, kTolXxx)` against `constexpr` references; add H₂ PES scan parameterized over ten R values; add test-file header comment documenting how the PySCF references were produced.
- `python/examples/03_uks.py` — replace "KNOWN OFFSET" framing with "validation gate" framing; add `sys.exit(0 if passed else 1)`.

## Tasks

- [ ] Write failing UKS-vs-PySCF reference tests for H / He / Li / H₂@0.74 with hardcoded `constexpr` PySCF energies and per-system Ha tolerances (`tests/cpu/scf/test_uks_solve.cpp`)
- [ ] Write failing parameterized H₂ PES scan test over R ∈ {0.5, 0.6, 0.74, 0.85, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5} Å with `constexpr` PySCF references and 1e-4 Ha per-R tolerance
- [ ] Add test-file header comment documenting PySCF version, `mol.basis`, `dft.UKS(xc="lda,vwn")` invocation, `atom_grid=(50, 110), prune=None`, `conv_tol=1e-12` reproducer
- [ ] Remove the loose `EXPECT_LT / EXPECT_GT` energy-bound assertions superseded by the new hard reference comparisons; tighten `kEnergyTol`
- [ ] Flip `03_uks.py` docstring from "KNOWN OFFSET" framing to "validation gate" framing; add `sys.exit(0 if passed else 1)` driven by `atom_max_rel < 1.0%` AND `h2_conv_max_rel < 1.0%`
- [ ] Run full ctest suite and verify all tests pass

## Testing strategy

- **Per-system reference energy**: for each of {H, He, Li, H₂@0.74}, `|E_atv − E_PySCF_const|` below per-system Ha tolerance.
- **H₂ PES scan**: value-parameterized, fails on first violating R with name.
- **Cross-method sanity** (kept): UKS ≠ UHF for He; LSDA ≠ PBE for He; density-matrix shapes; integrated electron counts; α/β symmetry on closed-shell He.
- **Python example gate**: `python python/examples/03_uks.py` exits 0 when budget met, non-zero otherwise; PASS/FAIL line printed.
- **Full suite regression**: `ctest --output-on-failure` green in default configuration.

## Out of scope

- Native↔libxc end-to-end runtime parity. Removed deliberately: primitive parity in spec 03 already guarantees integrated-quantity parity; including this would imply libxc is a production dependency.
- `ATV_BUILD_XC_ORACLE` ON/OFF build-matrix gating — orthogonal to UKS validation.
- meta-GGA / hybrid / range-separated — LSDA only.
- AIMD energy-conservation — separate concern.
- CUDA SCF — legacy, tracked separately.
- Pruned / SG-1 grid scoring — fixed `(50, 110)` per spec 04.
- Replacing `constexpr` reference with JSON/Zarr file.
