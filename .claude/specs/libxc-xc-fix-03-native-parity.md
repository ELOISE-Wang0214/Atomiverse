---
title: Native XC parity vs. libxc to 1e-12
slug: libxc-xc-fix-03-native-parity
status: approved
created: 2026-05-19
---

# Native XC parity vs. libxc to 1e-12

## Summary

With libxc CMake integration (sub-spec 01) and the `evaluate_xc_libxc()` wrapper (sub-spec 02) in place, drive the native `evaluate_xc()` implementation in `src/cpu/scf/xc_functional.cpp` to numerically agree with libxc to within 1e-12 in `(exc, vrho_alpha, vrho_beta, vsigma_aa, vsigma_ab, vsigma_bb)` at every test point of a curated battery, for both LSDA (Slater + VWN5) and PBE. The work is test-driven: a new GoogleTest binary `test_xc_vs_libxc` is written first, must fail against the current native code, and then the four suspected native bugs (Slater `kCx`, VWN5 parameter signs/decimals, spin-interpolation `f(ζ)` and `f'(ζ)`, vrho convention `∂exc/∂ρ_σ` vs. libxc's `∂(ρ·exc)/∂ρ_σ`) are corrected one at a time until the test passes cleanly.

## Domain basis

LSDA exchange (Dirac-Slater):
- `ε_x^σ(ρ_σ) = −C_x · (2ρ_σ)^{1/3}` per electron of spin σ, with `C_x = (3/4)(3/π)^{1/3} ≈ 0.7385587663820224` (Hartree atomic units).
- Spin combination: `ε_x(ρ_α, ρ_β) = (ρ_α · ε_x^α + ρ_β · ε_x^β) / ρ`.

LSDA correlation (Vosko-Wilk-Nusair 1980, parameterization V — "VWN5"):
- Padé form `ε_c(r_s, ζ) = ε_c^P(r_s) + α_c(r_s)·f(ζ)/f''(0) · (1 − ζ^4) + [ε_c^F(r_s) − ε_c^P(r_s)] · f(ζ) · ζ^4`, with `f(ζ) = ((1+ζ)^{4/3} + (1−ζ)^{4/3} − 2) / (2^{4/3} − 2)` and `f''(0) = 4/(9·(2^{1/3} − 1))`.
- Parameters in Hartree (A_Ha = A_Ry / 2): paramagnetic `{A=0.0310907, x0=−0.10498, b=3.72744, c=12.9352}`, ferromagnetic `{A=0.01554535, x0=−0.32500, b=7.06042, c=18.0578}`, spin stiffness `{A=−1/(6π²), x0=−0.0047584, b=1.13107, c=13.0045}`.
- Reference: Vosko, Wilk, Nusair, *Can. J. Phys.* **58**, 1200 (1980), DOI 10.1139/p80-159.

PBE exchange + correlation: Perdew, Burke, Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996), DOI 10.1103/PhysRevLett.77.3865; erratum DOI 10.1103/PhysRevLett.78.1396.

libxc output convention (ground truth):
- `vrho[σ] = ∂(ρ·ε_xc)/∂ρ_σ` (the XC potential `v_xc^σ`), not `∂ε_xc/∂ρ_σ`.
- `vsigma[k] = ∂(ρ·ε_xc)/∂σ_k` for `k ∈ {αα, αβ, ββ}`.
- Reference: Lehtola et al., *SoftwareX* **7**, 1 (2018), DOI 10.1016/j.softx.2017.11.002.
- Functional IDs: `XC_LDA_X` (1), `XC_LDA_C_VWN` (7, VWN5 — *not* `XC_LDA_C_VWN_3`=5), `XC_GGA_X_PBE` (101), `XC_GGA_C_PBE` (130).

All formulas live in Hartree atomic units.

## Design

The native `evaluate_xc()` (in `src/cpu/scf/xc_functional.cpp`) is the only production code modified. Signature, file location, header (`include/atomiverse/cpu/scf/xc_functional.h`), and `XCOutput` struct preserved; only internal arithmetic and any constants/parameters proven divergent from libxc are touched. No new symbols.

A new GoogleTest binary `test_xc_vs_libxc` is added under `tests/cpu/scf/` and registered with CTest. Gated by `#ifdef ATV_HAS_LIBXC`; when libxc is unavailable the TU compiles to an empty `int main()`. The test calls both `evaluate_xc()` (native) and `evaluate_xc_libxc()` (wrapper from sub-spec 02) at each tabulated `(ρ_α, ρ_β, σ_αα, σ_αβ, σ_ββ)` point, comparing all six output channels with `EXPECT_NEAR(..., 1e-12)`. Failure messages include functional, point label, channel name, native value, libxc value, and absolute difference.

The existing `tests/cpu/scf/test_xc_functional.cpp` is kept as a self-consistency unit test; its loose `kTol = 1e-6` is tightened where assertions are absolute equalities against libxc-derived expectations.

Bug fixes proceed in the order most likely to drive convergence:
1. vrho/vsigma convention — convert all native outputs to `∂(ρ·ε_xc)/∂ρ_σ` form, the single largest semantic gap.
2. `f(ζ)` and its derivative `df/dζ` — recompute with `(2^{4/3} − 2)` denominator and `f''(0) = 4/(9·(2^{1/3} − 1))` exactly.
3. VWN5 parameter table — verify each A is in Hartree (A_Ry/2), every `x0`, `b`, `c` sign and decimal vs. Table 4.4 of the VWN paper.
4. Slater `kCx` — recompute as `(3.0/4.0) * std::cbrt(3.0 / M_PI)` and confirm `2^{1/3}` spin factor.

After each fix the test is rerun; an item is only "done" when no point regresses any previously-passing channel.

## Files to create or modify

- `tests/cpu/scf/test_xc_vs_libxc.cpp` (new)
- `tests/cpu/scf/CMakeLists.txt` — register the new test target under `ATV_HAS_LIBXC` guard
- `src/cpu/scf/xc_functional.cpp` — fix vrho/vsigma convention, VWN5 spin-interp, VWN5 parameter table, Slater kCx
- `tests/cpu/scf/test_xc_functional.cpp` — update stale expected values once conventions changed

## Tasks

- [ ] Write failing parity test enumerating LSDA and PBE at all seven density/σ points (`tests/cpu/scf/test_xc_vs_libxc.cpp`)
- [ ] Register `test_xc_vs_libxc` binary under `ATV_HAS_LIBXC` guard (`tests/cpu/scf/CMakeLists.txt`)
- [ ] Fix vrho/vsigma output convention to `∂(ρ·exc)/∂ρ_σ` in `evaluate_xc` (`src/cpu/scf/xc_functional.cpp`)
- [ ] Correct VWN5 spin-interpolation `f(ζ)`, `f'(ζ)`, and parameter table in `evaluate_xc` (`src/cpu/scf/xc_functional.cpp`)
- [ ] Update stale expected values in `test_xc_functional.cpp` now that conventions changed
- [ ] Verify against libxc `XC_LDA_X+XC_LDA_C_VWN` and `XC_GGA_X_PBE+XC_GGA_C_PBE` at 1e-12 across all points
- [ ] Run full check + test suite

## Testing strategy

- **Happy path (LSDA)**: native vs. libxc agreement to 1e-12 on `(exc, vrho_alpha, vrho_beta)` at unpolarized `(ρ_α=ρ_β=0.5)`, partially polarized `(0.8, 0.2)`, fully polarized `(1.0, 0)`.
- **Edge cases (LSDA)**: low density `(ρ_α=1e-6, ρ_β=1e-8)` — must not be silently clamped above libxc `dens_threshold`; high density `(ρ_α=100, ρ_β=50)`.
- **Happy path (PBE)**: same three polarization points with σ sweep `σ_αα = σ_ββ ∈ {0.0, 0.01, 0.1, 1.0}` and `σ_αβ = (σ_αα + σ_ββ)/2`; all six channels to 1e-12.
- **Edge case (PBE)**: `σ = 0` collapses to LSDA-like density-only evaluation.
- **Domain validation**: the parity gate *is* the domain validation; every published functional identity is implicitly checked through libxc at all enumerated points. Per `science.required`, scientist-agent review at `/mol:impl` finalize confirms references and Hartree-unit consistency.

## Out of scope

- Grid-level XC integration (Lebedev / Treutler-Ahlrichs replacement) — chain sub-spec 04.
- End-to-end UKS vs. PySCF total-energy validation — chain sub-spec 05.
- XC functionals beyond LSDA (Slater + VWN5) and PBE; alternative parameterisations (VWN3, PW92, LYP).
- meta-GGA, hybrid, range-separated, non-collinear-spin.
- Performance optimisation — parity only; optimisation without a parity baseline risks re-introducing the bugs being removed.
