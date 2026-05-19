---
title: Native XC parity vs. libxc to 1e-12
slug: libxc-xc-fix-03-native-parity
status: approved
created: 2026-05-19
revised: 2026-05-19
depends_on:
  - libxc-xc-fix-01-cmake-integration
  - libxc-xc-fix-02-libxc-wrapper
---

# Native XC parity vs. libxc to 1e-12

## Summary

This sub-spec is where Atomiverse's native `evaluate_xc()` in `src/cpu/scf/xc_functional.cpp` is corrected to numerically reproduce libxc to 1e-12 in `(exc, vrho_alpha, vrho_beta, vsigma_aa, vsigma_ab, vsigma_bb)` for LSDA (Slater + VWN5) and PBE at every point of a curated battery. The work is an **honest reimplementation**: we read libxc's source — the Maple-derived primitives in `libxc/maple/` and the runtime routines in `libxc/src/` — understand each primitive (Slater exchange constant `C_x`, VWN5 Padé parameters, spin-interpolation `f(ζ)` and `f''(0)`, PBE κ/μ/β/γ, the `∂(ρ·exc)/∂ρ_σ` vrho convention), and copy the equivalent arithmetic into `xc_functional.cpp`. The libxc-derived **test-only oracle** built in sub-spec 02 (`tests/cpu/scf/oracle/xc_libxc_oracle.h`, namespace `ATV::test::xc_oracle`) serves as ground-truth comparator; libxc is **not wrapped or called at runtime** by any production code path. After this sub-spec ships, the native code stands alone.

## Domain basis

LSDA exchange (Dirac–Slater):
- `ε_x^σ(ρ_σ) = −C_x · (2ρ_σ)^{1/3}` per electron of spin σ, with `C_x = (3/4)(3/π)^{1/3} ≈ 0.7385587663820224` (Hartree atomic units).
- Spin combination: `ε_x(ρ_α, ρ_β) = (ρ_α · ε_x^α + ρ_β · ε_x^β) / ρ`.

LSDA correlation (Vosko-Wilk-Nusair 1980, parameterization V — "VWN5"):
- Padé form `ε_c(r_s, ζ) = ε_c^P(r_s) + α_c(r_s) · f(ζ)/f''(0) · (1 − ζ^4) + [ε_c^F(r_s) − ε_c^P(r_s)] · f(ζ) · ζ^4`,
  with `f(ζ) = ((1+ζ)^{4/3} + (1−ζ)^{4/3} − 2) / (2^{4/3} − 2)` and `f''(0) = 4 / (9·(2^{1/3} − 1))`.
- Parameters in Hartree (A_Ha = A_Ry/2): paramagnetic `{A=0.0310907, x0=−0.10498, b=3.72744, c=12.9352}`, ferromagnetic `{A=0.01554535, x0=−0.32500, b=7.06042, c=18.0578}`, spin-stiffness `{A=−1/(6π²), x0=−0.0047584, b=1.13107, c=13.0045}`.
- Reference: Vosko, Wilk, Nusair, *Can. J. Phys.* **58**, 1200 (1980), DOI 10.1139/p80-159.

PBE: `F_x(s) = 1 + κ − κ/(1 + μ s²/κ)`, κ=0.804, μ=0.21951. Correlation β=0.066725, γ=(1−ln 2)/π² ≈ 0.031091. Perdew, Burke, Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996), DOI 10.1103/PhysRevLett.77.3865.

libxc convention (the ground truth we adopt):
- `vrho[σ] = ∂(ρ·ε_xc)/∂ρ_σ`, **not** `∂ε_xc/∂ρ_σ`.
- `vsigma[k] = ∂(ρ·ε_xc)/∂σ_k`.
- Lehtola et al., *SoftwareX* **7**, 1 (2018), DOI 10.1016/j.softx.2017.11.002.
- Functional IDs (used by oracle): `XC_LDA_X` (1), `XC_LDA_C_VWN` (7, VWN5 — not VWN3=5), `XC_GGA_X_PBE` (101), `XC_GGA_C_PBE` (130).

All formulas in Hartree atomic units.

## Design

Native `evaluate_xc()` in `src/cpu/scf/xc_functional.cpp` is the only production code modified. Signature, file location, header, and `XCOutput` struct preserved; only internal arithmetic and any constants/parameters proven divergent from libxc are touched. No new public symbols, no new namespaces, no new build flags in production. **Crucially, no `#ifdef ATV_HAS_LIBXC` or `#ifdef ATV_USE_LIBXC` appears anywhere in production sources** — after this sub-spec, the native implementation has no compile-time dependency on libxc whatsoever.

The work proceeds test-first. A GoogleTest binary `test_cpu_scf_xc_vs_libxc` (rewired to live under `tests/cpu/scf/` and link the test-only `xc_libxc_oracle` target from sub-spec 02) calls both `ATV::cpu::scf::evaluate_xc()` (native) and `ATV::test::xc_oracle::evaluate_lsda` / `evaluate_pbe` (test-only oracle) at each tabulated point. The oracle is the **only** site that includes libxc headers; production code remains libxc-free.

Bug fixes proceed in the order most likely to drive convergence — the four root causes diagnosed during the sub-spec 02 wrapper work:

1. **vrho/vsigma convention** — convert all native outputs to `∂(ρ·ε_xc)/∂ρ_σ` form; the single largest semantic gap.
2. **VWN5 `f(ζ)` and `f'(ζ)`** — recompute with denominator `(2^{4/3} − 2)` and `f''(0) = 4/(9·(2^{1/3} − 1))` exactly.
3. **VWN5 parameter table** — verify each `A` is in Hartree (A_Ry/2), every `x0`, `b`, `c` sign and decimal vs. Table 4.4 of the VWN paper.
4. **Slater `kCx`** — recompute as `(3.0/4.0) * std::cbrt(3.0 / M_PI)` and confirm `2^{1/3}` spin factor.

Each primitive copied from libxc is documented in-place: a comment block at the top of the modified region of `xc_functional.cpp` cites the libxc source file and line range we studied (`libxc/maple/lda_x.mpl`, `libxc/maple/lda_c_vwn.mpl`, `libxc/src/work_lda.c`, `libxc/maple/gga_x_pbe.mpl`, `libxc/maple/gga_c_pbe.mpl`). Honesty marker — proving honest reimplementation — and an audit aid for future XC additions.

After this spec lands and the parity tests are green, the oracle directory `tests/cpu/scf/oracle/` is in principle deletable: it has served its diagnostic role and the native code is self-sufficient. A follow-up retirement spec removes it once parity has remained green in CI for a stabilization period; that retirement is not part of this sub-spec.

## Files to create or modify

- `tests/cpu/scf/test_xc_vs_libxc.cpp` — fill in the parity battery; includes `tests/cpu/scf/oracle/xc_libxc_oracle.h` from sub-spec 02
- `tests/cpu/scf/CMakeLists.txt` (or `tests/CMakeLists.txt` depending on layout) — re-gate `test_cpu_scf_xc_vs_libxc` registration on `ATV_BUILD_XC_ORACLE` (currently gated on the obsolete `ATV_USE_LIBXC`); link the `xc_libxc_oracle` target
- `src/cpu/scf/xc_functional.cpp` — fix vrho/vsigma convention, VWN5 spin-interp + parameter table, Slater `kCx`; add libxc-citation comment block documenting each copied primitive
- `tests/cpu/scf/test_xc_functional.cpp` — update stale expected values once conventions changed

## Tasks

- [ ] Write failing parity tests in `tests/cpu/scf/test_xc_vs_libxc.cpp` enumerating LSDA and PBE at all battery points, calling `ATV::test::xc_oracle::evaluate_lsda` / `evaluate_pbe`
- [ ] Re-gate `test_cpu_scf_xc_vs_libxc` registration on `ATV_BUILD_XC_ORACLE` and link the `xc_libxc_oracle` target (test CMakeLists)
- [ ] Implement vrho/vsigma output-convention fix to `∂(ρ·exc)/∂ρ_σ` in `evaluate_xc` (`src/cpu/scf/xc_functional.cpp`)
- [ ] Implement VWN5 `f(ζ)`, `f'(ζ)`, and parameter-table corrections in `evaluate_xc`
- [ ] Implement Slater `kCx` correction and confirm `2^{1/3}` spin factor in `evaluate_xc`
- [ ] Document each primitive copied from libxc — file-line citation block at the top of the modified region of `xc_functional.cpp` (Maple sources in `libxc/maple/` and runtime helpers in `libxc/src/`)
- [ ] Update stale expected values in `tests/cpu/scf/test_xc_functional.cpp`
- [ ] Verify native vs. oracle agreement to 1e-12 across the full battery for LSDA and PBE
- [ ] Run full check + test suite

## Testing strategy

- **Happy path (LSDA)**: native vs. oracle agreement to 1e-12 on `(exc, vrho_alpha, vrho_beta)` at unpolarized `(0.5, 0.5)`, partially polarized `(0.8, 0.2)`, fully polarized `(1.0, 0)`.
- **Edge cases (LSDA)**: low density `(1e-6, 1e-8)` and high density `(100, 50)`.
- **Happy path (PBE)**: same polarization points × σ sweep `σ_αα = σ_ββ ∈ {0.0, 0.01, 0.1, 1.0}`, `σ_αβ = (σ_αα + σ_ββ)/2`; all six channels to 1e-12.
- **Edge case (PBE)**: `σ = 0` collapses to LSDA-like density-only evaluation.
- **Build-off matrix**: with `-DATV_BUILD_XC_ORACLE=OFF`, the parity test is not registered and no production translation unit references the oracle; full build + ctest still pass.
- **Domain validation**: the parity gate *is* the domain validation; per `science.required`, scientist-agent review at `/mol:impl` finalize confirms cited references, Hartree-unit consistency, and that each in-source libxc citation points at the actual primitive being copied.

## Out of scope

- Retiring `tests/cpu/scf/oracle/` — follow-up spec.
- Grid-level XC integration — sub-spec 04.
- End-to-end UKS vs. PySCF — sub-spec 05.
- XC functionals beyond LSDA + PBE; alternative parameterisations (VWN3, PW92, LYP).
- meta-GGA, hybrid, range-separated, non-collinear-spin.
- Performance optimisation — parity only.
