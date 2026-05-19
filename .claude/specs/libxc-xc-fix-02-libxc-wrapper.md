---
title: libxc XC wrapper (LSDA + PBE)
slug: libxc-xc-fix-02-libxc-wrapper
status: approved
created: 2026-05-19
---

# libxc XC wrapper (LSDA + PBE)

## Summary

Stand up `evaluate_xc_libxc()` — a libxc-backed evaluator that mirrors the signature of the existing native `evaluate_xc()` in `include/atomiverse/cpu/scf/xc_functional.h`. The wrapper supports LSDA (Slater + VWN5) and PBE (PBE exchange + PBE correlation) in spin-polarized mode, returning Atomiverse's `XCOutput` struct unchanged so callers do not need to know whether the backend is native or libxc. All libxc-touching code lives behind `#ifdef ATV_HAS_LIBXC` (defined by sub-spec 01); when libxc is disabled, the translation unit compiles to an empty object. This sub-spec does **not** validate libxc output against the native implementation — that comparison is sub-spec 03's job.

## Domain basis

XC functional identifiers (libxc 6.x):

- `XC_LDA_X` — Dirac–Slater exchange. Dirac, *Proc. Cambridge Phil. Soc.* **26**, 376 (1930); Slater, *Phys. Rev.* **81**, 385 (1951). Energy density per electron: ε_x = −C_x ρ^{1/3} with C_x = (3/4)(3/π)^{1/3}.
- `XC_LDA_C_VWN` (functional id **7**, VWN5 — **not** id 5 which is VWN3 RPA). Vosko, Wilk, Nusair, *Can. J. Phys.* **58**, 1200 (1980), parameterization V.
- `XC_GGA_X_PBE` + `XC_GGA_C_PBE`. Perdew, Burke, Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996), DOI: 10.1103/PhysRevLett.77.3865.

libxc output conventions (must be translated to Atomiverse's `XCOutput`):

| libxc array      | meaning                                  | Atomiverse field        |
| ---------------- | ---------------------------------------- | ----------------------- |
| `exc`            | energy per electron (Hartree)            | `XCOutput.exc`          |
| `vrho[0]`        | ∂(ρ·exc)/∂ρ_α                            | `XCOutput.vrho_alpha`   |
| `vrho[1]`        | ∂(ρ·exc)/∂ρ_β                            | `XCOutput.vrho_beta`    |
| `vsigma[0]`      | ∂(ρ·exc)/∂σ_αα  (GGA only)               | `XCOutput.vsigma_aa`    |
| `vsigma[1]`      | ∂(ρ·exc)/∂σ_αβ  (GGA only)               | `XCOutput.vsigma_ab`    |
| `vsigma[2]`      | ∂(ρ·exc)/∂σ_ββ  (GGA only)               | `XCOutput.vsigma_bb`    |

Units: densities in electrons/Bohr³, σ in (electrons/Bohr³)²/Bohr², energies in Hartree — identical between libxc and Atomiverse, no scale conversion required.

Spin-polarized initialization (libxc C API):

```c
xc_func_type func;
xc_func_init(&func, XC_LDA_X, XC_POLARIZED);
/* ... xc_lda_exc_vxc(&func, 1, rho, exc, vrho); ... */
xc_func_end(&func);
```

Atomiverse always calls libxc in `XC_POLARIZED` mode so the wrapper can serve both closed-shell and open-shell paths uniformly.

## Design

New namespace-level free function in `ATV::cpu::scf`:

```cpp
#ifdef ATV_HAS_LIBXC
XCOutput evaluate_xc_libxc(XCType type,
                           double rho_a, double rho_b,
                           double sigma_aa = 0, double sigma_ab = 0,
                           double sigma_bb = 0);
#endif
```

The header `xc_libxc.h` already exists with this prototype; this sub-spec fleshes out the implementation in `xc_libxc.cpp`. No new public types are introduced — `XCOutput` and `XCType` continue to live in `xc_functional.h`, and the libxc header is included only from the `.cpp`.

Lifecycle. Each call performs a local `xc_func_init` / `xc_func_end` pair. Simplest correct lifetime; keeps the wrapper thread-safe and stateless. A batched / cached-handle variant is deferred until the production XC loop adopts libxc (out of scope).

Internal flow per call:

1. Apply `kRhoMin = 1e-12` total-density cutoff (return zero-initialized `XCOutput` below cutoff) so the two backends agree on the floor before libxc is consulted.
2. `XCType::LSDA`: initialize `XC_LDA_X`, `XC_LDA_C_VWN`, evaluate each, **sum** `exc` and `vrho`. `vsigma_*` zero (LDA has no σ dependence).
3. `XCType::PBE`: initialize `XC_GGA_X_PBE`, `XC_GGA_C_PBE`, evaluate each, sum `exc`, `vrho`, `vsigma`.
4. Pack libxc's interleaved arrays into the named fields of `XCOutput`.

When `ATV_HAS_LIBXC` is undefined the TU reduces to no symbols — no fallback, no stub returning zeros.

## Files to create or modify

- `include/atomiverse/cpu/scf/xc_libxc.h` — header already exists; verify Doxygen mentions `XC_POLARIZED` and VWN5 id=7.
- `src/cpu/scf/xc_libxc.cpp` — implementation: replace any placeholder body with LSDA/PBE branches.
- `tests/cpu/scf/test_xc_libxc.cpp` (new) — gated by `#ifdef ATV_HAS_LIBXC`; covers init/teardown, finiteness of `exc`, Slater-exchange textbook value.
- `tests/CMakeLists.txt` — register `test_xc_libxc` only when `ATV_HAS_LIBXC` is defined.

## Tasks

- [ ] Write failing tests for `evaluate_xc_libxc` LSDA + PBE sanity and Slater textbook value (`tests/cpu/scf/test_xc_libxc.cpp`)
- [ ] Register `test_xc_libxc` in the test build, gated on `ATV_HAS_LIBXC` (`tests/CMakeLists.txt`)
- [ ] Implement `evaluate_xc_libxc` LSDA branch (`XC_LDA_X` + `XC_LDA_C_VWN`, id=7) with `XC_POLARIZED` init/end and convention mapping (`src/cpu/scf/xc_libxc.cpp`)
- [ ] Implement `evaluate_xc_libxc` PBE branch (`XC_GGA_X_PBE` + `XC_GGA_C_PBE`) including `vsigma_{aa,ab,bb}` packing (`src/cpu/scf/xc_libxc.cpp`)
- [ ] Confirm `xc_libxc.h` Doxygen documents libxc id=7 (VWN5, not VWN3) and that the TU is empty without `ATV_HAS_LIBXC`
- [ ] Run full check + test suite (`ctest --test-dir build --output-on-failure`)

## Testing strategy

All tests guarded by `#ifdef ATV_HAS_LIBXC` and skip when libxc is unavailable.

- **Happy-path LSDA.** `evaluate_xc_libxc(XCType::LSDA, 0.5, 0.5)` returns finite `exc`, `vrho_alpha`, `vrho_beta`; `exc < 0`; `vsigma_*` exactly zero.
- **Happy-path PBE.** `evaluate_xc_libxc(XCType::PBE, 0.5, 0.5, 0.1, 0.1, 0.1)` returns finite values everywhere; `exc < 0`; all three `vsigma_*` finite.
- **Spin-polarized smoke.** `(0.8, 0.2)` produces `vrho_alpha != vrho_beta`.
- **Density cutoff.** `(1e-20, 1e-20)` returns zero-initialized `XCOutput`.
- **Domain validation — Slater exchange textbook value.** At ρ_α = ρ_β = 0.5 (total ρ = 1, ζ = 0), evaluating libxc `XC_LDA_X` alone in a focused sub-test must yield `|exc − (−0.7385587663820224)| < 1e-10 Ha`, matching the closed form `−(3/4)(3/π)^{1/3}`.
- **Build matrix.** A test build with `-DATV_USE_LIBXC=OFF` must still link cleanly: `xc_libxc.cpp` produces no symbols.

## Out of scope

- Comparing libxc output against native `evaluate_xc()` to 1e-12 — sub-spec 03.
- Batched evaluation, cached libxc handles, thread-local functional pools — re-init per call is acceptable until profiling shows otherwise.
- Meta-GGA, hybrid, range-separated, anything outside {LSDA, PBE}.
- CUDA backend XC — out of scope for the whole chain.
- The runtime `use_libxc` switch in `uks_solve` — handled, if at all, by a later sub-spec.
