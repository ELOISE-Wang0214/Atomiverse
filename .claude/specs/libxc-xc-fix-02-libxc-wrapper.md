---
title: libxc oracle for native XC parity (test-only scaffolding)
slug: libxc-xc-fix-02-libxc-wrapper
status: approved
created: 2026-05-19
revised: 2026-05-19
depends_on:
  - libxc-xc-fix-01-cmake-integration
---

# libxc oracle for native XC parity (test-only scaffolding)

## Summary

Stand up a libxc-backed oracle that lives entirely under `tests/cpu/scf/oracle/` and is linked **only** into test targets — never into `libatomiverse`. The oracle exposes two free functions in the test-only namespace `ATV::test::xc_oracle` — `evaluate_lsda(rho_a, rho_b)` and `evaluate_pbe(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb)` — both returning Atomiverse's existing production `ATV::cpu::scf::XCOutput` struct so the parity test in sub-spec 03 can compare native and oracle field-by-field without translation. A new CMake option `ATV_BUILD_XC_ORACLE` (introduced in sub-spec 01) gates the entire `oracle/` directory: when OFF, no oracle source, header, or test target is configured. **This oracle is temporary scaffolding; a follow-up spec will delete it after sub-spec 03's primitive parity gate stays green.**

## Domain basis

XC functional identifiers (libxc 6.x) used by the oracle:

- `XC_LDA_X` — Dirac–Slater exchange. Dirac, *Proc. Cambridge Phil. Soc.* **26**, 376 (1930); Slater, *Phys. Rev.* **81**, 385 (1951). `ε_x = −C_x ρ^{1/3}`, `C_x = (3/4)(3/π)^{1/3} ≈ 0.7385587663820224 Hartree`.
- `XC_LDA_C_VWN` — id **7**, VWN parameterization V ("VWN5"). **Not** `XC_LDA_C_VWN_3` (id 5). Vosko, Wilk, Nusair, *Can. J. Phys.* **58**, 1200 (1980), DOI 10.1139/p80-159.
- `XC_GGA_X_PBE` (id 101) + `XC_GGA_C_PBE` (id 130). Perdew, Burke, Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996), DOI 10.1103/PhysRevLett.77.3865.
- libxc reference: Lehtola et al., *SoftwareX* **7**, 1 (2018), DOI 10.1016/j.softx.2017.11.002.

libxc output convention (the oracle preserves this exactly; sub-spec 03 fixes native to match):

| libxc array | meaning                          | XCOutput field        |
| ----------- | -------------------------------- | --------------------- |
| `exc`       | energy per electron (Hartree)    | `XCOutput.exc`        |
| `vrho[0]`   | `∂(ρ·exc)/∂ρ_α`                  | `XCOutput.vrho_alpha` |
| `vrho[1]`   | `∂(ρ·exc)/∂ρ_β`                  | `XCOutput.vrho_beta`  |
| `vsigma[0]` | `∂(ρ·exc)/∂σ_αα` (GGA)           | `XCOutput.vsigma_aa`  |
| `vsigma[1]` | `∂(ρ·exc)/∂σ_αβ` (GGA)           | `XCOutput.vsigma_ab`  |
| `vsigma[2]` | `∂(ρ·exc)/∂σ_ββ` (GGA)           | `XCOutput.vsigma_bb`  |

All calls use `XC_POLARIZED` mode. Same `kRhoMin = 1e-12` cutoff as native.

Units: densities in electrons/Bohr³, σ in (electrons/Bohr³)²/Bohr², energies in Hartree.

## Design

**Lifetime statement.** This oracle is test-only scaffolding. The `tests/cpu/scf/oracle/` directory, the `xc_libxc_oracle` static library, the `ATV_BUILD_XC_ORACLE` CMake option, and any libxc dependency it introduces are slated for removal by a follow-up "retire libxc oracle" spec once sub-spec 03's primitive parity gate has stayed green for an agreed soak period. Nothing in `include/atomiverse/`, `src/`, or any production-linked target may include the oracle header. The oracle never appears in `libatomiverse`.

**Namespace.** All oracle symbols live in `ATV::test::xc_oracle`. The production namespace `ATV::cpu::scf` is *not* extended. Reusing `XCOutput` from `include/atomiverse/cpu/scf/xc_functional.h` is intentional: the oracle populates that production struct from libxc's output buffers so the parity test compares apples to apples without translation. The oracle header includes the production header; the production header has no awareness of the oracle.

**Public surface** (header `tests/cpu/scf/oracle/xc_libxc_oracle.h`):

```cpp
#pragma once
#include "atomiverse/cpu/scf/xc_functional.h"  // for ATV::cpu::scf::XCOutput

namespace ATV { namespace test { namespace xc_oracle {

ATV::cpu::scf::XCOutput evaluate_lsda(double rho_a, double rho_b);
ATV::cpu::scf::XCOutput evaluate_pbe(double rho_a, double rho_b,
                                     double sigma_aa, double sigma_ab,
                                     double sigma_bb);

}}}  // namespace ATV::test::xc_oracle
```

The header is unconditional (no `#ifdef ATV_HAS_LIBXC`); the entire `oracle/` directory is added to the build only when `ATV_BUILD_XC_ORACLE=ON`, so this header is only ever compiled in that configuration.

**Per-call flow.**

1. If `rho_a + rho_b < kRhoMin = 1e-12`, return zero-initialized `XCOutput`.
2. `evaluate_lsda`: local `xc_func_init(&fx, XC_LDA_X, XC_POLARIZED)` and `xc_func_init(&fc, XC_LDA_C_VWN, XC_POLARIZED)`, evaluate via `xc_lda_exc_vxc`, **sum** `exc` and `vrho[σ]` across the two functionals, leave `vsigma_*` zero, `xc_func_end` on both.
3. `evaluate_pbe`: local `xc_func_init(&fx, XC_GGA_X_PBE, XC_POLARIZED)` and `xc_func_init(&fc, XC_GGA_C_PBE, XC_POLARIZED)`, evaluate via `xc_gga_exc_vxc`, sum `exc`, `vrho[σ]`, `vsigma[k]`, `xc_func_end` on both.
4. Pack libxc's interleaved arrays into the named fields of `XCOutput`.

Each call performs its own init/teardown; this is test scaffolding, not a hot path.

**Build wiring.** A new `tests/cpu/scf/oracle/CMakeLists.txt` declares:

```cmake
# Entire file body is gated by sub-spec 01's add_subdirectory(),
# which only runs when ATV_BUILD_XC_ORACLE=ON.
add_library(xc_libxc_oracle STATIC xc_libxc_oracle.cpp)
target_include_directories(xc_libxc_oracle PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(xc_libxc_oracle PUBLIC Libxc::xc)

# Smoke test for the oracle itself
atv_add_cpp_test(test_xc_libxc_oracle test_xc_libxc_oracle.cpp)
target_link_libraries(test_xc_libxc_oracle PRIVATE xc_libxc_oracle)
```

`libatomiverse` is not modified.

**Production purity invariant.** `libatomiverse.{a,so}` contains zero references to `xc_func_init`, `xc_lda_exc_vxc`, `xc_gga_exc_vxc`, or `xc_oracle::*` regardless of `ATV_BUILD_XC_ORACLE`. Asserted by ac-002.

## Files to create or modify

- `tests/cpu/scf/oracle/xc_libxc_oracle.h` (new) — declares `ATV::test::xc_oracle::evaluate_lsda` and `evaluate_pbe`.
- `tests/cpu/scf/oracle/xc_libxc_oracle.cpp` (new) — implementation; `<xc.h>` included only here.
- `tests/cpu/scf/oracle/test_xc_libxc_oracle.cpp` (new) — GoogleTest binary: smoke, spin polarization, density cutoff, Slater textbook value.
- `tests/cpu/scf/oracle/CMakeLists.txt` (new) — declares `xc_libxc_oracle` STATIC target + oracle test registration.

## Tasks

- [ ] Write failing oracle tests covering LSDA/PBE smoke, spin polarization, density cutoff, and the Slater textbook value (`tests/cpu/scf/oracle/test_xc_libxc_oracle.cpp`)
- [ ] Add `tests/cpu/scf/oracle/CMakeLists.txt` declaring the `xc_libxc_oracle` STATIC target linking `Libxc::xc` and the oracle test registration
- [ ] Implement `ATV::test::xc_oracle::evaluate_lsda` summing `XC_LDA_X + XC_LDA_C_VWN` (id 7, VWN5) in `XC_POLARIZED` mode with `kRhoMin = 1e-12` cutoff (`tests/cpu/scf/oracle/xc_libxc_oracle.{h,cpp}`)
- [ ] Implement `ATV::test::xc_oracle::evaluate_pbe` summing `XC_GGA_X_PBE + XC_GGA_C_PBE` with `vsigma_{aa,ab,bb}` packing (`tests/cpu/scf/oracle/xc_libxc_oracle.cpp`)
- [ ] Verify production purity — `nm` on `libatomiverse.{a,so}` reports zero `xc_func` / `xc_oracle::*` symbols regardless of `ATV_BUILD_XC_ORACLE` setting
- [ ] Run full check + test suite (`cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_XC_ORACLE=ON && cmake --build build -j && ctest --test-dir build --output-on-failure`)

## Testing strategy

- **Happy-path LSDA**: `xc_oracle::evaluate_lsda(0.5, 0.5)` returns finite `exc`, `vrho_alpha`, `vrho_beta`; `exc < 0`; `vsigma_*` exactly zero.
- **Happy-path PBE**: `xc_oracle::evaluate_pbe(0.5, 0.5, 0.1, 0.1, 0.1)` returns finite values on all six channels; `exc < 0`.
- **Spin polarization**: `xc_oracle::evaluate_lsda(0.8, 0.2)` yields `vrho_alpha != vrho_beta`.
- **Density cutoff**: `xc_oracle::evaluate_lsda(1e-20, 1e-20)` returns zero-initialized `XCOutput`.
- **Domain validation — Slater textbook**: at ρ_α=ρ_β=0.5 (ζ=0), an LDA-exchange-only sub-test calling `XC_LDA_X` alone asserts `|exc − (−0.7385587663820224)| < 1e-10 Ha`, matching `−(3/4)(3/π)^{1/3}` independent of VWN.
- **Production purity build matrix**: with `ATV_BUILD_XC_ORACLE=OFF` (default), `nm libatomiverse.a` shows zero `xc_func` symbols. With `ATV_BUILD_XC_ORACLE=ON`, count remains zero (oracle is its own static library, never linked into `atomiverse`).

## Out of scope

- Comparing oracle output against native `evaluate_xc()` — sub-spec 03's parity gate.
- Any libxc-touching code under `include/atomiverse/`, `src/`, or `python/` — forbidden by the test-only invariant.
- Cached `xc_func_type` handles, batched evaluators, thread-local pools.
- Meta-GGA, hybrid, range-separated, anything beyond {LSDA, PBE}.
- Deletion of the oracle directory and the `ATV_BUILD_XC_ORACLE` option — a separate "retire libxc oracle" spec performs this.
- CUDA backend XC.
