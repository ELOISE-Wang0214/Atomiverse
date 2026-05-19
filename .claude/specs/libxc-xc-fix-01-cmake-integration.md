---
title: libxc test-oracle CMake scaffolding
slug: libxc-xc-fix-01-cmake-integration
status: approved
created: 2026-05-19
revised: 2026-05-19
---

# libxc test-oracle CMake scaffolding

## Summary

libxc is **not** an Atomiverse dependency. It is a test-only oracle used to validate the native XC implementation. After native parity is locked in (sub-spec 03), the entire CMake option, the oracle code, and every reference to libxc must be deletable from the repository in a single commit without touching `libatomiverse`. This sub-spec rewrites the build wiring so that:

- A single option `ATV_BUILD_XC_ORACLE` (default OFF) gates everything libxc-related.
- libxc is linked **only into specific test targets**, never into `libatomiverse`.
- No public compile macro (no `ATV_HAS_LIBXC` on `atomiverse`).
- The wrongly-placed files `include/atomiverse/cpu/scf/xc_libxc.h` and `src/cpu/scf/xc_libxc.cpp` (leaked into production by earlier commits) are deleted; sub-spec 02 will re-create the oracle under `tests/cpu/scf/oracle/`.

This spec's existence is temporary. Once sub-spec 03 locks primitive parity, a follow-up retirement spec removes `ATV_BUILD_XC_ORACLE` and the oracle entirely.

## Design

**Root `CMakeLists.txt`:**
- Remove the obsolete `option(ATV_USE_LIBXC ...)` declaration and the entire `if(ATV_USE_LIBXC) ... endif()` block (libatomiverse-link path + any FetchContent fallback).
- Introduce `option(ATV_BUILD_XC_ORACLE "Build libxc-backed XC oracle (test-only scaffolding; deletable after native parity is locked in)" OFF)`.
- When ON: `find_package(Libxc 6 CONFIG REQUIRED)`. On failure, emit `FATAL_ERROR` containing both `brew install libxc` and `apt install libxc-dev`.
- **Do not** modify the `atomiverse` target. No `target_link_libraries(atomiverse ...)`, `target_include_directories(atomiverse ...)`, or `target_compile_definitions(atomiverse ...)` referencing libxc.

**`tests/CMakeLists.txt`:**
- Replace the existing `if(ATV_USE_LIBXC) atv_add_cpp_test(test_cpu_scf_xc_vs_libxc ...)` block with `if(ATV_BUILD_XC_ORACLE) add_subdirectory(cpu/scf/oracle) endif()`. Guard with `if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cpu/scf/oracle/CMakeLists.txt")` so this sub-spec lands cleanly on its own (sub-spec 02 creates the oracle subdirectory).

**File deletions:**
- `include/atomiverse/cpu/scf/xc_libxc.h` — was a public header leaked into the production tree.
- `src/cpu/scf/xc_libxc.cpp` — was a libatomiverse TU.

Compile-time macro policy: no public macro on `libatomiverse`. If sub-spec 02's oracle needs a guard, it lives as a PRIVATE compile definition on the oracle target, scoped to the test build only.

Failure-mode contract: `ATV_BUILD_XC_ORACLE=ON` + missing libxc → one clear `FATAL_ERROR` from `find_package(... REQUIRED)` + our install-hint message. No silent fallback, no FetchContent.

## Files to create or modify

- `CMakeLists.txt` — remove the wrong-arch libxc block; introduce `ATV_BUILD_XC_ORACLE`; add gated `find_package(Libxc 6 CONFIG REQUIRED)` with install-hint FATAL_ERROR.
- `tests/CMakeLists.txt` — replace the old `if(ATV_USE_LIBXC)` block with a guarded `if(ATV_BUILD_XC_ORACLE) add_subdirectory(cpu/scf/oracle) endif()`.
- `include/atomiverse/cpu/scf/xc_libxc.h` — **DELETE**.
- `src/cpu/scf/xc_libxc.cpp` — **DELETE**.

## Tasks

- [ ] Delete `include/atomiverse/cpu/scf/xc_libxc.h` and `src/cpu/scf/xc_libxc.cpp` from disk; remove the source file from any source list it appears in
- [ ] Remove the obsolete `option(ATV_USE_LIBXC ...)` declaration and the entire `if(ATV_USE_LIBXC) ... endif()` block from `CMakeLists.txt`
- [ ] Add `option(ATV_BUILD_XC_ORACLE "..." OFF)` and a gated `find_package(Libxc 6 CONFIG REQUIRED)` invocation with `FATAL_ERROR` install hints mentioning `brew install libxc` and `apt install libxc-dev` — without touching the `atomiverse` target
- [ ] Update `tests/CMakeLists.txt` to replace the old `if(ATV_USE_LIBXC)` block with `if(ATV_BUILD_XC_ORACLE) add_subdirectory(cpu/scf/oracle) endif()`, guarded by an `EXISTS` check on the oracle CMakeLists
- [ ] Verify three behaviors: default-OFF (no libxc traces; `libatomiverse` link line contains no `xc` symbol), ON-with-libxc (configure clean; `libatomiverse` still has no libxc include paths in `compile_commands.json`), ON-without-libxc (FATAL_ERROR with both install hints)
- [ ] Run full check + test suite (`cmake --build build -j && ctest --test-dir build --output-on-failure`) with `ATV_BUILD_XC_ORACLE` left OFF

## Testing strategy

- **Default-OFF regression invariant**: clean configure with no `ATV_BUILD_XC_ORACLE` flag produces a cache identical to current behavior; cache contains `ATV_BUILD_XC_ORACLE:BOOL=OFF`; `compile_commands.json` has no occurrence of `ATV_HAS_LIBXC`, `ATV_HAS_XC_ORACLE`, `libxc`, or `Libxc` in any `atomiverse` TU; full ctest count unchanged.
- **ON, libxc present**: configure with `-DATV_BUILD_XC_ORACLE=ON` succeeds on a host with libxc installed; `compile_commands.json` for any `atomiverse` TU contains no libxc include path.
- **ON, no libxc**: configure exits non-zero; stderr contains both install hints; no FetchContent fallback.
- **Production-leak invariant**: `grep -r "xc_libxc" include/ src/` and `grep -r "ATV_HAS_LIBXC\|ATV_USE_LIBXC" include/ src/` both return zero matches.

## Out of scope

- The oracle wrapper itself (`evaluate_lsda`, `evaluate_pbe`) — sub-spec 02.
- The native parity test body — sub-spec 03.
- Native XC bug fixes — sub-spec 04.
- libxc version-range probing beyond `find_package(Libxc 6 ...)`.
- FetchContent fallback — explicitly rejected by user policy and reaffirmed by the test-only architecture.
- User-facing documentation of `ATV_BUILD_XC_ORACLE` — by design, this option is developer-facing only.
