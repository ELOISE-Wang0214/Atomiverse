---
title: libxc CMake Integration (find_package only)
slug: libxc-xc-fix-01-cmake-integration
status: approved
created: 2026-05-19
---

# libxc CMake Integration (find_package only)

## Summary

Wire libxc into the Atomiverse build as a system-only optional dependency. When the user configures with `-DATV_USE_LIBXC=ON`, CMake resolves libxc via `find_package(Libxc 6 CONFIG REQUIRED)` and links `Libxc::xc` into the `atomiverse` target, defining the public compile macro `ATV_HAS_LIBXC`. When `find_package` fails, the configure step aborts with an actionable message naming the supported install paths (`brew install libxc` on macOS, `apt install libxc-dev` on Linux). When `ATV_USE_LIBXC=OFF` (default), the build behaves identically to today: no libxc lookup, no macro defined, no source file added. This sub-spec ships build infrastructure only — it does not yet add any C++ translation unit that includes `<xc.h>`; that work lands in sub-spec 02.

## Domain basis

No physics is exercised by this sub-spec. libxc itself implements XC functionals from Lehtola et al., *SoftwareX* 7 (2018) 1 (DOI: 10.1016/j.softx.2017.11.002); functional-level validation is covered in later chain sub-specs. Units of energy density returned by `xc_func.h` are Hartree/Bohr³ (libxc convention) and must be converted at the call site downstream — out of scope here.

## Design

Touched entities:

- `CMakeLists.txt` (root) — replace the existing `if(ATV_USE_LIBXC)` block with a find_package-only path. Keep the `option(ATV_USE_LIBXC ... OFF)` declaration unchanged.

New CMake symbols emitted by this sub-spec:

- Public compile definition `ATV_HAS_LIBXC` on the `atomiverse` target, conditional on `ATV_USE_LIBXC=ON` AND `Libxc_FOUND`.
- Public link dependency `Libxc::xc` on the `atomiverse` target under the same condition.
- No new C++ source files. The current `src/cpu/scf/xc_libxc.cpp` is **dropped from the source list** in this sub-spec (it will be re-added in sub-spec 02 once a libxc-using API exists).

Lifecycle: dependency resolution happens during CMake configure. No runtime cost, no `build()`/`launch()` impact on any component. The macro is visible to all C++ translation units once the `atomiverse` target is depended on, so downstream sub-specs can `#ifdef ATV_HAS_LIBXC` cleanly.

Failure mode: if the user sets `-DATV_USE_LIBXC=ON` and libxc is missing, `find_package(... REQUIRED)` plus a custom `message(FATAL_ERROR ...)` produce a single clear error. No FetchContent path — explicit user policy decision.

## Files to create or modify

- `CMakeLists.txt` — replace the `if(ATV_USE_LIBXC) ... endif()` block with a find_package-only block; drop the `src/cpu/scf/xc_libxc.cpp` entry inside it.
- `tests/CMakeLists.txt` — add a tiny configure-time check test that compiles a stub TU verifying `ATV_HAS_LIBXC` is defined when `ATV_USE_LIBXC=ON`.
- `tests/build/test_libxc_macro.cpp` (new) — one-line program asserting the macro is defined; built only when `ATV_USE_LIBXC=ON AND ATV_BUILD_TESTS=ON`.

## Tasks

- [ ] Write failing test `tests/build/test_libxc_macro.cpp` (new) that `static_assert`s `defined(ATV_HAS_LIBXC)` and register it in `tests/CMakeLists.txt` behind `if(ATV_USE_LIBXC)`
- [ ] Replace the `if(ATV_USE_LIBXC)` block in `CMakeLists.txt` with `find_package(Libxc 6 CONFIG REQUIRED)`, `target_link_libraries(atomiverse PUBLIC Libxc::xc)`, `target_compile_definitions(atomiverse PUBLIC ATV_HAS_LIBXC)`, and a `FATAL_ERROR` message guiding users to `brew install libxc` / `apt install libxc-dev`
- [ ] Remove `src/cpu/scf/xc_libxc.cpp` from the libxc block source list (sub-spec 02 will re-introduce it)
- [ ] Verify default-OFF behavior: run `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON` (no `ATV_USE_LIBXC`) and confirm CMake cache shows `ATV_USE_LIBXC:BOOL=OFF` and `ATV_HAS_LIBXC` does NOT appear in `compile_commands.json`
- [ ] Verify ON-with-system-libxc path: `cmake -B build-libxc -DATV_USE_LIBXC=ON -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON` succeeds on a host with libxc installed and the new build-macro test passes
- [ ] Run full check + test suite: `cmake --build build -j$(nproc) && ctest --test-dir build --output-on-failure`

## Testing strategy

- **Happy path (default OFF):** clean configure with no `ATV_USE_LIBXC` flag completes; cache value is OFF; no `Libxc` symbols appear in CMake configure log; `ctest` passes the existing suite untouched.
- **Happy path (ON, system libxc present):** configure with `-DATV_USE_LIBXC=ON` succeeds; `tests/build/test_libxc_macro.cpp` compiles and links; the `static_assert(defined(ATV_HAS_LIBXC))` proves the macro propagates as a `PUBLIC` definition.
- **Failure path (ON, no libxc):** `find_package(Libxc 6 CONFIG REQUIRED)` fails with a `FATAL_ERROR` message containing the substrings `brew install libxc` and `apt install libxc-dev`.
- **Regression:** the full `ctest` run with `ATV_USE_LIBXC=OFF` shows zero test count delta vs. the parent branch.

## Out of scope

- Any `.cpp` or `.h` change that includes `<xc.h>` (sub-spec 02).
- A C++ wrapper class around `xc_func_type` (sub-spec 02).
- Switching `xc_grid_engine` over to libxc for ground-truth comparison (sub-spec 03).
- The native-XC bug fixes in `grid.cpp` (sub-spec 04).
- libxc version-range probing beyond `find_package(Libxc 6 ...)`.
- FetchContent fallback — explicitly rejected by user policy on this revision.
