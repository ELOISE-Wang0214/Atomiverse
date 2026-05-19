---
slug: libxc-xc-fix-01-cmake-integration
criteria:
  - id: ac-001
    summary: Default-OFF configure produces no libxc traces
    type: code
    pass_when: |
      Running `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON`
      on a clean checkout (no `-DATV_USE_LIBXC` flag) completes with
      exit code 0; `build/CMakeCache.txt` contains
      `ATV_USE_LIBXC:BOOL=OFF`; grepping `build/compile_commands.json`
      for `ATV_HAS_LIBXC` returns zero matches; and grepping the
      CMake configure stdout for `Libxc` returns zero matches.
    status: pending

  - id: ac-002
    summary: ON-path defines ATV_HAS_LIBXC and links Libxc::xc
    type: code
    pass_when: |
      On a host with libxc 6+ installed (e.g. `brew install libxc`),
      `cmake -B build-libxc -DATV_USE_LIBXC=ON -DATV_BUILD_TESTS=ON
      -DATV_BUILD_RHF=ON` exits 0; `cmake --build build-libxc -j`
      succeeds; `tests/build/test_libxc_macro.cpp` (which contains
      `static_assert(defined(ATV_HAS_LIBXC), "macro not propagated");`)
      compiles; and `ctest --test-dir build-libxc -R libxc_macro
      --output-on-failure` reports 1/1 passing.
    status: pending

  - id: ac-003
    summary: Missing libxc fails loudly with install hints
    type: code
    pass_when: |
      Running `cmake -B build-fail -DATV_USE_LIBXC=ON
      -DLibxc_DIR=/definitely/does/not/exist` (or equivalent
      environment scrub on a host without libxc) exits non-zero;
      stderr contains both the substring `brew install libxc` and
      the substring `apt install libxc-dev`; no FetchContent
      download of libxc is attempted (no `libxc-src` directory
      appears under the build tree).
    status: pending
---

# Acceptance criteria

- **ac-001** establishes the no-regression invariant: existing users who do not opt in see byte-for-byte identical builds.
- **ac-002** proves the macro and link dependency reach the `atomiverse` target's public interface so that sub-spec 02 can rely on `#ifdef ATV_HAS_LIBXC` in any downstream TU.
- **ac-003** locks in the user policy choice: no silent FetchContent fallback. The error message must guide the user to a one-line fix on either of the two supported platforms.
