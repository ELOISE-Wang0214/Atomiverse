---
slug: libxc-xc-fix-01-cmake-integration
criteria:
  - id: ac-001
    summary: Default-OFF build is libxc-free; production headers/sources gone
    type: code
    pass_when: |
      After this sub-spec lands, on a clean checkout with no
      `-DATV_BUILD_XC_ORACLE` flag, `cmake -B build -DATV_BUILD_TESTS=ON
      -DATV_BUILD_RHF=ON` exits 0; `build/CMakeCache.txt` contains
      `ATV_BUILD_XC_ORACLE:BOOL=OFF` and no `ATV_USE_LIBXC` entry;
      `grep -rE "xc_libxc|ATV_HAS_LIBXC|ATV_USE_LIBXC" include/ src/`
      returns zero matches; `grep -E "Libxc|libxc"
      build/compile_commands.json` returns zero matches; and the files
      `include/atomiverse/cpu/scf/xc_libxc.h` and
      `src/cpu/scf/xc_libxc.cpp` do not exist on disk.
    status: pending

  - id: ac-002
    summary: ON-path scopes libxc to tests; libatomiverse stays libxc-free
    type: code
    pass_when: |
      On a host with libxc 6+ installed, `cmake -B build-xc
      -DATV_BUILD_XC_ORACLE=ON -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON`
      exits 0; inspecting `build-xc/compile_commands.json` for every
      command whose `file` starts with `src/` or whose `output` is under
      `CMakeFiles/atomiverse.dir/` shows zero occurrences of any libxc
      include path or `-lxc` flag; if `tests/cpu/scf/oracle/CMakeLists.txt`
      exists (added by sub-spec 02), `add_subdirectory(cpu/scf/oracle)`
      runs and contributes targets to the test build; if it does not
      exist (this spec landing first), configure still succeeds with
      a clear STATUS or FATAL_ERROR.
    status: pending

  - id: ac-003
    summary: ON without libxc fails loudly with install hints
    type: code
    pass_when: |
      Running `cmake -B build-fail -DATV_BUILD_XC_ORACLE=ON
      -DLibxc_DIR=/definitely/does/not/exist -DCMAKE_PREFIX_PATH=""`
      (or equivalent scrub on a host without libxc) exits non-zero;
      combined stdout+stderr contains both the substring
      `brew install libxc` and the substring `apt install libxc-dev`;
      no FetchContent download of libxc is attempted (no `libxc-src`
      or `_deps/libxc-*` directory appears under the build tree); and
      `libatomiverse` is never configured for libxc linkage during the
      failed run.
    status: pending
---

# Acceptance criteria

- **ac-001** locks in the deletability invariant: with the option off (the default), the repo and configured build are clean of every libxc trace, and the production headers/sources leaked by earlier commits are gone. Grepping `include/` and `src/` is the load-bearing check.
- **ac-002** is the architectural guarantee: even when the oracle is built, `libatomiverse` is not linked, included, or compiled against libxc. `compile_commands.json` is the third-party-verifiable evidence.
- **ac-003** preserves no-silent-fallback policy and ensures the failure mode never partially configures `libatomiverse` against libxc.
