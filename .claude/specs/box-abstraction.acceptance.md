---
slug: box-abstraction
criteria:
  - id: ac-001
    summary: cpu::Box public API present at backend root
    type: code
    pass_when: |
      `include/atomiverse/cpu/box.h` defines class `ATV::cpu::Box` with
      these public members (declared, possibly inline-defined):
        - default ctor; explicit ctor from `const Mat3&`
        - static factory `orthorhombic(double a, double b, double c)`
        - static factory `cubic(double l)`
        - `double length(int d) const noexcept`
        - `double volume() const noexcept`
        - `bool is_orthorhombic() const noexcept`
        - `void min_image(double& dx, double& dy, double& dz) const noexcept`
        - `template <class B> void min_image_batch(B&, B&, B&) const noexcept`
        - `double wrap_axis(double q, int d) const noexcept`
        - `void wrap_into(double& x, double& y, double& z) const noexcept`
        - `void scale_diagonal(double, double, double)`
        - `void scale_isotropic(double)`
        - `void scale_axis(int d, double s)`
        - `void set_cell(const Mat3&)`
        - `const Mat3& cell() const noexcept`
        - `const Mat3& inv_cell() const noexcept`
      Each public member carries Doxygen `@brief` and (where applicable)
      `@param` annotations per CLAUDE.md. The ortho-fast-path members
      are marked `[[gnu::always_inline]]` (or platform equivalent).
    evidence: grep + Doxygen run

  - id: ac-002
    summary: cuda::Box at backend root with parity API
    type: code
    pass_when: |
      `include/atomiverse/cuda/box.cuh` defines `ATV::cuda::Box` as a
      trivially-copyable plain struct (verified by
      `static_assert(std::is_trivially_copyable_v<ATV::cuda::Box>)`),
      and provides `__host__ __device__` free functions matching the
      CPU API conceptually:
        - `min_image(T&, T&, T&, const Box&)`
        - `wrap_axis(T q, int d, const Box&)`
        - `length(const Box&, int)`
        - `volume(const Box&)`
      The header compiles in both `nvcc` and `clang++` (host-only)
      modes via `__host__ __device__` qualifiers.
    evidence: grep + static_assert compiles

  - id: ac-003
    summary: All duplicated min_image and wrap_to_box are removed
    type: code
    pass_when: |
      After landing,
        `git grep -nE 'static[[:space:]]+inline[[:space:]]+void[[:space:]]+minimum_image' src/cpu`
      returns zero matches, and
        `git grep -nE 'wrap_to_box' src/cpu`
      returns zero matches. The 4 duplicated copies of
      `static inline void minimum_image` (in lj.cpp, harmonic_bond.cpp,
      harmonic_angle.cpp, multi_harmonic_dihedral.cpp) and the 2
      duplicated copies of `wrap_to_box` (in mpi/rcb.cpp,
      mpi/atom_migration.cpp) are gone. Their callsites use
      `frame.box.min_image(...)` / `frame.box.wrap_axis(...)`.
    evidence: git grep before/after

  - id: ac-004
    summary: Raw Mat3 box reads outside Box are eliminated
    type: code
    pass_when: |
      `git grep -nE '\\.box\\([0-9],[[:space:]]*[0-9]\\)' src/cpu include/atomiverse/cpu`
      returns zero matches outside `include/atomiverse/cpu/box.h` and
      `src/cpu/box.cpp` themselves. (The implicit-conversion compat
      shim from Box to `const Mat3&` is removed by the closeout task,
      so any remaining raw read is a real bug.)
    evidence: git grep

  - id: ac-005
    summary: cpu::Frame::box field type migrated
    type: code
    pass_when: |
      `include/atomiverse/cpu/frame.h` declares the box field as
      `cpu::Box box{};` (not `Mat3 box{};`). Same for
      `include/atomiverse/cpu/frame_snapshot.h`.
      **CUDA Frame note (deferred):** `include/atomiverse/cuda/frame.cuh`
      keeps `prec_t* box` (device pointer to a 9-element row-major cell)
      because changing it to a host-resident `cuda::Box` value would
      require updating every kernel that today does `frame.box[i]` to
      take a `cuda::Box` by value or `__constant__`-memory upload. Tracked
      as a follow-up; the `cuda::Box` struct + `__host__ __device__`
      helpers are landed in this spec (see ac-002) so future CUDA
      components can adopt the unified API incrementally.
    evidence: grep

  - id: ac-006
    summary: Ortho min_image numerical parity to 1e-12
    type: numeric_oracle
    pass_when: |
      `tests/cpu/test_box.cpp` runs ≥ 10000 random `(L, dx)` cases per
      axis over orthorhombic boxes with L ∈ [5, 50] Å, dx ∈ [-2L, 2L],
      and asserts `|new - reference| < 1e-12` for each lane. The
      reference is the canonical branchy form `if (dx > L/2) dx -= L;
      else if (dx < -L/2) dx += L;` written inline in the test
      (independent of `Box::min_image`'s implementation). Test passes
      via `ctest -R test_box`.
    evidence: ctest output

  - id: ac-007
    summary: Triclinic min_image roundtrip identity
    type: numeric_oracle
    pass_when: |
      For two test cells with off-diagonal entries (one with a single
      shear in xy, one with shears in xy + xz), `box.min_image(dx, dy,
      dz)` produces a result `d'` such that
      `max_d |fractional(d')_d| ≤ 0.5 + 1e-12` (i.e. the fractional
      components fall in `[-0.5, 0.5]`), and the cache invariant
      `cell · inv_cell − I` is < 1e-12 in Frobenius norm. Verified in
      `tests/cpu/test_box.cpp` triclinic case. Test passes via
      `ctest -R test_box`.
    evidence: ctest output

  - id: ac-008
    summary: xsimd batch parity is lane-wise bit-exact
    type: numeric_oracle
    pass_when: |
      `tests/cpu/test_box.cpp` constructs an `xsimd::batch<double>`
      from a random vector, calls `min_image_batch`, calls scalar
      `min_image` per lane on the same input, and asserts lane-wise
      equality with **zero** tolerance (both paths reduce to the same
      IEEE ops). Repeated 1000× with fresh random inputs.
    evidence: ctest output

  - id: ac-009
    summary: Cache coherence after every mutation
    type: numeric_oracle
    pass_when: |
      In `tests/cpu/test_box.cpp`, after each of `set_cell`,
      `scale_diagonal`, `scale_isotropic`, `scale_axis`, the test
      asserts:
        - `||cell · inv_cell − I||_F < 1e-12`
        - `|cell(d,d) · inv_l[d] − 1| < 1e-12` for each ortho axis
        - `|2 · half_l[d] − cell(d,d)| < 1e-12` for each ortho axis
      All assertions pass.
    evidence: ctest output

  - id: ac-010
    summary: Full ctest suite still green
    type: code
    pass_when: |
      `cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
       cmake --build build -j$(nproc)
       ctest --test-dir build --output-on-failure`
      completes with no failures. Specific must-pass tests: test_md,
      test_bonded (any subtests), test_lj (any subtests), test_rdf,
      test_scf_evaluator, test_aimd_conservation, test_mpi_*
      (when MPI enabled), test_neighbor_list_cuda (when CUDA enabled).
    evidence: ctest log

  - id: ac-011
    summary: bench_md_step within 1% of pre-refactor baseline
    type: benchmark
    pass_when: |
      `bench_md_step` (LJ + bonded NVE, default 1k atoms / 100 steps)
      ns/step within ±1% of the pre-refactor baseline captured in
      `.claude/specs/box-abstraction.artifacts/perf.md` at 1 thread
      and 8 threads. Apple silicon M-series or x86-64 Zen3+. Three
      consecutive runs; report median.
    evidence: bench JSON in `.claude/specs/box-abstraction.artifacts/`

  - id: ac-012
    summary: bench_box throughput floors
    type: benchmark
    pass_when: |
      `benchmarks/bench_box` reports:
        - scalar `min_image` ns/call ≥ baseline static-inline form
          (≥ −5% margin)
        - `min_image_batch` ns/call ≥ 1.5× scalar at the platform's
          native `xsimd::batch<double>::size`
      Three consecutive runs; report median.
    evidence: bench JSON

  - id: ac-013
    summary: CUDA neighbor-list bench no regression
    type: benchmark
    pass_when: |
      `bench_neighbor_list` (or the project's CUDA NL bench)
      ns/call within ±1% of pre-refactor baseline. Required when
      ATV_BUILD_CUDA=ON; skipped otherwise. Skipped in this delivery
      because the dev environment has ATV_BUILD_CUDA=OFF (no nvcc).
      Triggered by the kernel-signature migration tracked under ac-005.
    evidence: bench JSON
    optional_if: ATV_BUILD_CUDA=OFF

  - id: ac-014
    summary: CLAUDE.md mentions cpu::Box / cuda::Box at backend root
    type: docs
    pass_when: |
      `CLAUDE.md` Architecture section names `cpu::Box` and
      `cuda::Box` and states: (a) they live at backend root as shared
      infrastructure, (b) they are the single source of truth for
      cell, periodic flags, min_image, and wrap_axis, (c) scattered
      `Mat3 box` reads and duplicated `static inline minimum_image`
      have been migrated.
    evidence: grep CLAUDE.md

  - id: ac-015
    summary: Doxygen clean for new public symbols
    type: docs
    pass_when: |
      Running the project's documentation build (`doxygen` or
      `make doc`) emits no `warning: ...` lines mentioning
      `ATV::cpu::Box` or `ATV::cuda::Box` symbols. Every public
      member of both classes has a `@brief`, every parameter that
      carries physical units has a `@param` line stating the unit.
    evidence: doxygen log
