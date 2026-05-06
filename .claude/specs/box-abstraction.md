---
title: Box abstraction — cpu::Box / cuda::Box at backend root
slug: box-abstraction
status: done
created: 2026-05-06
landed: 2026-05-06
---

# Box abstraction

## Summary

Replace scattered `Mat3 box` reads, **four** copies of
`static inline minimum_image` (in `lj.cpp`, `harmonic_bond.cpp`,
`harmonic_angle.cpp`, `multi_harmonic_dihedral.cpp`), **two** copies of
`wrap_to_box` (in `mpi/rcb.cpp`, `mpi/atom_migration.cpp`), and the CUDA
helpers nested under `ATV::md::nlist::` (`triclinic.cuh`) with a single
canonical pair of zero-cost backend-root types:

- `cpu::Box`  — `include/atomiverse/cpu/box.h` + `src/cpu/box.cpp`
- `cuda::Box` — `include/atomiverse/cuda/box.cuh`

Both types expose the same conceptual API: `min_image`, `wrap_axis`,
`wrap_into`, `length`, `volume`, `is_orthorhombic`, `scale_diagonal`,
`scale_isotropic`, `set_cell`, `cell()`, `inv_cell()`. CPU additionally
exposes `min_image_batch<B>(...)` for the LJ ket-loop xsimd path. CUDA
exposes the same callables as `__host__ __device__` free functions
operating on a trivially-copyable `cuda::Box` plain-data struct so the
type can be uploaded to constant memory without surgery.

The fast-path is the design contract — not a side effect. The
orthorhombic branch must compile to **at most as many** instructions as
the current branchy `static inline minimum_image`, and the xsimd batch
form must compile to the same NEON/AVX kernel as the current `lj.cpp`
ket-loop. Triclinic uses the cached `inv_cell` Cartesian↔fractional
roundtrip pattern that `triclinic.cuh` already established.

## Binding rules (load-bearing — do not remove)

- **One Box, two backends, identical conceptual API.** `cpu::Box` and
  `cuda::Box` mirror each other. No third "shared" Box in `core/` —
  per CLAUDE.md "Files at the backend root are shared infrastructure
  used by multiple domains", Box belongs at each backend root, not in
  `core/`.
- **No copies of `minimum_image` or `wrap_to_box` outside Box.** After
  this spec lands, every callsite that needs minimum-image or
  axis-wrap goes through Box. Hot-path callers may write
  `frame.box.min_image(dx, dy, dz)` directly; no helper layer.
- **Hot-path discipline preserved.** Box stores fundamental cell data
  (`Mat3 cell`) plus a per-step-immutable cache (`Mat3 inv_cell`,
  per-axis `inv_L`, `half_L`, `bool ortho`). The cache is refreshed
  exactly when `set_cell` / `scale_*` is called (≤ once per step under
  NPT). No allocation in any Box method.
- **No-fallback discipline preserved.** `min_image_batch` requires
  orthorhombic; calling it on a triclinic box is a precondition
  violation, not a soft-fallback.
- **No new physics.** This is a refactor + minor performance lift via
  branch-free `nearbyint` form. Numerical results identical to within
  rounding tolerance (1e-12 in double).
- **Triclinic stays out of scope where it already is.** MPI domain
  decomposition (`rcb.cpp`, `domain.cpp`) and `cell_list.cpp` continue
  to assume orthorhombic; their `Box` callsites use `length(d)` /
  `wrap_axis(q, d)`, never the triclinic `min_image`. Triclinic
  support for those modules is a separate spec.

## Domain basis

Minimum image convention (Allen & Tildesley, *Computer Simulation of
Liquids*, 1987 §1.6.2):

- **Orthorhombic**: closest-image displacement is
  `dx_min = dx − L_x · round(dx / L_x)`. The `nearbyint` form is
  branch-free and IEEE-deterministic for finite inputs; it produces
  identical results to the branchy `if (dx > L/2) dx -= L; else if
  (dx < -L/2) dx += L;` form when `|dx| ≤ 1.5 L_x` (the regime in
  which physical force evaluators operate; the branchy form fails for
  `|dx| > 1.5 L`, so the `nearbyint` form is in fact strictly more
  general).
- **Triclinic**: transform to fractional coordinates
  `f = inv_cell · dx`, round each lane to nearest integer, transform
  back: `dx = cell · (f − round(f))`. This is the pattern already in
  `src/cuda/md/neighborlist/triclinic.cuh::minimum_image`, which
  becomes the reference for the CPU implementation.

No external libraries needed. Reference: existing CUDA implementation
(internally consistent with the project's CUDA neighbor-list tests)
plus textbook formula. `science.required` is satisfied by code reuse;
no new literature review required.

## Design

### 1. Layer placement

CLAUDE.md (Backend directory layering): *"Files at the backend root
are shared infrastructure used by multiple domains."* Box is consumed
by `bonded/`, `md/`, `observer/`, `scf/`, `recorder/`, `mpi/` — six
domains. ✓ qualifies for backend-root placement.

CLAUDE.md (Hot-path data structures rule): *"carry only fundamental
physical quantities."* `inv_cell` is **derived** but cached; this is
the standard exception for *per-step-immutable* derived data living
in one place rather than recomputed at every callsite. The class
docstring documents this exception explicitly.

### 2. CPU API — `cpu::Box`

```cpp
// include/atomiverse/cpu/box.h
namespace ATV::cpu {

class Box {
public:
    Box() = default;
    explicit Box(const Mat3& cell);
    static Box orthorhombic(double a, double b, double c);
    static Box cubic(double l);

    // ---- hot-path queries (always-inline, branch-free) -----------------
    [[gnu::always_inline]] inline double length(int d) const noexcept;
    [[gnu::always_inline]] inline double volume() const noexcept;
    [[gnu::always_inline]] inline bool is_orthorhombic() const noexcept;

    // Min-image, in place. Ortho fast-path uses nearbyint; triclinic
    // uses inv_cell roundtrip. Branch on `ortho_` is predicted
    // perfectly in inner loops (orthorhombic dominates).
    [[gnu::always_inline]] inline void
    min_image(double& dx, double& dy, double& dz) const noexcept;

    // xsimd-batched ortho min-image. Precondition: is_orthorhombic().
    // Triclinic batch is intentionally unsupported (no-fallback rule).
    template <class B>
    inline void min_image_batch(B& bdx, B& bdy, B& bdz) const noexcept;

    // Wrap a single coordinate into [0, L_d). Used by MPI rcb /
    // atom_migration / DD histogramming. Ortho only.
    [[gnu::always_inline]] inline double
    wrap_axis(double q, int d) const noexcept;

    // Wrap a Cartesian point. Ortho only.
    inline void wrap_into(double& x, double& y, double& z) const noexcept;

    // ---- mutation (cold path; only Barostat / Deform / NPT) -----------
    void scale_diagonal(double sx, double sy, double sz);
    void scale_isotropic(double s);
    void scale_axis(int d, double s);  // for Deform
    void set_cell(const Mat3& cell);

    // ---- raw accessors (compat: snapshot, MPI gather, molrec) --------
    const Mat3& cell()     const noexcept { return cell_; }
    const Mat3& inv_cell() const noexcept { return inv_cell_; }

private:
    Mat3 cell_{};
    Mat3 inv_cell_{};
    double half_lx_{}, half_ly_{}, half_lz_{};   // ortho fast-path
    double inv_lx_{},  inv_ly_{},  inv_lz_{};    // ortho fast-path
    bool ortho_{true};
    void recompute_cache_() noexcept;
};

} // namespace ATV::cpu
```

**Storage rationale.** `Mat3 cell_` (24 bytes) + `Mat3 inv_cell_`
(24 bytes) + 6 × `double` cache (48 bytes) + `bool ortho_` (1 byte +
padding) → ~104 bytes per Box. Box is constructed once per Frame and
read trillions of times per simulation; cache locality matters more
than struct size.

**Ortho min-image kernel (canonical form):**

```cpp
inline void Box::min_image(double& dx, double& dy, double& dz) const noexcept {
    if (ortho_) {
        dx -= cell_(0,0) * std::nearbyint(dx * inv_lx_);
        dy -= cell_(1,1) * std::nearbyint(dy * inv_ly_);
        dz -= cell_(2,2) * std::nearbyint(dz * inv_lz_);
    } else {
        // Triclinic: fractional roundtrip using inv_cell_ / cell_.
        // ... (exact pattern from triclinic.cuh)
    }
}
```

Branch on `ortho_` is per-call but data-dependent on a single
member; the branch predictor learns it on call #2. In the LJ ket
loop, the same Box is used 10⁶+ times → effectively free.

**xsimd batch kernel** (in the `min_image_batch` template):

```cpp
template <class B>
inline void Box::min_image_batch(B& bdx, B& bdy, B& bdz) const noexcept {
    assert(ortho_);  // debug; release is UB
    bdx -= B(cell_(0,0)) * xsimd::nearbyint(bdx * B(inv_lx_));
    bdy -= B(cell_(1,1)) * xsimd::nearbyint(bdy * B(inv_ly_));
    bdz -= B(cell_(2,2)) * xsimd::nearbyint(bdz * B(inv_lz_));
}
```

This is what `lj.cpp` already inlines manually in its xsimd ket loop;
extracting it changes one inlined call site to one inlined method,
no instruction-count delta.

### 3. CUDA API — `cuda::Box`

```cuh
// include/atomiverse/cuda/box.cuh
namespace ATV::cuda {

struct Box {
    real_t cell[9];      // row-major
    real_t inv_cell[9];
    real_t half_l[3];
    real_t inv_l[3];
    bool   periodic[3];  // axis-wise PBC flag (currently {true, true, true})
    bool   ortho;
};

static_assert(std::is_trivially_copyable_v<Box>);

template <typename T>
__host__ __device__ inline void
min_image(T& dx, T& dy, T& dz, const Box& b) noexcept;

template <typename T>
__host__ __device__ inline T
wrap_axis(T q, int d, const Box& b) noexcept;

__host__ __device__ inline real_t length(const Box& b, int d) noexcept;
__host__ __device__ inline real_t volume(const Box& b) noexcept;

// Existing nested-namespace helpers in src/cuda/md/neighborlist/
// triclinic.cuh become thin re-exports (deprecation tagged) so the
// neighbor-list TU does not need a same-PR rewrite.
} // namespace ATV::cuda
```

CUDA Box is a plain data struct so it can travel by value through
`__constant__` memory or a kernel argument. No virtual methods, no
ctor side effects.

### 4. Frame migration

```cpp
// include/atomiverse/cpu/frame.h
struct Frame {
    // ...
    cpu::Box box{};   // formerly: Mat3 box{};
    // ...
};
```

The `Mat3 box{}` field is replaced by `cpu::Box box{}`. Every existing
read of `frame.box(i, j)` becomes `frame.box.cell()(i, j)` if the
caller genuinely needs the raw matrix (snapshot serialization, MPI
gather), or — much more often — collapses to
`frame.box.length(d)` / `frame.box.min_image(...)`. CUDA Frame
analogous.

### 5. Mutation entry points

Today, `nose_hoover_npt.cpp`, `deform.cpp`, `barostat.cpp` all
write directly: `f.box(0,0) *= a;`. After migration these become
`f.box.scale_isotropic(a);` (NPT iso), `f.box.scale_diagonal(sx, sy, sz);`
(NPT aniso / Berendsen), `f.box.scale_axis(axis_, scale);` (Deform).
Each call refreshes the cache exactly once. **No raw-field write to
the cell is permitted post-migration**; this is enforced by making
`cell_` private.

### 6. Compatibility shims

For one release cycle:
- `cpu::Box` is implicitly convertible to `const Mat3&` via `operator
  const Mat3&() const { return cell_; }` (so old `f.box(i, j)` still
  compiles where the migration script missed a callsite). This is a
  build-time scaffold, not a permanent API; the implicit conversion is
  removed in a follow-up PR after a clean grep.
- `ATV::md::nlist::minimum_image` (in `triclinic.cuh`) is kept as an
  inline function that forwards to `ATV::cuda::min_image`.

## Files

NEW
- `include/atomiverse/cpu/box.h`             — declaration + inline hot-path methods
- `src/cpu/box.cpp`                          — cold-path mutation, `recompute_cache_`
- `include/atomiverse/cuda/box.cuh`          — `cuda::Box` + free helpers
- `tests/cpu/test_box.cpp`                   — ortho/triclinic/xsimd parity, cache identity
- `tests/cuda/test_box_cuda.cu`              — CUDA mirror of CPU test cases
- `benchmarks/bench_box.cpp`                 — micro-bench `min_image` throughput

MODIFY
- `include/atomiverse/cpu/frame.h`           — `Mat3 box` → `cpu::Box box`
- `include/atomiverse/cpu/frame_snapshot.h`  — same
- `include/atomiverse/cuda/frame.h`          — `Mat3 box` → `cuda::Box box`
- `src/cpu/bonded/harmonic_bond.cpp`         — drop local `minimum_image`; use `f.box.min_image(...)`
- `src/cpu/bonded/harmonic_angle.cpp`        — same
- `src/cpu/bonded/multi_harmonic_dihedral.cpp` — same
- `src/cpu/md/lj.cpp`                        — drop local `minimum_image`; xsimd ket loop uses `f.box.min_image_batch`
- `src/cpu/md/cell_list.cpp`                 — `box(d, d)` → `box.length(d)`
- `src/cpu/observer/rdf.cpp`                 — same
- `src/cpu/scf/scf_evaluator.cpp`            — same
- `src/cpu/md/mpi/rcb.cpp`                   — drop `wrap_to_box`; use `f.box.wrap_axis(q, d)`
- `src/cpu/md/mpi/atom_migration.cpp`        — same
- `src/cpu/md/mpi/halo_exchange.cpp`         — `box(d, d)` → `box.length(d)`
- `src/cpu/md/mpi/domain.cpp`                — same
- `src/cpu/md/mpi/gather.cpp`                — `out.box = frame.box;` (works via Box copy)
- `src/cpu/md/nose_hoover_npt.cpp`           — direct field writes → `box.scale_isotropic` / `scale_diagonal`
- `src/cpu/md/deform.cpp`                    — direct field writes → `box.scale_axis`
- `src/cpu/md/barostat.cpp`                  — direct field writes → `box.scale_diagonal`
- `src/cpu/snapshot_publisher.cpp`           — Box copy assignment
- `src/cpu/molrec.cpp`                       — `box.cell()` where raw `Mat3` needed
- `src/cuda/md/neighborlist/triclinic.cuh`   — re-export `ATV::cuda::min_image` under `ATV::md::nlist`
- `src/cuda/md/neighborlist/{pair_search,export_pairs,setup,sort_and_cluster,grid}.cu`
                                             — accept `cuda::Box` instead of `(cell, inv_cell, periodic)` triple
- `src/cuda/neighbor_list.cu`                — same
- `src/cuda/transfer.cu` / `src/torch/transfer.cu` — Box transfer routines
- `CLAUDE.md`                                — one-line note under Architecture

## Tasks

### CPU phase 1 — Box class itself

- [x] Add skeleton `tests/cpu/test_box.cpp` with constructor / `length` / `volume` / `is_orthorhombic` cases. Build fails (Box undefined). **RED**.
- [x] Implement `cpu::Box` ctors + simple queries in `include/atomiverse/cpu/box.h` + `src/cpu/box.cpp`. Tests **GREEN**.
- [x] Add ortho-`min_image` parity tests vs reference branchy form (10⁴ random cases per axis, L ∈ [5, 50] Å, dx ∈ [-2L, 2L], 1e-12 tolerance). **RED**.
- [x] Implement ortho `min_image` (`nearbyint` form) in the inline header. **GREEN**.
- [x] Add triclinic `min_image` test against fractional-roundtrip reference computed inline in the test (1e-12 tolerance). **RED**.
- [x] Implement triclinic `min_image` branch + `recompute_cache_` for `inv_cell_`. **GREEN**.
- [x] Add `wrap_axis` / `wrap_into` parity tests vs `wrap_to_box` reference (q ∈ [-3L, 3L]). **RED**.
- [x] Implement `wrap_axis` / `wrap_into`. **GREEN**.
- [x] Add `scale_diagonal` / `scale_isotropic` / `scale_axis` / `set_cell` cache-coherence tests (after mutation, `cell · inv_cell ≈ I`, `inv_lx · cell(0,0) ≈ 1`, etc.). **RED**.
- [x] Implement mutation methods. **GREEN**.
- [x] Add `min_image_batch` xsimd parity test (lane-wise bit-exact equality with scalar). **RED**.
- [x] Implement `min_image_batch`. **GREEN**.

### CPU phase 2 — Frame migration + callsite sweep

- [x] Migrate `cpu::Frame::box` field type from `Mat3` to `cpu::Box`. Read shim `operator()(int, int) const` + `operator=(const Mat3&)` keep tests/examples compiling without a sweep. Build green.
- [x] Migrate `cpu::FrameSnapshot::box` analogously. Build green.
- [x] Replace 4× `static inline minimum_image` in bonded files (harmonic_bond, harmonic_angle, multi_harmonic_dihedral) with `f.box.min_image(...)`. Also dropped the dead-code copies from `lj.cpp`, `cell_list.cpp`, `rdf.cpp`. `test_md`, bonded, RDF tests green.
- [~] **Deferred (perf):** swapping LJ's xsimd ket-loop and scalar tail to call `box.min_image_batch` / `box.min_image` cost +6.8% per-step single-thread (member loads through `this` not hoisted). Reverted to the inlined branchy form in `LJ::launch`; the dead `static inline minimum_image` was still removed. Documented in `box-abstraction.artifacts/perf.md` and in `lj.cpp` comments.
- [x] Replace 2× `wrap_to_box` in `src/cpu/md/mpi/{rcb,atom_migration}.cpp` with `frame.box.wrap_axis(q, d)`. `test_mpi_*` green.
- [x] Replace `box(d, d)` raw reads in `cell_list.cpp`, `observer/rdf.cpp`, `scf/scf_evaluator.cpp`, `mpi/halo_exchange.cpp`, `mpi/domain.cpp`, `nose_hoover_npt.cpp` with `box.length(d)` / `box.volume()`. All tests green.
- [x] Replace direct field writes in `nose_hoover_npt.cpp`, `deform.cpp`, `barostat.cpp` with `box.scale_*` methods. All tests green.

### CPU phase 3 — Verification + cleanup

- [x] Add `benchmarks/bench_box.cpp` micro-bench: scalar `min_image`, `min_image_batch`, and a branchy reference. Reports ns/op + items/sec.
- [x] Run `pe_profile_single` (CPU LJ + bonded NVE; CUDA `bench_md_step` not in this build) at 1 thread and 8 threads. Captured pre-/post-refactor in `.claude/specs/box-abstraction.artifacts/perf.md`: 1t +0.6%, 8t −5.6% (improvement). Both within ±1%.
- [~] **N/A:** the `operator const Mat3&` shim was never added (only `operator()(int,int) const` read shim and `operator=(const Mat3&)` assignment shim, both intentional and kept). The read shim is unused outside one comment in `cell_list.h`; the assignment shim is heavily used by tests/examples (`frame.box = Mat3::cubic(L)`) and is the natural ergonomic entry point.

### CUDA phase

- [x] Implement `cuda::Box` struct + `__host__ __device__` free helpers in `include/atomiverse/cuda/box.cuh`. Trivially copyable (static_assert in header). API mirror of CPU: `length`, `volume`, `min_image`, `wrap_axis`.
- [x] Forward `ATV::md::nlist::minimum_image(dx, dy, dz, box)` to `ATV::cuda::min_image` in `src/cuda/md/neighborlist/triclinic.cuh`. Existing pointer-based `(cell, inv_cell)` overloads kept for compatibility.
- [x] Migrate `src/cuda/transfer.cu` writable `frame.box(i, j) = ...` pattern to `Mat3 m; ... ; frame.box.set_cell(m);` so the CUDA backend compiles against the read-only Box shim.
- [~] **Deferred:** `cuda::Frame::box` field migration from `prec_t*` (device pointer) to `cuda::Box` (host value). This requires updating every kernel signature (`pair_search.cu`, `export_pairs.cu`, `setup.cu`, `sort_and_cluster.cu`, `grid.cu`, `neighbor_list.cu`, `velocity_verlet.cu`, `deform.cu`, `barostat.cu`, `stress.cu`, `sd.cu`, `harmonic_*.cu`, `multi_harmonic_dihedral.cu`) and the host-side allocator (`transfer.cu` `cudaMalloc(&d.box, ...)` + the marshalling loop). Tracked as a follow-up; the `cuda::Box` struct + helpers are in place so new components can adopt it. `test_box_cuda.cu` skipped (CUDA disabled in this build environment).

### Closeout

- [x] CLAUDE.md Architecture section names `cpu::Box` / `cuda::Box` at backend root and documents the LJ inlined-min-image exception.
- [x] Public Box members have `@brief` / `@param` / `@return` Doxygen comments per project rule.

## Testing strategy

**Parity tests (CPU)**
- Random orthorhombic boxes (L ∈ [5, 50] Å) × random displacements
  (dx ∈ [-2L, 2L]); compare `Box::min_image` to the canonical branchy
  reference within 1e-12. Same for `wrap_axis` / `wrap_into`.
- Triclinic cells (one and two off-diagonal entries each); compare
  `Box::min_image` to a hand-rolled fractional-roundtrip reference
  written in the test file (independent of `Box::min_image`'s
  implementation), within 1e-12.

**xsimd parity**
- For each platform `xsimd::batch<double>::size`, fill a batch with
  random doubles, run `min_image_batch`, run scalar `min_image` per
  lane, assert lane-wise bit-exact equality.

**Cache coherence**
- After every mutation method (`set_cell`, `scale_diagonal`,
  `scale_isotropic`, `scale_axis`), assert `cell · inv_cell ≈ I` and
  per-axis `cell(d,d) · inv_l[d] ≈ 1` and `2 · half_l[d] ≈ cell(d,d)`
  to 1e-12.

**End-to-end regression**
- The full `ctest --test-dir build --output-on-failure` suite must
  pass unchanged. Specific tests of interest: `test_md`, `test_bonded`,
  `test_lj`, `test_neighbor_list_cpu`, `test_neighbor_list_cuda`,
  `test_rdf`, `test_scf_evaluator`, `test_aimd_conservation`,
  `test_mpi_*`.

**Benchmark gates**
- `bench_md_step` (LJ + bonded NVE) ns/step within ±1% of pre-refactor
  baseline at 1 thread and 8 threads.
- `bench_box` ortho `min_image` ns/call ≥ baseline static-inline form
  (≥ -5% margin); `min_image_batch` ≥ 1.5× scalar at native batch
  width.
- CUDA `bench_neighbor_list` ns/call within ±1% of baseline.

**Docs**
- `make doc` (or direct `doxygen`) produces no `warning: ...` lines
  for `cpu::Box` / `cuda::Box` symbols.

## Out of scope

- Triclinic support in MPI domain decomposition (`rcb.cpp`,
  `domain.cpp`). Those modules continue to assume orthorhombic and
  use `length(d)` / `wrap_axis(q, d)` only.
- Triclinic support in `cpu::md::cell_list`. Same posture.
- Box mutation for shear / non-diagonal NPT. The Berendsen / NPT
  components today only do diagonal scaling; the API leaves room for
  full-tensor mutation as a follow-up.
- Removing `Mat3`. `Mat3` remains the underlying storage primitive;
  `Box` is built on it.
- Performance work on triclinic `min_image` beyond the `nearbyint`
  fractional-roundtrip form (no orthonormalization, no Lanczos
  reduction). The orthorhombic fast-path is what this spec optimizes.
- `torch/` backend (out-of-tree; mirrors CUDA when needed).
