# Precision Migration — Atom/Frame to tiered types

status: Phase 1 + 2 + 3 (CPU) done; CUDA tier unification pending (GPU node)
date: 2026-07-01

## Phase 3 outcome (integers, done)

`cpu::Atom` `tag` (and the bonded partner arrays `bond_atom`/`angle_atom`/
`dihedral_atom`/`special`, and the `tag→local` map key) migrated `int64 → TagInt`.
`-DATV_ID_SIZE=smallbig` gives **32-bit IDs** (LAMMPS default + molrs's on-disk tag
width) — compiles + the tag/topology/comm/materialize tests pass; `bigbig` (64-bit,
**the Atomiverse default**) is byte-identical to pre-migration.

- **Default flipped to `bigbig`** so the default build preserves the legacy always-
  int64 tag width (byte-identical), with `smallbig` as the memory-saving opt-in —
  the same "wide default, narrow opt-in" shape as `float` storage. CI (which does
  not set `ATV_ID_SIZE`) therefore stays on int64 tags, unchanged.
- Same self-consistent recipe, no casts: `const auto*` tag pointers, templated
  `copy_range_i64`, generic `compact_1d_i64` (`auto&`), `TagInt` field/param/map
  types, one test static_assert made tier-aware (`is_same<tag_t, TagInt>`).
- `type`/counts left as `int` (== `LocalInt`, byte-identical); `image` left as the
  per-axis `int8` (not the LAMMPS packed `imageint`); global step counters were
  already `int64_t` (BigInt-equivalent). CUDA `real_t/work_t/prec_t` + int tiers
  unification is the only remaining piece — code-only, needs the GPU node.
depends on: [[project-precision-tiered-lammps]] (core/precision.h foundation, committed)

## Outcome (Phase 2 done)

`cpu::Atom` storage fields (x/y/z, v*, mass, charge) are now `StorageReal`; forces
are `AccReal`. **`-DATV_REAL_STORAGE=float` compiles and runs** classical MD with
physically-sound results at float precision (ε~1e-7). Per-atom real storage for
positions+velocities drops ~2×.

**The type system stayed self-consistent — no scattered casts.** The coherent
recipe that worked:
- Raw pointers deduce the tier: `const auto* x = atoms.x.data()` (not `const
  double*`), so `x[i]` promotes to double in the math automatically.
- Helpers are templated on the element type (`copy_range`, `copy_optional_f64`)
  or generic (`gather` widen-to-double lambda, `compact_1d` generic lambda);
  bodies use plain assignment — C++'s implicit float↔double conversion does the
  rest. Zero `static_cast`.
- Compute/comm/output buffers stay `double`; conversion happens implicitly at the
  storage-array read/write boundary only.

**Supported precision matrix** (mirrors the CI rows):
- `double` storage + `RHF=ON` — full suite, SCF/DFT.
- `float` storage + `RHF=OFF` — classical MD, memory-optimized. **SCF requires
  double** (Gaussian-center integral accuracy; float+RHF is out of scope — its
  only build blocker is `{double}` brace-init narrowing in SCF tests, left as-is).

**Precision-exact tests need double.** Under float, ~4 classical tests that assert
double-tight bounds (stationary-point forces < 1e-10, analytic energy to 1e-9,
NVE-limit exact `<=`) legitimately fail — they detect the reduced precision, which
is correct behavior, not a bug. CI's float row excludes them (they are double-only
correctness checks); a finer per-subtest `GTEST_SKIP when StorageReal != double`
guard is a follow-up.

## Goal

Make the tiered precision from `core/precision.h` actually *bite*: migrate the
`cpu::Atom` SoA fields and the global counts from hardcoded `double`/`int` to the
role-correct tiers, so `-DATV_REAL_STORAGE=float` genuinely halves per-atom
storage and `-DATV_ID_SIZE=bigbig` gives 64-bit IDs — while the default build
stays byte-identical (all tiers default to the wide type).

## Field → tier plan (cpu/frame/atom.h)

| field(s) | current | target tier | rationale |
|---|---|---|---|
| `fx/fy/fz` forces | `double` | **`AccReal`** (double, fixed) | force kernels do `fx[i] += …` — the array **is** the accumulator; never economize (LAMMPS double-accumulation rule) |
| `x/y/z` positions | `double` | `StorageReal` | per-atom storage; the biggest arrays. **SIMD hazard** (below) |
| `vx/vy/vz` velocities | `double` | `StorageReal` | per-atom storage |
| `mass`, `charge` | `double` | `StorageReal` | read-mostly per-atom storage |
| `type` | `int` | `LocalInt` (int32) | local type index |
| `tag` | `int64_t` | `TagInt` | **nuance**: `tag` is currently *always* int64; `TagInt` is int32 under `smallbig` (default) → this would NARROW tag to 32-bit unless `bigbig`. The FFI, `bond_atom`, `map` all use int64 tags — migrate carefully or default `ATV_ID_SIZE=bigbig` |
| `image` | `int8_t` | keep / `ImageInt`? | image is int8 per-axis here, not LAMMPS packed imageint — leave unless a packed form is adopted |

## The hazard — why this is not a one-liner (audit 2026-07-01)

Every real storage field has `double*`-typed consumers that break under `float`:
- **SIMD positions**: `pair_lj_cut.cpp` (22 `x.data()` sites), `box_simd.h`,
  `native_integrals.cpp` load `batch<double>` from `x.data()`. Under `float` this
  needs load-and-convert, and NEON/AVX have **different lane counts** for
  `batch<float>` vs `batch<double>` — the inner loop must be restructured or a
  scalar fallback selected. Highest-risk, highest-value.
- **`gather.cpp`**: `gather_double(frame.atoms.vx.data(), …)` and `out.vx.assign(
  data, data+n)` take/compare `const double*` — needs a `StorageReal` overload or
  templatization.
- **`halo_exchange.cpp:301`**: `const double* q = frame.atoms.charge.data()`.
- **`snapshot_publisher.cpp`**: `copy_range(snap.vx, atoms.vx.data(), …)` — `snap.*`
  are `std::vector<double>` (the async snapshot is deliberately double); the copy
  must convert.

## Phased plan

- **Phase 1 (done, safe): forces → `AccReal`.** `AccReal=double` always → zero
  behavior change, default and float builds both compile (no field uses
  `StorageReal` yet, so the CI `-DATV_REAL_STORAGE=float` row stays green). Makes
  the accumulator tier explicit in the hot data structure.
- **Phase 2 (the payoff): storage fields → `StorageReal`.** Coordinated change:
  (a) templatize/overload `gather_double` + `copy_range` + the halo `q` on the
  field type; (b) restructure the SIMD position loads to load `batch<StorageReal>`
  and convert to `batch<WorkReal>` (or scalar-fallback when `StorageReal != double`);
  (c) keep the async snapshot double (convert on copy). **Validation: an
  energy-conservation test** — NVE LJ over N steps must conserve energy to a
  documented drift bound under `StorageReal=float`, and be byte-identical under
  `double`. Only after this is `-DATV_REAL_STORAGE=float` a real, tested build.
- **Phase 3: integers.** `tag → TagInt` (resolve the int64→int32-under-smallbig
  nuance first — likely default `ATV_ID_SIZE=bigbig` or make the FFI/bond arrays
  tier-aware), local counts `→ LocalInt`, and audit global totals (natoms across
  ranks, ntimestep) for `→ BigInt`. Mirror the CUDA `real_t/work_t/prec_t`
  (`cuda/frame/types.cuh`) to alias the core tiers — validate on the GPU node.

## Acceptance

- Default (`StorageReal=double`) build: byte-identical, full ctest green.
- After Phase 2: `-DATV_REAL_STORAGE=float` compiles, runs, and passes the NVE
  energy-conservation test within the documented drift bound; per-atom memory for
  x/y/z/v drops ~2×.
- No analysis math or accumulation moves to single precision (AccReal stays double).
