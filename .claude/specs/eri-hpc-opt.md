---
title: ERI Shell-Block HPC Optimization
slug: eri-hpc-opt
status: in-progress
created: 2026-04-14
revised: 2026-05-19
---

# ERI Shell-Block HPC Optimization

**Priority:** high
**Progress:** phases 1-3 implemented; phase 4+ pending
**Depends on:** none (shell-block `eri_block()` already merged)

## Problem

Single Fock build for aspirin (STO-3G, 73 AOs): **~2460 ms** vs PySCF direct SCF **~952 ms** (2.5× gap).
Full RHF SCF solve: **~38 s** vs PySCF direct **~11.5 s**.

Our algorithm structure already mirrors libcint/PySCF (shell-block grouping, Boys amortization,
Cauchy-Schwarz screening, OpenMP at the quartet level). The gap is purely compute efficiency
inside `eri_block()`.

### Root cause analysis

| # | Issue | Impact |
|---|-------|--------|
| 1 | No `-march=native` — baseline scalar, no NEON/AVX | ~1.5× |
| 2 | No xsimd — xtensor reduction in `build_g_rhf_direct_block` is scalar | ~1.3× |
| 3 | Generic `eri_block()`: runtime branches over `l`, `ang_table` lookups, 6-nested loop with variable bounds | ~2× |
| 4 | Boys function called once per quartet but not vectorized across multiple quartets | ~1.3× |
| 5 | R-table: non-contiguous 5×5×5×5 array makes inner contraction loop non-vectorizable | ~1.5× |

## Alignment with PySCF / libcint

libcint approach (what we must replicate):
1. **Shell-type dispatch:** auto-generates one C function per `(l_a, l_b, l_c, l_d)` combination.
   No runtime branches in the hot path. Each specialization has fully unrolled angular loops.
2. **`-march=native` always on**: libcint is always compiled with native tuning; uses NEON on ARM64,
   AVX2/AVX-512 on x86-64. This is the single largest gap item.
3. **xsimd-backed xtensor**: `XTENSOR_USE_XSIMD` activates aligned allocators and SIMD expression
   evaluation for all xtensor ops (Fock reduction, density matrix, DIIS).
4. **Vectorized Boys / prefactor**: processes multiple ket primitive pairs simultaneously using
   SIMD `erf`, `sqrt`, and arithmetic.
5. **BLAS primitive contraction** (for `l ≥ 2`): expresses primitive-quartet contraction as DGEMM.
   Not needed for STO-3G (`l_max=1`), essential for cc-pVDZ (`l_max=2`).

## Strategy

Add xsimd and `-march=native`; add shell-type-specialized `eri_block` variants; update the HPC
agent and skill files to encode these patterns as the canonical approach.

## Phases

### Phase 1 — Build flags & xsimd dependency ✓ DONE

**CMakeLists.txt** actual changes:

1. `xsimd` FetchContent **tag 11.2.0** (NOT 13.0.0 — xsimd 13.x breaks xtensor 0.25.0's `reshape_view`).
2. **Do NOT** set `XTENSOR_USE_XSIMD` as a cmake cache variable or compile definition globally.
   Enabling it activates SIMD-backed xtensor storage, which makes `const_reference` a SIMD batch
   type — then `reshape_view` fails to compile (xstrided_view tries to take address of rvalue).
   Use xsimd only directly via `#include <xsimd/xsimd.hpp>` in .cpp files.
3. `-march=native` via language-guarded generator expression (language guard prevents NVCC poisoning):
   ```cmake
   target_compile_options(atomiverse PRIVATE
       $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
   ```
   Note: `-fopt-info-vec-missed` is GCC-only; Clang rejects it. Omit it.
4. Pass `-DXTENSOR_USE_XSIMD=OFF` at cmake invocation to prevent xtensor from calling
   `find_package(xsimd)` from a stale cache entry.

### Phase 2 — Compile-time Ltot dispatch ✓ DONE

**`src/cpu/scf/native_integrals.cpp`**

**Dispatch by Ltot** (NOT by bit-pattern `(A.l<<3)|...` which misses permutations like sp|ss):
```cpp
void eri_block(const ShellBlock& A, ..., double* out) {
    const int Ltot = A.l + B.l + C.l + D.l;
    if (Ltot == 0) return eri_blk_ss_ss(A, B, C, D, out);  // xsimd
    switch (Ltot) {
        case 1: return eri_blk_ltot<1>(A, B, C, D, out);
        case 2: return eri_blk_ltot<2>(A, B, C, D, out);
        case 3: return eri_blk_ltot<3>(A, B, C, D, out);
        case 4: return eri_blk_ltot<4>(A, B, C, D, out);
        default: return eri_blk_general(A, B, C, D, out);  // l≥2
    }
}
```

**Key lesson**: dispatching on `(A.l<<3)|(B.l<<2)|(C.l<<1)|D.l` only handles specific orderings
(e.g., ps|ss=0b1000) but NOT permutations (sp|ss=0b0100 falls to default with 9×9×9×9 table).
Ltot-based dispatch handles ALL orderings correctly.

`eri_blk_ltot<LTOT>` is a template identical to the general code but with `Ltot_max=LTOT` as a
compile-time constant, enabling the compiler to unroll Boys recurrence and R-table build loops.
`eri_blk_general` uses `double R_tab[9][9][9][9]` (for l≥2 future support).

### Phase 3 — xsimd SIMD vectorization of ss|ss ✓ DONE

**`src/cpu/scf/native_integrals.cpp`** — `eri_blk_ss_ss`

Precompute ket geometry before bra loop; for each bra pair vectorize ket loop:
```cpp
batch_t acc_v(0.0);
for (int kk = 0; kk + W <= nket; kk += W) {
    batch_t T_v    = batch_t::load_unaligned(T_k + kk);
    batch_t pref_v = batch_t::load_unaligned(pref_k + kk);
    auto small  = T_v < batch_t(1e-12);
    auto safe_T = xsimd::select(small, batch_t(1.0), T_v);
    auto F0_v   = xsimd::select(small, batch_t(1.0),
                      batch_t(0.5)*xsimd::sqrt(batch_t(PI)/safe_T)*xsimd::erf(xsimd::sqrt(safe_T)));
    acc_v += pref_v * F0_v;
}
// Horizontal reduce: store-to-array (xsimd::hadd NOT available in xsimd 11.2.0 as free function)
double arr[16]; acc_v.store_unaligned(arr); double acc=0; for(int i=0;i<W;i++) acc+=arr[i];
```

**Key lessons learned:**
- `xsimd::hadd` is NOT a free function in xsimd 11.2.0. Use `store_unaligned` + scalar sum.
- Use `safe_T = select(small, 1.0, T)` before sqrt/erf to guard against division by zero.
- Precompute ket K_CD, Q, gamma2 before the bra loop to avoid redundant `exp()` calls per bra pair.

For higher-l cases, `-march=native -O3` enables compiler auto-vectorization of the scalar loops.

### Phase 4 — BLAS angular contraction (future, l ≥ 2)

Applicable when supporting d-shells (cc-pVDZ, 6 Cartesian functions per shell).
For (dd|dd): 6×6×6×6 = 1296 angular combinations. The contraction becomes DGEMM:
```
out[na*nb, nc*nd] += E_bra[na*nb, n_hermite] * R_mat[n_hermite, n_hermite] * E_ket^T
```
Use `xt::linalg::dot()` (already available via xtensor-blas). Not needed for STO-3G.

### Phase 5 — Update av-hpc agent and skill

**`.claude/agents/av-hpc.md`** — add to SIMD section:
- xsimd API: `batch_t::load_aligned`, `xsimd::erf`, `xsimd::select`, `xsimd::hadd`
- `-march=native` CMake pattern with `$<$<COMPILE_LANGUAGE:CXX>:...>` language guard
- `XTENSOR_USE_XSIMD` compile definition pattern
- ERI shell-type dispatch as the canonical HPC ERI pattern
- Boys SIMD vectorization: `xsimd::erf(xsimd::sqrt(T_v))` with select for T<1e-12

**`.claude/skills/av-hpc/SKILL.md`** — add in Implementation mode:
- ERI `eri_block()` follows shell-type dispatch; specializations in `native_integrals.cpp`
- `-march=native` is always required for the CPU SCF backend
- xsimd ket-loop pattern for primitive quartet vectorization

**`.claude/skills/av-perf/SKILL.md`** — add profiling guidance:
- ERI profiling: time `build_g_rhf_direct_block` vs `build_g_rhf_direct` with `bench_native`
- Vectorization blockers in `eri_block()`: R-table non-contiguous access, ang_table branch

## Iron rules (added to CLAUDE.md and agent)

1. **`-march=native` is always on** for the CPU SCF backend. Never build without it.
2. **Do NOT set `XTENSOR_USE_XSIMD` globally** — it breaks `reshape_view` in `xc_grid_engine.h`.
   Use xsimd directly via `#include <xsimd/xsimd.hpp>` in hot .cpp files.
3. **xsimd version**: use 11.2.0. xsimd 13.x is incompatible with xtensor 0.25.0.
4. **Ltot-based dispatch in `eri_block()`** — dispatch on `A.l+B.l+C.l+D.l`, not bit-pattern.
   Bit-pattern dispatch misses permutations (sp|ss → general fallback → 9^4 table).
5. **Correctness gate**: all BlockG tests must pass to 1e-9 (`NativeIntegrals.BlockG_*`).
6. **Regression gate**: `NativeVsLibcint` tolerance is 1e-7; may tighten with native improvements.

## Test plan

```bash
# After each phase, run:
cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_TESTS=ON \
      -DATV_BUILD_BENCHMARKS=ON -DATV_BUILD_LIBCINT=ON
cmake --build build -j$(nproc)

# Correctness
ctest --test-dir build -R "BlockG|NativeIntegrals" --output-on-failure

# Performance (compare with PySCF baseline ~952 ms direct Fock build)
./build/benchmarks/bench_native \
    --benchmark_filter="BM_GBuild_Block_Aspirin|BM_GBuild_AO_Aspirin" \
    --benchmark_repetitions=5
```

Expected benchmark progression (wall time, Fock build, aspirin STO-3G):

| Phase | Expected | vs PySCF |
|-------|----------|----------|
| Baseline | 2460 ms | 2.6× slower |
| Phase 1 (-march=native + xsimd alloc) | ~1500 ms | 1.6× slower |
| Phase 2 (shell-type specialization) | ~800 ms | ~0.85× (parity+) |
| Phase 3 (xsimd ket loop) | ~500 ms | **faster** |

## Files

| File | Change |
|------|--------|
| `CMakeLists.txt` | xsimd FetchContent; `-march=native`; `XTENSOR_USE_XSIMD` |
| `src/cpu/scf/native_integrals.cpp` | Phase 2: 6 specializations + dispatcher; Phase 3: xsimd |
| `include/atomiverse/cpu/scf/native_integrals.h` | Declare `eri_blk_*` static helpers (internal only) |
| `.claude/agents/av-hpc.md` | xsimd API, CMake pattern, ERI dispatch knowledge |
| `.claude/skills/av-hpc/SKILL.md` | ERI pattern in Implementation mode |
| `.claude/skills/av-perf/SKILL.md` | ERI profiling guidance |
| `benchmarks/bench_native.cpp` | Already has `BM_GBuild_Block_Aspirin` (added in this session) |

## Tasks

- [x] Add xsimd 11.2.0 FetchContent and `-march=native` language-guarded flag in `CMakeLists.txt`
- [x] Implement Ltot-based dispatcher and `eri_blk_ltot<LTOT>` template in `src/cpu/scf/native_integrals.cpp`
- [x] Implement xsimd-vectorized `eri_blk_ss_ss` ket loop with safe-T guard in `src/cpu/scf/native_integrals.cpp`
- [x] Verify `BlockG|NativeIntegrals` ctest filters pass to 1e-9 / 1e-7 tolerance
- [x] Benchmark `BM_GBuild_Block_Aspirin` reaches Phase 3 target (~500 ms, faster than PySCF)
- [ ] Phase 4: Implement BLAS angular contraction path for `l ≥ 2` via `xt::linalg::dot` in `src/cpu/scf/native_integrals.cpp` (cc-pVDZ enablement)
- [ ] Phase 4: Add cc-pVDZ correctness test against libcint to 1e-7 in `tests/scf/test_native_integrals.cpp`
- [ ] Phase 4: Add `BM_GBuild_Block_Aspirin_ccPVDZ` benchmark in `benchmarks/bench_native.cpp` and record PySCF baseline
- [ ] Phase 5: Update `.claude/agents/av-hpc.md`, `.claude/skills/av-hpc/SKILL.md`, `.claude/skills/av-perf/SKILL.md` with xsimd / Ltot-dispatch / Boys-SIMD patterns
- [ ] Run full check + test suite (`ctest --test-dir build --output-on-failure`)
