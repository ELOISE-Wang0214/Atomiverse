---
name: av-perf
description: Performance review for Atomiverse code.
argument-hint: "[path or module]"
user-invocable: true
---

# Performance Review

Delegate to the `av-optimizer` agent with the user's argument (or all source files if none given).

The agent will:

## CUDA backend (`src/cuda/`, `include/atomiverse/cuda/`)
1. Read CLAUDE.md for GPU data residency and lifecycle rules
2. Scan hot-path code for anti-patterns (D2H in launch, malloc in launch, etc.)
3. Check grid/block sizing, memory access patterns, atomic contention

## CPU backend (`src/cpu/`, `include/atomiverse/cpu/`)
1. Check xtensor usage patterns (unnecessary copies, missing move semantics, redundant evaluations)
2. Look for cache-unfriendly access patterns in SoA data
3. Identify vectorization blockers and missing SIMD-friendly idioms
4. Check for unnecessary allocations in hot loops
5. **ERI performance baseline:** benchmark `BM_GBuild_Block_Aspirin` (shell-block) vs `BM_GBuild_AO_Aspirin` (AO-level) in `bench_native`. Target: Fock build wall time ≤ PySCF direct SCF (~950 ms for aspirin STO-3G). Check that `-march=native` and `XTENSOR_USE_XSIMD` are set.

## Cross-cutting
1. Review transfer functions (`to_cuda`, `to_cpu`) for redundant copies or missing async overlap
2. Check build-time vs launch-time work partitioning (expensive setup should be in build, not launch)
3. Report findings with impact estimates in `[SEVERITY] file:line` format
