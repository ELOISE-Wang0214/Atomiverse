---
name: av-hpc
description: CPU HPC review and optimization -- OpenMP, MPI, SIMD, cache optimization for the CPU backend.
argument-hint: "[path or module]"
user-invocable: true
---

# CPU HPC Review & Optimization

Delegate to the `av-hpc` agent with the user's argument (or all CPU backend files if none given).

The agent operates in two modes:

## Review mode (default)
1. Discover CPU backend files (`src/cpu/`, `include/atomiverse/cpu/`)
2. Identify hot paths (force computation, neighbor list, integrators)
3. Audit parallelism: OpenMP usage, race conditions, reduction patterns
4. Audit vectorization: xtensor ops vs scalar loops, SoA access, SIMD blockers
5. Audit memory: allocation timing, cache friendliness, false sharing
6. Report findings in severity-sorted `[SEVERITY] file:line` format

## Implementation mode (when asked to optimize or implement)
1. Profile before optimizing
2. Add OpenMP / SIMD / MPI as appropriate
3. Validate numerical correctness against serial baseline
4. Report speedup measurements
