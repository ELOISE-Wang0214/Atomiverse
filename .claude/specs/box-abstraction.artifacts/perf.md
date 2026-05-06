# box-abstraction performance gate

Hardware: Apple silicon, 10 cores. PE polymer benchmark
(`benchmarks/pe_profile_single --steps 1000`, 1152 atoms, LJ + bonded).
Three consecutive runs per cell; reported as median.

## pe_profile_single (LJ + bonded NVE)

| Threads | Pre-refactor (ms/step) | Post-refactor (ms/step) | Δ      |
|---------|------------------------|--------------------------|--------|
| 1       | 0.661                  | 0.665                    | +0.6%  |
| 8       | 0.844                  | 0.797                    | −5.6%  |

Both within the ±1% gate (8-thread is a net improvement).

### Raw timestamps (pre)

```
1t: 661, 659, 667  ms/run    median 661 ms / 1000 steps = 0.661 ms/step
8t: 844, 845, 841  ms/run    median 844 ms / 1000 steps = 0.844 ms/step
```

### Raw timestamps (post)

```
1t: 665, 664, 665  ms/run    median 665 ms / 1000 steps = 0.665 ms/step
8t: 797, 777, 797  ms/run    median 797 ms / 1000 steps = 0.797 ms/step
```

## bench_box

Apple silicon, NEON `xsimd::batch<double>::size = 2`,
`benchmark_min_time=1s`. Median of one run (microbench is steady).

| Kernel                        | ns/op  | M ops/s |
|-------------------------------|--------|---------|
| `BM_BoxMinImageScalar`        | 7.31   | 137     |
| `BM_BoxMinImageBatch` (W=2)   | 0.94   | 1066    |
| `BM_BranchyReference`         | 5.66   | 177     |

`BM_BoxMinImageBatch` ≥ 1.5× scalar: **7.78×** (gate pass).
`BM_BoxMinImageScalar` vs reference: −22% (member-load through `this`
not hoisted as cleanly as raw locals; flagged in lj.cpp comments —
LJ keeps the inlined branchy form for the inner loop).

## Decision: LJ keeps inlined min-image

Substituting the LJ scalar tail and ket loop with `Box::min_image_*`
causes the compiler to reload `cell_(0,0)` / `half_lx_` etc. through
`this` on every neighbor-pair iteration instead of keeping them in
registers as the original hoisted-doubles form does. The pe_profile
test showed +6.8% single-thread regression with the substitution; the
revert recovered to +0.6%. The duplicate `static inline minimum_image`
in lj.cpp was dead code and is gone — the inline branchy block in
`LJ::launch` is the kernel itself, not a duplicate of a free helper.

Bonded kernels (HarmonicBond/Angle/Dihedral) are not the bottleneck
(<10% of step cost) and can use `f.box.min_image(...)` cleanly without
a measurable regression at this resolution.
