---
title: Pair / Bond / Angle / Dihedral base-class design (LAMMPS-style)
slug: potential-base-classes-design
status: design
created: 2026-05-19
owed-to: .claude/specs/pe-potential-abstraction.md
---

# Pair / Bond / Angle / Dihedral bases — design

Resolves audit gaps 1, 2, 3, 6 (HIGH). Spec: independent
per-backend hierarchies, no shared CPU/CUDA ancestor, classic
LAMMPS virtual `compute()` once per step, no templates / functors.

## (a) CPU base headers

`include/atomiverse/cpu/potentials/pair.h`:

```cpp
#pragma once
#include "atomiverse/cpu/component.h"
namespace ATV::cpu {
class Pair : public ComponentBase {
public:
    ~Pair() override = default;
    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }
    bool needs_neighbor_list() const final { return true; }
    double cutoff() const { return cutoff_; }
    void launch(const LaunchContext& ctx) final { compute(ctx); }
protected:
    explicit Pair(double cutoff) : cutoff_(cutoff) {}
    virtual void compute(const LaunchContext& ctx) = 0;
    double cutoff_ = 0.0;
};
} // namespace ATV::cpu
```

`bond.h` is the same skeleton minus `cutoff()` /
`needs_neighbor_list()`, plus a protected helper that absorbs the
inline `local_bonds_` cache from `harmonic_bond.cpp:53-67`:

```cpp
/// Indices into frame.bonds owned by this rank; cached against
/// dd->neighbor_generation.
const std::vector<int>& select_owned_tuples(const LaunchContext& ctx);
```

Private state: `std::vector<int> owned_` + `std::uint64_t
cached_generation_`. `angle.h` / `dihedral.h` identical, reading
`frame.angles` / `frame.dihedrals`.

## (b) CUDA mirror

Same shape, derived from `cuda::ComponentBase`. CUDA bases ship
`role()` / `capabilities()` — **adds the two virtuals to
`cuda::ComponentBase`** with CPU defaults (resolves gap 3).
Compile-only. `Bond::select_owned_tuples()` is a
shim returning `0..n_bonds` until DD lands.

## (c) Port sequence + risks

1. **Bases (CPU).** `CPUDriver::set_pair`'s `static_assert` is
   keyed to `cpu::PairEvaluator` (`cpu_driver.h:191`). Make `Pair`
   derive from `PairEvaluator` one cycle.
2. **PairLJCut canary.** Move `src/cpu/md/lj.cpp` (504 lines) →
   `src/cpu/potentials/`. Special-bonds CSR cache +
   `pair_cutoff()` must survive verbatim. Gate: LJ::launch
   154.79 ms ±2 %.
3. **BondHarmonic.** Cache + `multi_rank ? local_bonds_[ib] : ib`
   dispatch (`harmonic_bond.cpp:96-98`) must move into
   `Bond::select_owned_tuples()` *without changing semantics*.
   Owner rule (`owner_tag = min(ti, tj)`, lines 59-62) is
   bond-specific; angle/dihedral use `min` over 3/4 tags. Helper
   takes a lambda `owner_tag_of(idx)` — base override count zero.
4. **AngleHarmonic / DihedralMultiHarmonic.** Same pattern. Gate:
   `test_pe_chains_8rank_n8`, `test_mpi_pe_chain_drifting_n2`
   green.
5. **CUDA mirror.** Compile-only. Inline the kernel in `.cu`
   (skip `_kernel.cuh`).
6. **PairMorse demo.** Risk: Morse math correctness only.

## (d) `final` strategy / de-virtualisation

`launch() final → compute() virtual`. `launch` is `final`;
`comp->launch(ctx)` (Driver indirect call) unchanged. The direct call to `compute(ctx)` is the only virtual call per step — identical cost to today. No per-pair virtuals. Compilers inline `launch → compute` only when static type known (rare); even when not, one indirect branch per component per step, not per pair. Spec ±2 % gates (154.79 ms, 117.36 ms) catch regressions.

## (e) Backwards-compat aliases

Old headers become one-line redirects:

```cpp
// include/atomiverse/cpu/md/lj.h
#pragma once
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
namespace ATV::cpu::md { using LJ = ATV::cpu::PairLJCut; }
```

Bonded shim same shape into bare `ATV::` (matching gap 2):
`namespace ATV { using HarmonicBond = ATV::cpu::BondHarmonic; }`.
User code (e.g. `sim.set_pair(cpu::md::LJ{...})`) compiles unchanged. Aliases zero-cost. Classes move to `ATV::cpu::` (partial
fix for gap 2).

## (f) Live extension demo — `PairMorse` file count

**Status: confirmed 2026-05-19** — prediction matched exactly.

Prediction: **2 files** plus one uncounted test:

  - `include/atomiverse/cpu/potentials/pair_morse.h`
  - `src/cpu/potentials/pair_morse.cpp`
  - `tests/cpu/md/test_pair_morse_smoke.cpp` (uncounted)

Zero CMake edits if `src/cpu/potentials/*.cpp` is a `GLOB` (today
explicit at `CMakeLists.txt:179-183` — convert in the bases-land
commit). If GLOB rejected, count becomes 3. Either is a 75 % cut
vs today's 8.

**Verified result**: GLOB landed in the bases commit; PairMorse added
with exactly 2 production files (`pair_morse.h` 45 lines +
`pair_morse.cpp` 55 lines), 0 CMakeLists.txt edits to production
code. Smoke test `test_pair_morse_smoke.cpp` (uncounted; verifies
the demo, not part of the "add a style" workflow) passes a 2-atom
NVE energy-conservation gate in 19 ms.

## Gap-coverage footer

  - Gap 1 (no shared base) — (a)/(b).
  - Gap 2 (CPU bare-`ATV::` bonded) — (a)/(e).
  - Gap 3 (CUDA missing `role`/`capabilities`) — (b).
  - Gap 6 (8-file Morse add) — (f): target 2.

No HIGH gap unaddressed. MEDIUM gap 4 (no CUDA `set_pair`)
deferred to `cuda-driver-taxonomy.md`. LOW gap 5 (string-keyed
registry) out of scope per spec.
