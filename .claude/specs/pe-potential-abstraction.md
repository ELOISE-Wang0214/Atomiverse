---
title: Pair / Bond / Angle / Dihedral base classes — LAMMPS-style potential hierarchy on CPU and CUDA
slug: pe-potential-abstraction
status: in-progress
created: 2026-05-19
last-updated: 2026-05-19
parent: pe-ff-migrate-dreiding-ua.md
---

# Pair / Bond / Angle / Dihedral base classes (LAMMPS-style)

## Summary

Regularise the Atomiverse potential layer with four LAMMPS-named base classes on each backend: `cpu::Pair`, `cpu::Bond`, `cpu::Angle`, `cpu::Dihedral` (deriving from `cpu::ComponentBase`), and `cuda::Pair` / `cuda::Bond` / `cuda::Angle` / `cuda::Dihedral` (deriving from `cuda::ComponentBase`). The CPU and CUDA hierarchies are **independent** — different base classes — but use the **same names** and provide the **same shape** so user code and reviewers see one mental model. Rewrite the existing potentials bottom-up, renaming to LAMMPS terminology: `LJ` → `PairLJCut`, `HarmonicBond` → `BondHarmonic`, `HarmonicAngle` → `AngleHarmonic`, `MultiHarmonicDihedral` → `DihedralMultiHarmonic`. Each derived class owns its own inner loop (LAMMPS pattern, same as today) — there is no template / functor / "walker" machinery. The one virtual call per `compute()` per step is unchanged from the current cost; zero new overhead. Adding a new potential form (e.g. `PairMorse`, `BondFENE`) becomes "write one new derived class in the matching directory, register it in the user API." CUDA mirrors the rewrite (same base class shape, same derived class names) but ships **compile-only**, no tests, per user direction.

## Domain basis

No physics change in this spec. The naming hierarchy and the rewrite are organisational — every potential currently shipping in Atomiverse continues to compute the same energies and forces it does today (Dreiding-UA after `pe-ff-migrate-dreiding-ua` lands, OPLS-like before; either way unchanged by this spec). The win is in **architectural regularity**: a future contributor adding a new pair style writes `class PairMorse : public cpu::Pair` in `include/atomiverse/cpu/potentials/pair_morse.h` and the system knows how to compose it because every existing pair style follows the same shape.

LAMMPS has lived with exactly this hierarchy since 2003 and remains performance-competitive on modern hardware; the pattern is sufficient. We adopt the names verbatim so the mental cost of switching between LAMMPS and Atomiverse for a developer is near-zero.

## Design

### Naming and layering

Backend-rooted base classes (per backend, separate hierarchies):

```
cpu::ComponentBase                          (existing — `cpu::Pair` derives from this)
  ├── cpu::Pair                             (new — non-bonded pair styles)
  │     └── cpu::PairLJCut                  (new — was `cpu::md::LJ`)
  │     └── cpu::PairMorse                  (new, smoke-test only — demonstrates extension)
  ├── cpu::Bond                             (new — bonded 2-body styles)
  │     └── cpu::BondHarmonic               (new — was `cpu::HarmonicBond`)
  ├── cpu::Angle                            (new — bonded 3-body styles)
  │     └── cpu::AngleHarmonic              (new — was `cpu::HarmonicAngle`)
  └── cpu::Dihedral                         (new — bonded 4-body styles)
        └── cpu::DihedralMultiHarmonic      (new — was `cpu::MultiHarmonicDihedral`)
```

Mirror on the GPU side:

```
cuda::ComponentBase                         (existing)
  ├── cuda::Pair
  │     └── cuda::PairLJCut                 (was `cuda::LJ`)
  ├── cuda::Bond
  │     └── cuda::BondHarmonic              (was `cuda::HarmonicBond`)
  ├── cuda::Angle
  │     └── cuda::AngleHarmonic             (was `cuda::HarmonicAngle`)
  └── cuda::Dihedral
        └── cuda::DihedralMultiHarmonic     (was `cuda::MultiHarmonicDihedral`)
```

The two backends share no common ancestor for potentials — they share **shape** (same four base names, same derived names, same role/capability semantics) but not types. This is the regularity the user asked for: a developer reading either side immediately recognises the other.

### Base-class shape

All four bases follow the same skeleton. Example for `cpu::Pair`:

```cpp
// include/atomiverse/cpu/potentials/pair.h
namespace ATV::cpu {

/// @brief LAMMPS-style abstract base for non-bonded pair styles.
///
/// Every concrete pair style (PairLJCut, PairMorse, ...) derives from
/// this class, overrides compute(), and provides its own neighbor-list
/// inner loop. The base class itself adds NO virtual call inside the
/// per-pair loop — the only virtual dispatch is the one call to
/// compute() per launch(), which the Driver invokes once per step.
class Pair : public ComponentBase {
public:
    ~Pair() override = default;

    /// Pair styles are Evaluators that write Energy, Force, and Stress.
    ComponentRole role() const final { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }

    /// The framework calls launch() once per step; derived classes
    /// override compute() to do the actual neighbor walk + energy/
    /// force math.
    void launch(const LaunchContext& ctx) final { compute(ctx); }

    /// Per-type pair cutoff (Angstrom). Used by the NeighborList
    /// component to size the skin and bins. Derived classes report
    /// their cutoff via the constructor; the base stores it.
    double cutoff() const { return cutoff_; }

protected:
    explicit Pair(double cutoff) : cutoff_(cutoff) {}

    /// The hot-path override. Derived classes implement this — the
    /// neighbor-list walk, the SIMD ket loop, the per-thread force
    /// accumulation are all here. No virtual call inside the inner
    /// loop because `compute` is non-virtual at the call site
    /// (final dispatch via launch).
    virtual void compute(const LaunchContext& ctx) = 0;

    double cutoff_ = 0.0;
};

} // namespace ATV::cpu
```

`cpu::Bond` / `cpu::Angle` / `cpu::Dihedral` follow the same shape, with a few differences:

- They do not declare `cutoff()` (cutoff is implicit in the topology block — atoms with no bond, no force).
- They default `capabilities()` to `Energy | Force` (no Stress for bonded by default; derived class may extend if it writes virial — `BondHarmonic` does, so it overrides).
- They expose a protected helper `select_owned_tuples(const LaunchContext&, std::vector<int>& out)` that runs the multi-rank ownership filter (the `local_bonds_` cache currently duplicated in each bonded `launch()`). Derived classes call this once per step before iterating their topology block.

The CUDA base classes have the same shape but their `compute()` launches a CUDA kernel rather than running an OpenMP loop.

### Bottom-up rewrite plan

The migration runs **bottom-up** — bases first, then port one potential end-to-end as the canary, then port the remaining three. Each port is a single commit that the test suite gates on:

1. **Bases.** Land `cpu::Pair`, `cpu::Bond`, `cpu::Angle`, `cpu::Dihedral` headers + the multi-rank ownership helper. No derived class yet — only the suite of existing tests builds against the new code (they don't use the bases yet).

2. **First port: `PairLJCut`.** Move the body of `src/cpu/md/lj.cpp` into `src/cpu/potentials/pair_lj_cut.cpp` with `class PairLJCut : public cpu::Pair`. Constructor signature `PairLJCut(double eps, double sigma, double cutoff)` matches today's `LJ(double, double, double)`. Add a one-line backwards-compat alias `using LJ = PairLJCut` in the **old header** `include/atomiverse/cpu/md/lj.h` so existing user code (`sim.set_pair(cpu::md::LJ{...})`) keeps working unchanged. Run the full PE + LJ test suite — must remain green.

3. **`BondHarmonic`.** Same shape — move `harmonic_bond.cpp` body into `pair_bond_harmonic.cpp` with `class BondHarmonic : public cpu::Bond`. Use the `select_owned_tuples()` helper instead of the inlined `local_bonds_` cache. Add `using HarmonicBond = BondHarmonic` alias. Run suite.

4. **`AngleHarmonic`** and **`DihedralMultiHarmonic`** — same pattern. Each is one commit with the move + alias + suite re-run.

5. **CUDA mirror.** Land `cuda::Pair` / `cuda::Bond` / `cuda::Angle` / `cuda::Dihedral` bases and the four derived classes (`cuda::PairLJCut`, etc.) on the CUDA side. **Compile-only check** — `cmake --build build --target atomiverse_cuda` must succeed. No CUDA tests added (per user direction).

6. **Live demo: `PairMorse`.** Write `cpu::PairMorse : public cpu::Pair` from scratch as the proof that "add a new potential = one new class". The implementation reuses the same `neighbor-list walk → r² → energy/force` skeleton as `PairLJCut`, swapping the LJ math for Morse (`V = D(1 − e^{−α(r−r₀)})²`). One smoke test `test_pair_morse_smoke.cpp` runs a 2-atom Morse-bonded NVE for 100 steps, asserts E conservation at 1e-4 and equilibrium r* ≈ r₀.

7. **Cleanup.** After all four potentials are ported and tests are green, delete the empty old `.cpp` files (their bodies now live in `pair_*.cpp` / `bond_*.cpp` / `angle_*.cpp` / `dihedral_*.cpp`). Leave the old header aliases in place for backwards compatibility — they cost nothing.

### Zero-overhead argument

There is **no template machinery and no new virtual call** introduced by this spec. The cost profile is identical to today:

- `Driver::run()` calls `Pair::launch(ctx)` once per step (one virtual call, already there).
- `Pair::launch` is `final` and inlines a single direct call to `compute(ctx)` (also virtual, but pinned to the derived class — modern compilers de-virtualise this through devirtualisation of the `launch → compute` two-step when the static type is known at instantiation; even when they don't, that's still one indirect call per step, not per pair).
- Inside `compute()`, the derived class runs its own neighbor walk / SIMD ket loop, exactly as today. No new code in the inner loop.

The "same as LAMMPS" pattern is sufficient for performance; LAMMPS itself runs at near-peak with this hierarchy on production polymer simulations.

The performance gate (ac-005, ac-006 below) measures LJ::launch and HarmonicBond::launch ms-per-1000-steps against the pre-refactor `profile_1.csv` numbers. Target ±2%. Any larger regression indicates the refactor introduced a real cost (most likely a missed `final` annotation or a copy that wasn't elided), and the spec halts until it's resolved.

## Files to create or modify

### New (CPU)

```
include/atomiverse/cpu/potentials/pair.h
include/atomiverse/cpu/potentials/bond.h
include/atomiverse/cpu/potentials/angle.h
include/atomiverse/cpu/potentials/dihedral.h
include/atomiverse/cpu/potentials/pair_lj_cut.h
include/atomiverse/cpu/potentials/bond_harmonic.h
include/atomiverse/cpu/potentials/angle_harmonic.h
include/atomiverse/cpu/potentials/dihedral_multi_harmonic.h
include/atomiverse/cpu/potentials/pair_morse.h            (live extension demo)
src/cpu/potentials/pair.cpp                                (base — ownership helper, defaults)
src/cpu/potentials/bond.cpp                                (base — select_owned_tuples helper)
src/cpu/potentials/angle.cpp
src/cpu/potentials/dihedral.cpp
src/cpu/potentials/pair_lj_cut.cpp                         (body moved from src/cpu/md/lj.cpp)
src/cpu/potentials/bond_harmonic.cpp                       (body moved from harmonic_bond.cpp)
src/cpu/potentials/angle_harmonic.cpp
src/cpu/potentials/dihedral_multi_harmonic.cpp
src/cpu/potentials/pair_morse.cpp                          (live extension demo)
tests/cpu/md/test_pair_lj_cut.cpp                          (regression vs pre-refactor LJ)
tests/cpu/md/test_bond_harmonic.cpp                        (regression vs pre-refactor HarmonicBond)
tests/cpu/md/test_angle_harmonic.cpp
tests/cpu/md/test_dihedral_multi_harmonic.cpp
tests/cpu/md/test_pair_morse_smoke.cpp                     (demonstrates "one new class" workflow)
```

### Modified (CPU — backwards-compat aliases)

```
include/atomiverse/cpu/md/lj.h                             → `using LJ = PairLJCut;` + include redirect
include/atomiverse/cpu/potentials/harmonic_bond.h          → `using HarmonicBond = BondHarmonic;`
include/atomiverse/cpu/potentials/harmonic_angle.h         → `using HarmonicAngle = AngleHarmonic;`
include/atomiverse/cpu/potentials/multi_harmonic_dihedral.h → `using MultiHarmonicDihedral = DihedralMultiHarmonic;`
src/cpu/md/lj.cpp                                          → delete (body moved)
src/cpu/potentials/harmonic_bond.cpp                       → delete
src/cpu/potentials/harmonic_angle.cpp                      → delete
src/cpu/potentials/multi_harmonic_dihedral.cpp             → delete
```

### New (CUDA mirror, compile-only)

```
include/atomiverse/cuda/potentials/pair.cuh
include/atomiverse/cuda/potentials/bond.cuh
include/atomiverse/cuda/potentials/angle.cuh
include/atomiverse/cuda/potentials/dihedral.cuh
include/atomiverse/cuda/potentials/pair_lj_cut.cuh
include/atomiverse/cuda/potentials/bond_harmonic.cuh
include/atomiverse/cuda/potentials/angle_harmonic.cuh
include/atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh
src/cuda/potentials/pair_lj_cut.cu                         (body moved from src/cuda/potentials/lj.cu)
src/cuda/potentials/bond_harmonic.cu
src/cuda/potentials/angle_harmonic.cu
src/cuda/potentials/dihedral_multi_harmonic.cu
```

### Modified (CUDA — backwards-compat aliases)

```
include/atomiverse/cuda/potentials/lj.cuh                   → `using LJ = PairLJCut;`
include/atomiverse/cuda/potentials/harmonic_bond.cuh        → `using HarmonicBond = BondHarmonic;`
include/atomiverse/cuda/potentials/harmonic_angle.cuh       → `using HarmonicAngle = AngleHarmonic;`
include/atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh → `using MultiHarmonicDihedral = DihedralMultiHarmonic;`
src/cuda/potentials/{lj,harmonic_bond,harmonic_angle,multi_harmonic_dihedral}.cu → delete
```

### Notes

- `.claude/notes/potential-abstraction-audit.md` — consumed (produced by `pe-ff-migrate-dreiding-ua`'s audit task).
- `.claude/notes/potential-base-classes-design.md` — written by this spec's first task; one-pass design pass before code lands.

## Tasks

- [ ] **Design pass.** Read `.claude/notes/potential-abstraction-audit.md`. Author `.claude/notes/potential-base-classes-design.md` covering: (a) the four base-class headers as drafted in this spec, (b) the bottom-up port sequence (PairLJCut canary first), (c) the multi-rank ownership helper signature on `Bond` / `Angle` / `Dihedral` and how it replaces the current `local_bonds_` cache in each derived class, (d) the `final` annotation strategy on `launch()` for de-virtualisation, (e) backwards-compat alias policy in the old headers so user code does not break, (f) CUDA mirror layout. Get explicit sign-off in the note before any source file lands.
- [ ] **Land bases.** Create the four `cpu::Pair` / `cpu::Bond` / `cpu::Angle` / `cpu::Dihedral` headers + their `.cpp` files containing only the protected helpers (`select_owned_tuples` etc.) and trivial constructors. Build green; no existing tests should break (they don't use the bases yet).
- [ ] **Port `PairLJCut` (canary).** Write `tests/cpu/md/test_pair_lj_cut.cpp` first — a RED test asserting `cpu::PairLJCut` exists and matches the pre-refactor `cpu::md::LJ` output bit-for-bit (or at 5e-12 FMA-tolerance). Move the body of `src/cpu/md/lj.cpp` into `src/cpu/potentials/pair_lj_cut.cpp` with the new class. Add the `using LJ = PairLJCut` alias in `include/atomiverse/cpu/md/lj.h`. Delete `src/cpu/md/lj.cpp`. Verify the new test is GREEN, all existing `test_lj_*` tests are still green, and `pe_profile_single` reports LJ::launch within ±2% of the 154.79 ms baseline.
- [ ] **Port `BondHarmonic`.** RED test → move body into `src/cpu/potentials/bond_harmonic.cpp` → alias the old name → delete the old `.cpp` → verify `test_pe_*` tests still pass and `pe_profile_single` HarmonicBond::launch within ±2% of 117.36 ms.
- [ ] **Port `AngleHarmonic`.** Same pattern.
- [ ] **Port `DihedralMultiHarmonic`.** Same pattern.
- [ ] **CUDA mirror (compile-only).** Add the four CUDA bases + four CUDA derived classes mirroring the CPU layout. Add backwards-compat aliases on the old CUDA headers. Delete old `.cu` files after moves. `cmake --build build --target atomiverse_cuda` must succeed with no new warnings. **No CUDA tests added.**
- [ ] **Live demo: `PairMorse`.** Write `cpu::PairMorse : public cpu::Pair` from scratch — the smallest reasonable new pair-style. Implementation: neighbor-list walk identical in shape to `PairLJCut`, swap the LJ math for `V = D(1 − e^{−α(r−r₀)})²` and its derivative. Smoke test `test_pair_morse_smoke.cpp` runs a 2-atom Morse-bonded NVE for 100 steps; asserts energy conservation at 1e-4 relative and equilibrium r* matches `r₀` after a small initial perturbation. Document in the design note the total **count of files touched** to add this new style — target: 2 (one `.h`, one `.cpp`, plus the test which doesn't count as part of the "adding a style" workflow itself).
- [ ] **Cleanup.** Run full `ctest --test-dir build --output-on-failure`; verify no regression vs the post-`pe-ff-migrate-dreiding-ua` baseline. Delete any straggler old `.cpp` files. Update `.claude/notes/potential-base-classes-design.md` with the final touched-file count for the Morse demo.

## Testing strategy

**Zero-overhead validation:**

- `pe_profile_single` LJ::launch total_ms within ±2% of `profile_1.csv` pre-refactor baseline (154.79 ms).
- `pe_profile_single` HarmonicBond::launch total_ms within ±2% of pre-refactor baseline (117.36 ms).
- `test_lj_nve` (parent ac-012, 5e-12 FMA tolerance) still passes after the LJ port.

**Happy path:**

- All existing `test_pe_*`, `test_lj_*`, `test_*multi_chain`, `test_mpi_pe_*`, `test_mpi_lj_*` continue green.
- Full ctest count after this spec lands: pre-baseline + 4 regression tests (one per port) + 1 Morse smoke = +5 tests.

**Edge cases:**

- A user who today writes `cpu::md::LJ{eps, sigma, cutoff}` must still compile and link without any change. The alias `using LJ = PairLJCut` carries this.
- Multi-rank bonded tests: `BondHarmonic` must continue to use the multi-rank ownership filter. The `select_owned_tuples` helper in the `Bond` base replaces the inline `local_bonds_` cache; verify by re-running `test_mpi_pe_chains_8rank_n8` and `test_mpi_pe_chain_drifting_n2`.

**Domain validation** (`$META.science.required: true`):

- `PairMorse` reaches its analytical equilibrium r* = r₀ in NVE within a few oscillation periods — verifies the new derived class actually implements Morse correctly (not just compiles).

**Performance:**

- LJ::launch and HarmonicBond::launch profile_1.csv numbers within ±2%.
- xsimd ket-loop assembly inspection (one-off, documented in the design note): verify the inner pair-distance / energy / force math still appears inside the SIMD unroll, not behind a `call` instruction.

## Out of scope

- **Runtime style selection via string keys.** LAMMPS exposes `pair_style lj/cut` at the input-script level. Atomiverse stays compile-time selected; a string-keyed registry can land later if Python-bindings users need it. Pure follow-up.
- **torch backend.** Out of scope per user direction.
- **CUDA functional tests.** Compile-only mirror per user direction.
- **Improper, KSpace, Fix base classes.** LAMMPS has more base classes than these four. We only port what Atomiverse already has analogues for — Improper (out-of-plane bonded) and KSpace (long-range) don't exist in tree, so no base for them in this spec.
- **Refactoring `cpu::md::NeighborList`** to belong under the Pair hierarchy. Today `NeighborList` is a standalone component and pair styles consume its output. That stays the same; making `Pair` own a neighbor-list reference is a separate ergonomic refactor.
- **Sharing a base class across CPU and CUDA.** Explicitly rejected by the user — the two backends have different `ComponentBase` ancestors and there is no value in forcing a common ancestor on the potentials.
- **Performance optimisation beyond preserving today's numbers.** The xsimd ket loop, OMP per-thread buffers, etc. all carry over verbatim. Further perf work is downstream.
- **`pe-ff-migrate-dreiding-ua` FF constant migration.** Orthogonal — that spec changes values, this one changes structure. They can run in either order.
