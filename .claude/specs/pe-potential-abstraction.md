---
title: Zero-overhead Potential abstraction shared by CPU and CUDA
slug: pe-potential-abstraction
status: approved
created: 2026-05-19
parent: pe-ff-migrate-dreiding-ua.md
---

# Zero-overhead Potential abstraction (CPU + CUDA)

## Summary

Introduce a two-layer `Potential` abstraction shared by `cpu::` and `cuda::` backends so adding a new force-field form (e.g. Morse bond, FENE bond, Buckingham pair) becomes a **single new header** with one functor instead of the current 8-file copy-edit ritual. Layer 1 is a header-only inline **Form** functor capturing the per-pair / per-tuple energy + derivative math — `static constexpr inline std::pair<double,double> eval(...)`. Layer 2 is a templated **Walker** component (`PairWalker<Form>`, `BondWalker<Form>`, `AngleWalker<Form>`, `DihedralWalker<Form>`) that owns the iteration infrastructure: neighbor-list walk for pairs, topology iteration for bonded terms, OpenMP per-thread force buffers, xsimd ket-loop for pair styles, multi-rank halo / ownership filtering, virial accumulation. Existing `LJ`, `HarmonicBond`, `HarmonicAngle`, `MultiHarmonicDihedral` collapse to single-line type aliases (`using LJ = PairWalker<LJForm>`). CUDA mirrors the same shape (one Form + one templated kernel) but does **not** ship tests (per user direction). Zero virtual-call overhead in the inner loop is the hard constraint — `Form::eval` is `static inline` so the compiler can fully inline it into the SIMD batch loop and into the CUDA per-thread kernel.

This spec consumes the audit findings produced by `pe-ff-migrate-dreiding-ua`'s architecture-audit task (delivered at `.claude/notes/potential-abstraction-audit.md`). It does **not** depend on the FF-constants migration landing first; the two specs are independent and can run in either order.

## Domain basis

No new physics. The abstraction wraps the same Dreiding-UA force field already implemented (harmonic bond/angle, multi-harmonic dihedral, LJ 12-6 with cutoff) plus opens the door for canonical alternates (Morse, FENE, Buckingham, Stillinger-Weber) without re-implementing the iteration machinery. The Form/Walker split mirrors LAMMPS-Kokkos's pair/bond style hierarchy, which has been validated for performance portability across CPU+GPU on production polymer simulations.

The zero-overhead constraint is critical because the inner loops are called O(N · ⟨n_neigh⟩) ≈ 10⁵–10⁷ times per step on a 10k–10⁶-atom system; any virtual call inside the ket loop would add tens of nanoseconds per pair × millions of pairs × thousands of steps and dominate the runtime. CRTP and template policy patterns achieve compile-time dispatch with no indirection.

## Design

### Layer 1 — Form (header-only inline functor)

A `Form` is a `struct` exposing one or two `static constexpr inline` member functions and a parameter type. Two flavours:

**Pair forms** (called per atom-pair inside the neighbor-list walk):

```cpp
// include/atomiverse/potentials/forms/lj.h  (backend-neutral, header-only)
namespace ATV::potentials {
struct LJForm {
    struct Params { double eps, sigma, rcut_sq, sig6, eps4, eps24, vshift; };
    // returns {E, dE_over_r2}  so the caller computes Fx_ij = dE_over_r2 * dx, etc.
    // Branchless cutoff: caller passes a mask; this returns 0 contributions when outside.
    static constexpr inline std::pair<double,double> eval(double r2, const Params& p) noexcept {
        double inv_r2 = 1.0 / r2;
        double s2 = p.sigma * p.sigma * inv_r2;
        double s6 = s2 * s2 * s2;
        double E  = p.eps4 * (s6 * s6 - s6) - p.vshift;
        double dE_over_r2 = p.eps24 * (2.0 * s6 * s6 - s6) * inv_r2;
        return {E, dE_over_r2};
    }
};
} // namespace ATV::potentials
```

**Bonded forms** (called per topology tuple — bond/angle/dihedral):

```cpp
struct HarmonicBondForm {
    struct Params { double K, r0; };
    // V(r) = K · (r − r0)²,  dV/dr = 2K(r − r0)
    static constexpr inline std::pair<double,double> eval(double r, const Params& p) noexcept {
        double dr = r - p.r0;
        return {p.K * dr * dr, 2.0 * p.K * dr};
    }
};
struct HarmonicAngleForm {
    struct Params { double K, theta0; };
    static constexpr inline std::pair<double,double> eval(double theta, const Params& p) noexcept {
        double dt = theta - p.theta0;
        return {p.K * dt * dt, 2.0 * p.K * dt};
    }
};
struct MultiHarmonicDihedralForm {
    struct Params { double C[5]; };
    static constexpr inline std::pair<double,double> eval(double phi, const Params& p) noexcept {
        // Σ_n C[n] · cos(φ)^n, with explicit dE/dφ in terms of cos(φ), sin(φ).
        // [implementation detail; see existing src/cpu/potentials/multi_harmonic_dihedral.cpp]
    }
};
```

`Form::eval` must:
- be `static`, `constexpr`, `noexcept`, and `inline`-friendly (no virtual, no heap, no exceptions).
- take a `const Params&` (not separate scalars) so users can extend the parameter set without changing the signature.
- be **header-only** so the compiler can inline it into both the CPU xsimd ket-loop and the CUDA `__device__` kernel.
- accept and return only `double` (or vectorized form — see SIMD section below); no allocation, no callbacks.

A `static_assert` in each Walker template enforces that `Form::eval` is invocable with `(double, const Params&)` and returns `std::pair<double,double>` (or the SIMD-batch equivalent).

### Layer 2 — Walker (templated component)

Four Walker templates, one per topology shape:

```cpp
// include/atomiverse/cpu/potentials/pair_walker.h
namespace ATV::cpu {
template <class Form>
class PairWalker : public ComponentBase {
public:
    explicit PairWalker(std::vector<typename Form::Params> params_by_type);

    std::string name() const override { return Form::name; }   // each Form exposes a static name
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }
    void build(const BuildContext& ctx) override;     // alloc per-thread force buffers, NL handles
    void launch(const LaunchContext& ctx) override;   // OMP-parallel neighbor walk; xsimd ket loop calls Form::eval
    void destroy(const DestroyContext& ctx) override;

private:
    std::vector<typename Form::Params> params_;
    // ...per-thread fxt/fyt/fzt buffers, neighbor-list handles, etc.
};
} // namespace ATV::cpu
```

Inner-loop hot path inside `PairWalker<Form>::launch()`:

```cpp
for (int jj = 0; jj + W <= jnum; jj += W) {
    // gather W neighbors -> xsimd::batch<double>
    // compute r²_batch, mask = (r²_batch < rcut²)
    // call Form::eval lane-by-lane OR Form::eval_batch if Form exposes one
    for (std::size_t lane = 0; lane < W; ++lane) {
        if (mask[lane]) {
            auto [E, dE_over_r2] = Form::eval(r2_lane[lane], params_[type_ij]);
            // accumulate forces, virial, energy — same scalar pattern as today
        }
    }
}
```

The `Form::eval` call is `static constexpr inline`, so the compiler inlines it into the lane-by-lane unroll — **zero virtual dispatch, zero indirection**. The xsimd ket-loop survives unchanged (the existing 1.21× speedup is preserved).

If a Form additionally exposes `static inline void eval_batch(xsimd::batch<double>&, const xsimd::batch<...>&, ...) noexcept`, the Walker dispatches the full SIMD batch through it (no lane-by-lane fan-out). This is the **opt-in fast path**; absent it, the lane-by-lane scalar fallback above is used. The `static_assert` checks for both shapes via a SFINAE / concept.

Bonded Walkers (`BondWalker<Form>`, `AngleWalker<Form>`, `DihedralWalker<Form>`) follow the same pattern: topology iteration is in the Walker, energy/force math is in the Form. They reuse the multi-rank ownership filter (`tag_to_local`) the current `HarmonicBond::launch()` already uses, the OMP per-thread force buffers, and the virial accumulation.

### Existing-class collapse

After the Walker templates land, the existing CPU classes become type aliases:

```cpp
// include/atomiverse/cpu/potentials/lj.h
namespace ATV::cpu {
using LJ = PairWalker<potentials::LJForm>;
} // namespace ATV::cpu
```

```cpp
// include/atomiverse/cpu/potentials/harmonic_bond.h
namespace ATV::cpu {
using HarmonicBond = BondWalker<potentials::HarmonicBondForm>;
} // namespace ATV::cpu
```

Same for `HarmonicAngle` and `MultiHarmonicDihedral`. The user-facing API (`sim.add(cpu::md::LJ{eps, sigma, cutoff})`) is preserved by giving each Walker a constructor that forwards to `Form::Params`. The `.cpp` files in `src/cpu/{md,potentials}/` collapse to a few hundred lines of Walker boilerplate (one .cpp per Walker shape, not per Form) plus the inline Form headers under `include/atomiverse/potentials/forms/`.

### CUDA mirror (compile-only per user direction)

A mirror set of Walker templates lives in `include/atomiverse/cuda/potentials/`:

```cpp
template <class Form>
class PairWalker : public ComponentBase {
    // build() allocates device buffers
    // launch() launches a __global__ kernel that calls Form::eval inside the per-thread loop
};
```

The `Form` structs are reused as-is — they're `static constexpr inline`, callable from both host and device. `Form::eval` gets a `__host__ __device__` annotation either inline at the Form declaration site, or via a macro `ATV_HOST_DEVICE` that expands to `__host__ __device__` under CUDA and to nothing on the CPU-only path.

CUDA Walkers are not unit-tested in this spec — per user direction, the CUDA backend ships a compile-only mirror. The smoke gate is: `cmake --build build --target atomiverse_cuda` succeeds. Future work (`cuda-driver-taxonomy` follow-up) will land CUDA tests.

### What changes for users

The user-visible API is **unchanged**:

```cpp
// Before and after — identical:
cpu::md::MD sim;
sim.set_pair(cpu::md::LJ{eps, sigma, cutoff});
sim.add(cpu::potentials::HarmonicBond{{K_eV}, {r0}});
```

Adding a new force-field form is **one new file**:

```cpp
// include/atomiverse/potentials/forms/morse.h
namespace ATV::potentials {
struct MorseForm {
    struct Params { double D, alpha, r0; };
    static constexpr inline std::pair<double,double> eval(double r, const Params& p) noexcept {
        double e = std::exp(-p.alpha * (r - p.r0));
        double v = p.D * (1.0 - e) * (1.0 - e);
        double dv = 2.0 * p.D * p.alpha * (1.0 - e) * e;
        return {v, dv};
    }
};
} // namespace ATV::potentials

// And one type alias in include/atomiverse/cpu/potentials/morse_bond.h:
namespace ATV::cpu {
using MorseBond = BondWalker<potentials::MorseForm>;
} // namespace ATV::cpu
```

Done. No new `.cpp`, no test-file copy-paste, no CUDA kernel to write by hand.

## Files to create or modify

### New (CPU + backend-neutral)

```
include/atomiverse/potentials/forms/lj.h
include/atomiverse/potentials/forms/harmonic_bond.h
include/atomiverse/potentials/forms/harmonic_angle.h
include/atomiverse/potentials/forms/multi_harmonic_dihedral.h
include/atomiverse/cpu/potentials/pair_walker.h
include/atomiverse/cpu/potentials/bond_walker.h
include/atomiverse/cpu/potentials/angle_walker.h
include/atomiverse/cpu/potentials/dihedral_walker.h
src/cpu/potentials/pair_walker.cpp
src/cpu/potentials/bond_walker.cpp
src/cpu/potentials/angle_walker.cpp
src/cpu/potentials/dihedral_walker.cpp
tests/cpu/md/test_potential_walker_lj.cpp           (LJ regression vs pre-refactor)
tests/cpu/md/test_potential_walker_harmonic.cpp     (bond + angle regression)
tests/cpu/md/test_potential_walker_dihedral.cpp     (dihedral regression)
tests/cpu/md/test_potential_walker_morse_smoke.cpp  (demonstrates the "one new file" workflow)
```

### Modified (CPU)

```
include/atomiverse/cpu/md/lj.h                          → type alias to PairWalker<LJForm>
include/atomiverse/cpu/potentials/harmonic_bond.h       → type alias to BondWalker<HarmonicBondForm>
include/atomiverse/cpu/potentials/harmonic_angle.h      → type alias to AngleWalker<HarmonicAngleForm>
include/atomiverse/cpu/potentials/multi_harmonic_dihedral.h → type alias
src/cpu/md/lj.cpp                                       → collapse / delete (logic moves to pair_walker.cpp)
src/cpu/potentials/harmonic_bond.cpp                    → collapse / delete
src/cpu/potentials/harmonic_angle.cpp                   → collapse / delete
src/cpu/potentials/multi_harmonic_dihedral.cpp          → collapse / delete
```

### New (CUDA mirror, compile-only)

```
include/atomiverse/cuda/potentials/pair_walker.cuh
include/atomiverse/cuda/potentials/bond_walker.cuh
include/atomiverse/cuda/potentials/angle_walker.cuh
include/atomiverse/cuda/potentials/dihedral_walker.cuh
src/cuda/potentials/pair_walker.cu     (templated CUDA Walker .cu, ETI per Form)
src/cuda/potentials/bond_walker.cu
src/cuda/potentials/angle_walker.cu
src/cuda/potentials/dihedral_walker.cu
```

### Modified (CUDA mirror)

```
include/atomiverse/cuda/potentials/{harmonic_bond,harmonic_angle,multi_harmonic_dihedral,lj}.cuh
                                                        → type aliases on the same Walker<Form>
src/cuda/potentials/{harmonic_bond,harmonic_angle,multi_harmonic_dihedral,lj}.cu  → collapse
```

### Macros / common utilities

```
include/atomiverse/potentials/host_device.h   → ATV_HOST_DEVICE macro (expands to __host__ __device__ under NVCC, empty otherwise)
```

### Other

- `.claude/notes/potential-abstraction-audit.md` — written by `pe-ff-migrate-dreiding-ua`'s audit task; this spec consumes it.
- `.claude/notes/potential-abstraction-design.md` — written by this spec's first task; captures the design decisions (CRTP vs template-policy, what `Form::eval` shape, SIMD opt-in, CUDA host/device annotation strategy). One pass; signed off before code lands.

## Tasks

- [ ] **Design pass.** Read `.claude/notes/potential-abstraction-audit.md` (produced by `pe-ff-migrate-dreiding-ua`). Author `.claude/notes/potential-abstraction-design.md` covering: (a) Form contract (`static constexpr inline`, Params struct, eval signature, optional `eval_batch`); (b) Walker class skeletons (CPU and CUDA); (c) `ATV_HOST_DEVICE` macro strategy; (d) SIMD opt-in fast path via SFINAE or C++20 concepts; (e) how multi-rank halo / ghost filtering survives the refactor (it lives in the Walker, the Form is unaware); (f) name-resolution: where the Form headers live (`include/atomiverse/potentials/forms/`) and how user code finds them.
- [ ] Write failing test `tests/cpu/md/test_potential_walker_lj.cpp` — exercises the not-yet-existing `cpu::PairWalker<LJForm>` on a known geometry, compares against the pre-refactor LJ output bit-for-bit (or at the documented 5e-12 FMA-reorder tolerance per parent ac-012 refinement). Confirm RED.
- [ ] Implement `include/atomiverse/potentials/forms/lj.h` (the LJForm functor with Params, eval, optional eval_batch).
- [ ] Implement `include/atomiverse/cpu/potentials/pair_walker.h` + `src/cpu/potentials/pair_walker.cpp` (neighbor-list walk, OpenMP per-thread buffers, xsimd ket-loop with `Form::eval` lane fan-out, multi-rank halo). The existing `src/cpu/md/lj.cpp` body becomes the seed implementation, generalised to call `Form::eval` instead of inlined `s12 - s6` math.
- [ ] Collapse `cpu::md::LJ` to `using LJ = PairWalker<LJForm>` in `include/atomiverse/cpu/md/lj.h`; delete `src/cpu/md/lj.cpp`. Confirm `test_potential_walker_lj.cpp` GREEN; confirm existing `test_lj_*` tests still pass (LJ NVE bit-identical, multi-rank LJ correctness, etc.).
- [ ] Write failing test `tests/cpu/md/test_potential_walker_harmonic.cpp` (bond + angle); implement `HarmonicBondForm` / `HarmonicAngleForm`; implement `BondWalker<>` + `AngleWalker<>` (templated component, topology iteration + multi-rank ownership filter); collapse the existing `cpu::HarmonicBond` / `cpu::HarmonicAngle`; verify all `test_pe_*` tests still pass.
- [ ] Write failing test for dihedral; implement `MultiHarmonicDihedralForm`; implement `DihedralWalker<>`; collapse `cpu::MultiHarmonicDihedral`; verify all PE tests still pass.
- [ ] **Demonstrate "one new file" for adding Morse bond.** Write `include/atomiverse/potentials/forms/morse.h` + one-line `using MorseBond = BondWalker<MorseForm>` in a test-only header `tests/cpu/md/test_potential_walker_morse_smoke.h`. Test `test_potential_walker_morse_smoke.cpp` builds a 2-atom Morse-bonded system, runs 100 NVE steps, asserts (a) energy conservation at 1e-4 relative, (b) the equilibrium r* matches `r0`. **Count and document the total touched files for this addition** in the design note — should be exactly 2 (one form header + one test). This is the live proof of zero-paint-cost extension.
- [ ] CUDA mirror: implement `cuda::PairWalker<Form>` + `cuda::BondWalker<Form>` + `cuda::AngleWalker<Form>` + `cuda::DihedralWalker<Form>` against `__device__`-callable Form::eval (annotate with `ATV_HOST_DEVICE`). Collapse the existing CUDA potential classes to type aliases. **No CUDA tests required** (per user direction); the smoke gate is `cmake --build build --target atomiverse_cuda` compiles clean.
- [ ] Performance gate: re-run `pe_profile_single --steps 1000 --warmup 100` on the new code; compare `LJ::launch` total_ms to the pre-refactor profile_1.csv value (154.79 ms). The new value must be within ±2% (no perf regression from the templated abstraction). If it regresses, find the inlining failure and fix before merging.
- [ ] Run full `ctest --test-dir build --output-on-failure`; verify no regression vs the 88-test baseline after `pe-ff-migrate-dreiding-ua` (or 87 if that spec hasn't landed yet). New tests above bring the count to 91+.

## Testing strategy

**Zero-overhead validation** (the critical correctness gate):

- `pe_profile_single` LJ::launch total_ms before vs after refactor: within ±2% (if larger, the Form::eval inlining failed). This is the load-bearing performance gate — a templated abstraction whose runtime cost dominates the savings is a failure.
- Per-pair bit-identicality: `test_lj_nve` (the parent ac-012 5e-12 FMA-tolerance test) continues to pass — the refactor moves code around but should not change FP behaviour.
- Bonded per-term energies: `test_pe_single_chain` continues to pass at 1e-9 absolute (this verifies BondWalker / AngleWalker / DihedralWalker compute the same energies as the old classes).

**Happy path:**

- Single-rank: every existing `test_pe_*`, `test_lj_*`, `test_*_multi_chain` continues green.
- Multi-rank: every existing `test_mpi_*pe*`, `test_mpi_*lj*` continues green at 1, 2, 8 ranks.
- Build: `cmake --build build` clean on Release + Debug.

**Edge cases:**

- Form with **no** `eval_batch` falls back to lane-by-lane scalar fan-out from the xsimd batch — verify this path still SIMDs (the lane loop unrolls to W copies of the scalar formula, which clang vectorises).
- Form with **both** `eval` and `eval_batch`: SFINAE picks `eval_batch` automatically; assert via `static_assert` and a small "dispatch trace" test that flags which path was taken.
- Multi-rank ghost ownership: BondWalker / AngleWalker / DihedralWalker must reuse the existing `dd->tag_to_local` filter; verify by re-running `test_mpi_pe_chains_8rank_n8` and `test_mpi_pe_chain_drifting_n2`.

**Domain validation** (`$META.science.required: true`):

- The Morse-bond smoke test (`test_potential_walker_morse_smoke`) demonstrates a never-before-implemented FF form working end-to-end in <2 hours of editing. This is the live demonstration that "add a new FF = one file" — the very claim the abstraction makes.

**Performance:**

- Pre/post LJ::launch ms/1000-step: ±2% target.
- Pre/post HarmonicBond::launch ms/1000-step: ±2% target.
- xsimd ket-loop assembly inspection (one-off, written into the design note): verify the inlined `Form::eval` body appears inside the SIMD unroll, not behind a call instruction.

## Out of scope

- **Runtime FF selection via string keys** (LAMMPS-style `pair_style lj/cut`). The abstraction makes future addition of a `Registry<std::string, PotentialFactory>` cheap, but this spec only enables compile-time form dispatch. A string-keyed registry is the natural follow-up if Python-bindings users need it.
- **torch backend.** Per user direction, this spec covers CPU and CUDA only.
- **CUDA test coverage.** Compile-only smoke; CUDA functional tests land separately under the `cuda-driver-taxonomy` umbrella.
- **Hybrid CPU/CUDA Forms.** Each Form is one struct usable from either backend's Walker; this spec doesn't introduce a "this Form runs on the GPU half and the CPU half" mixed-mode kernel.
- **Long-range / k-space potentials (Ewald, PPPM).** The Walker contract assumes short-range cutoff; long-range methods need a different iteration shape and are out of scope.
- **SCF / QC potentials.** The Hartree-Fock / DFT machinery in `cpu::scf` is owned by `cpu::scf::SCFEvaluator` (an Evaluator component, but not a "potential" in the pair/topology sense). No change here.
- **`pe-ff-migrate-dreiding-ua` FF constant migration.** This spec doesn't touch the OPLS-like vs Dreiding-UA value debate; whichever values live in the existing classes will continue to live in `Form::Params` defaults after the refactor. The two specs are orthogonal.
