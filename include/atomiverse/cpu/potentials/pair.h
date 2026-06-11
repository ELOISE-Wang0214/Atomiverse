/// @file cpu/potentials/pair.h
/// @brief LAMMPS-style abstract base for non-bonded pair styles (CPU).
///
/// Concrete pair styles (PairLJCut, future PairMorse, PairBuckingham, ...)
/// derive from this class, override compute(), and provide their own
/// neighbor-list inner loop. The base class adds NO virtual call inside
/// the per-pair loop — `launch()` is `final` and dispatches once per step
/// to `compute()`. Same cost profile as today's direct evaluators.
///
/// Derives from `cpu::PairEvaluator` (one cycle) so the existing
/// `CPUDriver::set_pair` `static_assert` keeps working unchanged.
#pragma once

#include "atomiverse/cpu/component/pair_evaluator.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Abstract base for non-bonded pair styles.
///
/// Derived classes implement `compute(ctx)` containing the neighbor-list
/// walk + per-pair energy/force math. `pair_cutoff()` reports the
/// interaction cutoff for the framework's `NeighborList` integration.
class Pair : public PairEvaluator {
  public:
    ~Pair() override = default;

    /// @brief Structural role — pair styles are evaluators.
    ComponentRole role() const final { return ComponentRole::Evaluator; }

    /// @brief Frame fields a pair style may write. Override to narrow.
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }

    /// @brief Framework dispatch entry — `final`, calls `compute()`.
    void launch(const LaunchContext& ctx) final { compute(ctx); }

    /// @brief Per-type pair cutoff [Å]. Reports the value stored at
    /// construction; derived classes can override to compute it from
    /// per-type tables.
    Real pair_cutoff() const override { return cutoff_; }

  protected:
    explicit Pair(Real cutoff) : cutoff_(cutoff) {}

    /// @brief The hot-path override. Derived class implements the
    /// neighbor-list walk, SIMD ket loop, per-thread force buffers,
    /// virial accumulation. No virtual call inside the per-pair loop.
    virtual void compute(const LaunchContext& ctx) = 0;

    Real cutoff_ = 0.0;
};

} // namespace cpu
} // namespace ATV
