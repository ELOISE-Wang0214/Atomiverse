/// @file cpu/component/pair_evaluator.h
/// @brief PairEvaluator — backend-root base for cutoff-pair force fields.
///
/// Lives at `cpu/` root (not `cpu/md/`) because cutoff-pair forces are
/// shared infrastructure across multiple domains:
///   - MD (`cpu::md::MD`) — pair forces in the per-step force kernel.
///   - Minimizer (`cpu::minimizer::Minimize`) — same evaluators reused
///     in geometry optimisation.
///   - Evaluator (`cpu::Evaluate`) — single-point pair energy / force.
///   - Future Monte Carlo / WIDOM / replica-exchange drivers — same
///     pair-cutoff abstraction, different outer loop.
///
/// Concrete subclasses (today: `cpu::PairLJCut`) live in their domain
/// directory; only the abstraction is at root.
///
/// The `NeighborConsumer` protocol is implemented here so subclasses
/// don't have to remember it — every cutoff-pair force field needs
/// a `NeighborList` upstream. The CPU category drivers throw at
/// `build()` time when a consumer is registered and no `set_neighbor()`
/// was made.
#pragma once

#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/neighbor_protocol.h"

namespace ATV {
namespace cpu {

/// @brief Backend-root base for cutoff-pair force fields.
///
/// Subclasses must override `pair_cutoff()` to report the interaction
/// cutoff in Ångströms; the base wires it into the `NeighborConsumer`
/// protocol for every pair-cutoff evaluator.
class PairEvaluator : public cpu::ComponentBase, public cpu::NeighborConsumer {
  public:
    /// @brief Pair interaction cutoff [Å]. Used by category drivers
    /// to size MPI halo exchange and surfaced in error messages when
    /// a `NeighborList` is missing.
    virtual Real pair_cutoff() const = 0;

    /// @brief Pair-cutoff evaluators are always neighbor consumers — the
    /// driver scan probes this typed accessor (missing-NL gate + halo
    /// sizing), no `dynamic_cast`.
    const NeighborConsumer* as_neighbor_consumer() const override {
        return this;
    }

    /// @brief The neighbor list must cover the pair interaction cutoff.
    /// @return `pair_cutoff()` [Å].
    Real neighbor_cutoff() const override { return pair_cutoff(); }
};

} // namespace cpu
} // namespace ATV
